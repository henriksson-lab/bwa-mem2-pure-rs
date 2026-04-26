#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "bwa-mem2" / "src"
GENERATED_ROOT = REPO_ROOT / "src" / "generated"
CCC_BIN = Path("/data/henriksson/github/claude/code-complexity-comparator/target/debug/ccc-rs")
ORDER_CSV = REPO_ROOT / "ccc_order.csv"
MAPPING_TOML = REPO_ROOT / "ccc_mapping.toml"


CPP_EXTENSIONS = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}
RUST_KEYWORDS = {
    "as", "break", "const", "continue", "crate", "else", "enum", "extern", "false", "fn",
    "for", "if", "impl", "in", "let", "loop", "match", "mod", "move", "mut", "pub",
    "ref", "return", "self", "Self", "static", "struct", "super", "trait", "true", "type",
    "unsafe", "use", "where", "while", "async", "await", "dyn", "abstract", "become",
    "box", "do", "final", "macro", "override", "priv", "typeof", "unsized", "virtual",
    "yield", "try",
}


@dataclass(frozen=True)
class StructDef:
    name: str
    file_path: Path
    kind: str


@dataclass
class FunctionDef:
    other_name: str
    file_path: Path
    line_start: int
    rust_name: str
    rust_class: str | None
    arg_count: int
    returns_value: bool
    returns_self: bool


def run_ccc_analyze(path: Path) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".json") as tmp:
        subprocess.run(
            [str(CCC_BIN), "analyze", str(path), "-l", "cpp", "-o", tmp.name],
            check=True,
            cwd=REPO_ROOT,
            stdout=subprocess.DEVNULL,
        )
        return json.loads(Path(tmp.name).read_text())


def sanitize_module_name(path: Path) -> str:
    stem = path.name.replace(".", "_").replace("-", "_")
    return re.sub(r"[^0-9A-Za-z_]", "_", stem).lower()


def sanitize_type_name(name: str) -> str:
    value = re.sub(r"[^0-9A-Za-z_]", "_", name)
    if not value:
        value = "unnamed_type"
    if value[0].isdigit():
        value = f"_{value}"
    return value


def sanitize_fn_name(name: str) -> str:
    value = re.sub(r"[^0-9A-Za-z_]", "_", name)
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        value = "stub"
    if value[0].isdigit():
        value = f"_{value}"
    if value in RUST_KEYWORDS:
        value = f"r#{value}"
    return value


def rust_path_for(source_path: Path) -> str:
    return f"src/generated/{sanitize_module_name(source_path)}.rs"


def supplemental_structs(path: Path) -> list[StructDef]:
    text = path.read_text(errors="ignore")
    found: dict[str, StructDef] = {}
    typedef_pattern = re.compile(r"typedef\s+(struct|union)\s*(?:\w+\s*)?\{.*?\}\s*(\w+)\s*;", re.S)
    named_pattern = re.compile(r"\b(struct|class|union)\s+(\w+)\s*(?:[:\s\n\r\t]*\{)", re.S)
    for kind, name in typedef_pattern.findall(text):
        found[name] = StructDef(name=name, file_path=path, kind=kind)
    for kind, name in named_pattern.findall(text):
        found.setdefault(name, StructDef(name=name, file_path=path, kind=kind))
    return sorted(found.values(), key=lambda item: item.name)


def map_function(report_fn: dict) -> FunctionDef:
    other_name = report_fn["name"]
    parts = other_name.split("::")
    rust_class = None
    raw_name = parts[-1]
    if len(parts) > 1:
        rust_class = sanitize_type_name(parts[-2])
        class_name = parts[-2]
        if raw_name == class_name:
            rust_name = "ctor"
            returns_self = True
        elif raw_name == f"~{class_name}":
            rust_name = "dtor"
            returns_self = False
        else:
            rust_name = sanitize_fn_name(raw_name)
            returns_self = False
    else:
        rust_name = sanitize_fn_name(raw_name)
        returns_self = False
    outputs = report_fn.get("signature", {}).get("outputs", [])
    returns_value = returns_self or bool(outputs and outputs[0].get("text", "").strip() != "void")
    arg_count = len(report_fn.get("signature", {}).get("inputs", []))
    return FunctionDef(
        other_name=other_name,
        file_path=Path(report_fn["location"]["file"]),
        line_start=report_fn["location"]["line_start"],
        rust_name=rust_name,
        rust_class=rust_class,
        arg_count=arg_count,
        returns_value=returns_value,
        returns_self=returns_self,
    )


def format_struct(item: StructDef) -> str:
    return (
        f'#[doc = "Original {item.kind}: {item.name} ({item.file_path.relative_to(REPO_ROOT).as_posix()})"]\n'
        "#[derive(Debug, Default, Clone, Copy)]\n"
        f"pub struct {sanitize_type_name(item.name)} {{\n"
        "    pub _opaque: (),\n"
        "}\n"
    )


def format_args(count: int) -> str:
    return ", ".join(f"_arg{i}: crate::support::Opaque" for i in range(count))


def format_stub_call(func: FunctionDef, return_ty: str) -> str:
    quoted = json.dumps(func.other_name)
    return f"crate::support::stub::<{return_ty}>({quoted})"


def format_function(func: FunctionDef, indent: str = "") -> str:
    args = format_args(func.arg_count)
    if func.returns_self:
        return_ty = "Self"
    elif func.returns_value:
        return_ty = "crate::support::Opaque"
    else:
        return_ty = "()"
    signature = f"pub fn {func.rust_name}({args}) -> {return_ty}" if return_ty != "()" else f"pub fn {func.rust_name}({args})"
    lines = [
        f'{indent}#[doc = "Original function: {func.other_name}:{func.line_start}"]',
        f"{indent}{signature} {{",
        f"{indent}    {format_stub_call(func, return_ty)}",
        f"{indent}}}",
    ]
    return "\n".join(lines)


def write_module(path: Path, structs: list[StructDef], functions: list[FunctionDef]) -> None:
    lines = [
        "#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]",
        "",
        f'//! Generated scaffold for `{path.relative_to(REPO_ROOT).as_posix()}`.',
        "",
    ]
    emitted_types: set[str] = set()
    all_types = sorted(structs, key=lambda item: (item.name, item.file_path.as_posix()))
    class_names = sorted({fn.rust_class for fn in functions if fn.rust_class})
    for class_name in class_names:
        if class_name not in {sanitize_type_name(item.name) for item in all_types}:
            all_types.append(StructDef(name=class_name, file_path=path, kind="class"))
    all_types.sort(key=lambda item: sanitize_type_name(item.name))
    for item in all_types:
        rust_name = sanitize_type_name(item.name)
        if rust_name in emitted_types:
            continue
        emitted_types.add(rust_name)
        lines.append(format_struct(item))
    free_functions = [fn for fn in functions if fn.rust_class is None]
    methods_by_class: dict[str, list[FunctionDef]] = defaultdict(list)
    for func in functions:
        if func.rust_class is not None:
            methods_by_class[func.rust_class].append(func)
    for func in free_functions:
        lines.append(format_function(func))
        lines.append("")
    for class_name in sorted(methods_by_class):
        lines.append(f"impl {class_name} {{")
        for func in sorted(methods_by_class[class_name], key=lambda item: (item.line_start, item.other_name)):
            lines.append(format_function(func, indent="    "))
            lines.append("")
        if lines[-1] == "":
            lines.pop()
        lines.append("}")
        lines.append("")
    content = "\n".join(lines).rstrip() + "\n"
    (GENERATED_ROOT / f"{sanitize_module_name(path)}.rs").write_text(content)


def write_generated_mod(files: list[Path]) -> None:
    lines = ["#![allow(dead_code)]", ""]
    for path in files:
        lines.append(f"pub mod {sanitize_module_name(path)};")
    lines.append("")
    (GENERATED_ROOT / "mod.rs").write_text("\n".join(lines))


def write_mapping(functions: list[FunctionDef]) -> None:
    lines: list[str] = []
    for func in sorted(functions, key=lambda item: (item.file_path.as_posix(), item.line_start, item.other_name)):
        lines.append("[[entries]]")
        lines.append(f'rust = "{func.rust_name}"')
        lines.append(f'rust_path = "{rust_path_for(func.file_path)}"')
        lines.append(f"rust_line = {func.line_start}")
        if func.rust_class is not None:
            lines.append(f'rust_class = "{func.rust_class}"')
        lines.append(f'other = "{func.other_name}"')
        lines.append(f'other_path = "{func.file_path.relative_to(REPO_ROOT).as_posix()}"')
        lines.append(f"other_line = {func.line_start}")
        lines.append("")
    MAPPING_TOML.write_text("\n".join(lines).rstrip() + "\n")


def write_order_csv() -> None:
    subprocess.run(
        [
            str(CCC_BIN),
            "order",
            str(SRC_ROOT),
            "--recurse",
            "-l",
            "cpp",
            "-o",
            str(ORDER_CSV),
        ],
        check=True,
        cwd=REPO_ROOT,
    )


def main() -> None:
    if not CCC_BIN.exists():
        raise SystemExit(f"missing ccc binary: {CCC_BIN}")
    GENERATED_ROOT.mkdir(parents=True, exist_ok=True)
    source_files = sorted(path for path in SRC_ROOT.iterdir() if path.suffix in CPP_EXTENSIONS)
    structs_by_file: dict[Path, list[StructDef]] = defaultdict(list)
    functions_by_file: dict[Path, list[FunctionDef]] = defaultdict(list)
    all_functions: list[FunctionDef] = []
    for path in source_files:
        report = run_ccc_analyze(path)
        seen_structs: set[str] = set()
        for entry in report.get("structs", []):
            name = entry["name"]
            if name in seen_structs:
                continue
            seen_structs.add(name)
            structs_by_file[path].append(StructDef(name=name, file_path=path, kind=entry.get("kind", "struct")))
        for entry in supplemental_structs(path):
            if entry.name not in seen_structs:
                seen_structs.add(entry.name)
                structs_by_file[path].append(entry)
        for entry in report.get("functions", []):
            func = map_function(entry)
            functions_by_file[path].append(func)
            all_functions.append(func)
    for path, items in functions_by_file.items():
        seen: dict[tuple[str | None, str], list[FunctionDef]] = defaultdict(list)
        for func in items:
            seen[(func.rust_class, func.rust_name)].append(func)
        for duplicates in seen.values():
            if len(duplicates) < 2:
                continue
            for func in duplicates:
                func.rust_name = f"{func.rust_name}__L{func.line_start}"
    for path in source_files:
        write_module(path, structs_by_file[path], functions_by_file[path])
    write_generated_mod(source_files)
    write_mapping(all_functions)
    write_order_csv()
    print(
        f"generated {len(source_files)} modules, "
        f"{sum(len(items) for items in structs_by_file.values())} structs, "
        f"{len(all_functions)} functions"
    )


if __name__ == "__main__":
    main()
