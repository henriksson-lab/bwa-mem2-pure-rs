#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

pub mod bwa_mem2;
pub mod mem_api;
pub mod output;
mod support;

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;

    #[derive(Clone, Debug, Eq, PartialEq)]
    struct OriginalFunctionMapping {
        module_path: String,
        original: String,
        original_line: u32,
        rust: String,
    }

    #[test]
    fn translated_stubs_are_not_public_api() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/bwa_mem2");
        let mut files = Vec::new();
        collect_rs_files(&root, &mut files);

        let mut exposed = Vec::new();
        for path in files {
            let text = fs::read_to_string(&path).expect("read source file");
            let mut current_fn = String::new();
            let mut collecting_fn = false;

            for line in text.lines() {
                let trimmed = line.trim_start();
                if trimmed.starts_with("pub fn ")
                    || trimmed.starts_with("pub(crate) fn ")
                    || trimmed.starts_with("pub(super) fn ")
                {
                    current_fn.clear();
                    current_fn.push_str(trimmed);
                    collecting_fn = !trimmed.contains('{') && !trimmed.ends_with(';');
                } else if collecting_fn {
                    current_fn.push(' ');
                    current_fn.push_str(trimmed);
                    collecting_fn = !trimmed.contains('{') && !trimmed.ends_with(';');
                }

                if trimmed.contains("crate::support::stub")
                    || trimmed.contains("inside #if 0 in upstream")
                {
                    assert!(
                        !current_fn.starts_with("pub fn "),
                        "stub placeholder exposed as public API in {}: {}",
                        path.display(),
                        current_fn
                    );
                    if current_fn.is_empty() {
                        exposed.push(path.display().to_string());
                    }
                }
            }
        }
        assert!(
            exposed.is_empty(),
            "stub without tracked function: {exposed:?}"
        );
    }

    #[test]
    fn placeholder_only_modules_are_not_public_api() {
        let mod_rs =
            fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("src/bwa_mem2/mod.rs"))
                .expect("read bwa_mem2 mod.rs");

        for module in ["khash", "kthread", "memcpy_bwamem", "profiling"] {
            assert!(
                !mod_rs.contains(&format!("pub mod {module};")),
                "placeholder-only module unexpectedly exposed as public API: {module}"
            );
            assert!(
                mod_rs.contains(&format!("pub(crate) mod {module};")),
                "placeholder-only module should remain crate-private: {module}"
            );
        }
    }

    #[test]
    fn translated_original_functions_use_canonical_snake_case() {
        let entries = collect_original_function_mappings();
        for entry in &entries {
            let canonical = canonical_rust_function_name(
                &entry.original,
                entry.original_line,
                duplicate_base_name_count(&entries, &entry.module_path, &entry.original),
            );
            assert_eq!(
                entry.rust, canonical,
                "non-canonical Rust name for original function {}:{} in {}",
                entry.original, entry.original_line, entry.module_path
            );
        }
    }

    fn collect_rs_files(dir: &Path, out: &mut Vec<std::path::PathBuf>) {
        for entry in fs::read_dir(dir).expect("read source dir") {
            let path = entry.expect("read source entry").path();
            if path.is_dir() {
                collect_rs_files(&path, out);
            } else if path.extension().is_some_and(|ext| ext == "rs") {
                out.push(path);
            }
        }
    }

    fn collect_original_function_mappings() -> Vec<OriginalFunctionMapping> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/bwa_mem2");
        let mut files = Vec::new();
        collect_rs_files(&root, &mut files);
        files.sort();

        let mut entries = Vec::new();
        for path in files {
            let text = fs::read_to_string(&path).expect("read source file");
            let lines: Vec<_> = text.lines().collect();
            for (line_idx, line) in lines.iter().enumerate() {
                let Some((original, original_line)) = parse_original_function_doc(line) else {
                    continue;
                };
                let rust = find_next_fn_name(&lines[line_idx + 1..]).unwrap_or_else(|| {
                    panic!(
                        "Original function tag without following fn in {}:{}",
                        path.display(),
                        line_idx + 1
                    )
                });
                let module_path = path
                    .strip_prefix(Path::new(env!("CARGO_MANIFEST_DIR")).join("src"))
                    .expect("strip src prefix")
                    .to_string_lossy()
                    .replace('\\', "/");
                entries.push(OriginalFunctionMapping {
                    module_path,
                    original,
                    original_line,
                    rust,
                });
            }
        }
        entries
    }

    fn parse_original_function_doc(line: &str) -> Option<(String, u32)> {
        let trimmed = line.trim();
        let prefix = "#[doc = \"Original function: ";
        let suffix = "\"]";
        let body = trimmed.strip_prefix(prefix)?.strip_suffix(suffix)?;
        let (name, line) = body.rsplit_once(':')?;
        Some((
            name.to_string(),
            line.parse().expect("original function source line"),
        ))
    }

    fn find_next_fn_name(lines: &[&str]) -> Option<String> {
        for line in lines.iter().take(16) {
            let trimmed = line.trim_start();
            let Some(fn_pos) = trimmed.find("fn ") else {
                continue;
            };
            let after_fn = &trimmed[fn_pos + 3..];
            let end = match (after_fn.find('('), after_fn.find('<')) {
                (Some(paren), Some(generic)) => paren.min(generic),
                (Some(paren), None) => paren,
                (None, Some(generic)) => generic,
                (None, None) => continue,
            };
            if end == 0 {
                continue;
            };
            return Some(after_fn[..end].trim().to_string());
        }
        None
    }

    fn duplicate_base_name_count(
        entries: &[OriginalFunctionMapping],
        module_path: &str,
        original: &str,
    ) -> usize {
        entries
            .iter()
            .filter(|entry| {
                entry.module_path == module_path
                    && canonical_base_rust_function_name(&entry.original)
                        == canonical_base_rust_function_name(original)
            })
            .count()
    }

    fn canonical_rust_function_name(
        original: &str,
        original_line: u32,
        base_count: usize,
    ) -> String {
        let base = canonical_base_rust_function_name(original);
        if base_count > 1 {
            format!("{base}_l{original_line}")
        } else {
            base
        }
    }

    fn canonical_base_rust_function_name(original: &str) -> String {
        let mut member = original.rsplit("::").next().expect("function name");
        if member.starts_with('~') {
            return "dtor".to_string();
        }
        if let Some((class, method)) = original.rsplit_once("::") {
            if class.rsplit("::").next() == Some(method) {
                return "ctor".to_string();
            }
        }

        member = member.trim_start_matches('_');
        let normalized = member
            .replace("SMEMs", "Smems")
            .replace("LMSsort", "LmsSort")
            .replace("LMSpostproc", "LmsPostproc")
            .replace("BWT", "Bwt")
            .replace("SA", "Sa")
            .replace("HT", "Ht");
        to_snake_case(&normalized)
    }

    fn to_snake_case(name: &str) -> String {
        let mut out = String::new();
        let chars: Vec<char> = name.chars().collect();
        for (idx, &ch) in chars.iter().enumerate() {
            let prev = idx.checked_sub(1).and_then(|i| chars.get(i)).copied();
            let next = chars.get(idx + 1).copied();
            if ch == '_' {
                if !out.ends_with('_') {
                    out.push('_');
                }
            } else if ch.is_ascii_uppercase() {
                if !out.is_empty()
                    && !out.ends_with('_')
                    && (prev.is_some_and(|p| p.is_ascii_lowercase() || p.is_ascii_digit())
                        || (prev.is_some_and(|p| p.is_ascii_uppercase())
                            && next.is_some_and(|n| n.is_ascii_lowercase())))
                {
                    out.push('_');
                }
                out.push(ch.to_ascii_lowercase());
            } else {
                out.push(ch);
            }
        }
        out.trim_matches('_').to_string()
    }
}
