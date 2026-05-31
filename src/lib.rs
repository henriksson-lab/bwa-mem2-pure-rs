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
}
