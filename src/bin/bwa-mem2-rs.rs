fn main() {
    let argv: Vec<String> = std::env::args().collect();

    // Limit glibc malloc arenas. Default is `8 * ncpu` (320 on a 40-core box) — under memory
    // pressure this causes __malloc_trim to consume 25%+ of CPU as glibc walks fragmented arenas.
    // With our compute-bound workload, we only need one arena per worker thread. Parse `-t N`
    // (or `--threads N`) from argv; fall back to 4 if not specified.
    //
    // Also disable trim by raising M_TRIM_THRESHOLD to its maximum (effectively never returns
    // memory to the OS via sbrk, just keeps it for reuse). Trim shows up at ~1% on long runs
    // (700K reads) even after the arena fix; with a long-lived process there's no benefit to
    // returning memory to the OS, since we'll just allocate again.
    #[cfg(target_os = "linux")]
    unsafe {
        let n_threads = parse_thread_count(&argv).unwrap_or(4).max(1);
        // Skip if user already set MALLOC_ARENA_MAX so they can override.
        if std::env::var_os("MALLOC_ARENA_MAX").is_none() {
            libc::mallopt(libc::M_ARENA_MAX, n_threads as libc::c_int);
        }
        if std::env::var_os("MALLOC_TRIM_THRESHOLD_").is_none() {
            libc::mallopt(libc::M_TRIM_THRESHOLD, libc::c_int::MAX);
        }
    }

    std::process::exit(bwa_mem2_pure_rs::generated::main_cpp::main(&argv));
}

#[cfg(target_os = "linux")]
fn parse_thread_count(argv: &[String]) -> Option<usize> {
    let mut iter = argv.iter().peekable();
    while let Some(a) = iter.next() {
        match a.as_str() {
            "-t" | "--threads" => {
                if let Some(v) = iter.next() {
                    return v.parse::<usize>().ok();
                }
            }
            s if s.starts_with("-t") && s.len() > 2 => {
                return s[2..].parse::<usize>().ok();
            }
            _ => {}
        }
    }
    None
}
