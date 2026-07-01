#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    #[cfg(target_os = "linux")]
    if std::env::var_os("MALLOC_ARENA_MAX").is_none() {
        use std::os::unix::process::CommandExt;

        if let Ok(exe) = std::env::current_exe() {
            let err = std::process::Command::new(exe)
                .args(std::env::args_os().skip(1))
                .env("MALLOC_ARENA_MAX", "1")
                .exec();
            eprintln!("[W::bwa-mem2-rs] failed to re-exec with MALLOC_ARENA_MAX=1: {err}");
        }
    }

    let argv: Vec<String> = std::env::args().collect();

    // The re-exec above sets MALLOC_ARENA_MAX before glibc malloc initializes. That reliably
    // avoids high transient SAM-phase RSS from per-thread arenas; an in-process mallopt() call was
    // not early enough on the PE benchmark. Users can override by setting MALLOC_ARENA_MAX.
    //
    // Also disable trim by raising M_TRIM_THRESHOLD to its maximum (effectively never returns
    // memory to the OS via sbrk, just keeps it for reuse). Trim shows up at ~1% on long runs
    // (700K reads) even after the arena fix; with a long-lived process there's no benefit to
    // returning memory to the OS, since we'll just allocate again.
    #[cfg(target_os = "linux")]
    unsafe {
        if std::env::var_os("MALLOC_TRIM_THRESHOLD_").is_none() {
            libc::mallopt(libc::M_TRIM_THRESHOLD, libc::c_int::MAX);
        }
    }

    std::process::exit(bwa_mem2_pure_rs::bwa_mem2::main::main(&argv));
}
