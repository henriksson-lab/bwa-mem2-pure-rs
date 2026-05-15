#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/runsimd.cpp`.

use std::{env, fs};
use std::path::{Path, PathBuf};


// --- runsimd.cpp ---

#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
#[cfg(unix)]
use std::os::unix::process::CommandExt;

const SIMD_SSE: i32 = 0x1;
const SIMD_SSE2: i32 = 0x2;
const SIMD_SSE3: i32 = 0x4;
const SIMD_SSSE3: i32 = 0x8;
const SIMD_SSE4_1: i32 = 0x10;
const SIMD_SSE4_2: i32 = 0x20;
const SIMD_AVX: i32 = 0x40;
const SIMD_AVX2: i32 = 0x80;
const SIMD_AVX512F: i32 = 0x100;
const SIMD_AVX512BW: i32 = 0x200;

/// Run the x86 `cpuid` instruction at `func_id`/`subfunc_id`.
///
/// Mirrors the upstream `__cpuidex` shim (adapted from linux-sgx's
/// `cpuid_gnu.h`) so the SIMD-dispatch logic can probe CPU features.
/// Writes EAX/EBX/ECX/EDX into `cpuid`. Non-x86 targets receive zeros.
///
/// # Arguments
/// * `cpuid` - 4-element output buffer (EAX, EBX, ECX, EDX)
/// * `func_id` - leaf number passed in EAX
/// * `subfunc_id` - sub-leaf number passed in ECX
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
#[doc = "Original function: __cpuidex:58"]
pub fn cpuidex(cpuid: &mut [i32; 4], func_id: i32, subfunc_id: i32) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        #[cfg(target_arch = "x86")]
        use std::arch::x86::__cpuid_count;
        #[cfg(target_arch = "x86_64")]
        use std::arch::x86_64::__cpuid_count;

        // SAFETY: cpuid is a pure CPU feature query on x86/x86_64.
        let res = unsafe { __cpuid_count(func_id as u32, subfunc_id as u32) };
        cpuid[0] = res.eax as i32;
        cpuid[1] = res.ebx as i32;
        cpuid[2] = res.ecx as i32;
        cpuid[3] = res.edx as i32;
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = (func_id, subfunc_id);
        *cpuid = [0; 4];
    }
}

/// Probe the host's x86 SIMD feature set via cpuid.
///
/// Returns a bitmask combining `SIMD_SSE..SIMD_AVX512BW` flags
/// matching what the host CPU advertises. Returns 0 if cpuid is
/// unsupported (max leaf == 0).
#[doc = "Original function: x86_simd:72"]
pub fn x86_simd() -> i32 {
    let mut flag = 0;
    let mut cpuid = [0_i32; 4];
    cpuidex(&mut cpuid, 0, 0);
    let max_id = cpuid[0];
    if max_id == 0 {
        return 0;
    }
    cpuidex(&mut cpuid, 1, 0);
    if (cpuid[3] >> 25) & 1 != 0 {
        flag |= SIMD_SSE;
    }
    if (cpuid[3] >> 26) & 1 != 0 {
        flag |= SIMD_SSE2;
    }
    if (cpuid[2] >> 0) & 1 != 0 {
        flag |= SIMD_SSE3;
    }
    if (cpuid[2] >> 9) & 1 != 0 {
        flag |= SIMD_SSSE3;
    }
    if (cpuid[2] >> 19) & 1 != 0 {
        flag |= SIMD_SSE4_1;
    }
    if (cpuid[2] >> 20) & 1 != 0 {
        flag |= SIMD_SSE4_2;
    }
    if (cpuid[2] >> 28) & 1 != 0 {
        flag |= SIMD_AVX;
    }
    if max_id >= 7 {
        cpuidex(&mut cpuid, 7, 0);
        if (cpuid[1] >> 5) & 1 != 0 {
            flag |= SIMD_AVX2;
        }
        if (cpuid[1] >> 16) & 1 != 0 {
            flag |= SIMD_AVX512F;
        }
        if (cpuid[1] >> 30) & 1 != 0 {
            flag |= SIMD_AVX512BW;
        }
    }
    flag
}

fn exe_path_impl(exe: &str) -> Result<(String, usize), i32> {
    if exe.is_empty() {
        return Err(-1);
    }
    let last_slash = exe.rfind('/').map(|i| i as isize).unwrap_or(-1);
    let base_st = usize::try_from(last_slash + 1).expect("base_st");
    if exe.starts_with('/') {
        return Ok((exe[..base_st].to_string(), base_st));
    }
    // last_slash >= 0 means the path is relative with a slash (e.g. "dir/exe");
    // resolve against cwd. (Actually, can't be 0 since the bare-name branch handles
    // that case below.)
    if last_slash >= 0 {
        let cwd = env::current_dir().map_err(|_| -1)?;
        let mut buf = cwd.to_string_lossy().to_string();
        buf.push('/');
        buf.push_str(&exe[..base_st]);
        return Ok((buf, base_st));
    }

    let path_env = env::var("PATH").map_err(|_| -2)?;
    for entry in path_env.split(':') {
        let candidate = Path::new(entry).join(exe);
        if let Ok(meta) = fs::metadata(&candidate) {
            #[cfg(unix)]
            {
                if meta.permissions().mode() & 0o100 != 0 {
                    let candidate_str = candidate.to_string_lossy().to_string();
                    return exe_path_impl(&candidate_str);
                }
            }
            #[cfg(not(unix))]
            {
                if meta.is_file() {
                    let candidate_str = candidate.to_string_lossy().to_string();
                    return exe_path_impl(&candidate_str);
                }
            }
        }
    }
    Err(-2) // shouldn't happen!
}

/// Resolve the directory containing executable `exe`.
///
/// Mirrors C's `exe_path`: writes the directory portion (including the
/// trailing slash) into `buf` and the index of the basename into
/// `base_st`. Handles absolute paths, relative paths with slashes
/// (joined to `cwd`), and bare names (searched on `PATH`). Returns 0
/// on success, `-1` for malformed input / cwd failure, `-2` if not
/// found on `PATH`.
///
/// # Arguments
/// * `exe` - executable name or path
/// * `buf` - out buffer to receive the directory prefix
/// * `base_st` - optional out parameter for the basename offset
#[doc = "Original function: exe_path:95"]
pub fn exe_path(exe: &str, _max: i32, buf: &mut String, base_st: Option<&mut i32>) -> i32 {
    match exe_path_impl(exe) {
        Ok((resolved, base)) => {
            *buf = resolved;
            if let Some(slot) = base_st {
                *slot = i32::try_from(base).expect("base");
            }
            0
        }
        Err(code) => code,
    }
}

fn candidate_is_executable(path: &Path) -> bool {
    match fs::metadata(path) {
        Ok(meta) => {
            #[cfg(unix)]
            {
                meta.permissions().mode() & 0o100 != 0
            }
            #[cfg(not(unix))]
            {
                meta.is_file()
            }
        }
        Err(_) => false,
    }
}

/// Try to `exec` the SIMD-suffixed binary at `prefix + simd`.
///
/// Appends `simd` to `prefix`, stats the resulting path, and execs it
/// with the original `argv` if it exists and is executable. On
/// success this never returns. The trailing suffix is stripped from
/// `prefix` before returning so the caller can try the next variant.
///
/// We assume `prefix` is long enough to receive the suffix; on return the original
/// length is restored regardless of success.
///
/// # Arguments
/// * `argv` - full argv of the launcher (forwarded to the child)
/// * `prefix` - in/out path prefix (e.g. `/usr/bin/bwa-mem2`)
/// * `simd` - suffix such as `.avx2`, `.sse42`
#[doc = "Original function: test_and_launch:152"]
pub fn test_and_launch(argv: &[String], prefix: &mut String, simd: &str) {
    let prefix_len = prefix.len();
    prefix.push_str(simd);
    eprintln!(
        "Looking to launch executable \"{}\", simd = {}",
        prefix, simd
    );
    let path = PathBuf::from(prefix.as_str());
    if path.exists() {
        if candidate_is_executable(&path) {
            eprintln!("Launching executable \"{}\"", prefix);
            #[cfg(unix)]
            {
                let err = std::process::Command::new(&path).args(&argv[1..]).exec();
                panic!("exec failed for '{}': {err}", prefix);
            }
            #[cfg(not(unix))]
            {
                let status = std::process::Command::new(&path)
                    .args(&argv[1..])
                    .status()
                    .unwrap_or_else(|e| panic!("spawn failed for '{}': {e}", prefix));
                panic!("child exited with status {status}");
            }
        } else {
            #[cfg(unix)]
            eprintln!(
                "(st.st_mode & S_IXUSR) = {}, can not run executable: {}",
                fs::metadata(&path)
                    .map(|m| m.permissions().mode() & 0o100)
                    .unwrap_or(0),
                prefix
            );
            #[cfg(not(unix))]
            eprintln!("can not run executable: {}", prefix);
        }
    } else {
        eprintln!("stat(prefix, &st) = -1, can not run executable: {}", prefix);
    }
    prefix.truncate(prefix_len);
}

/// Entry point of the runsimd dispatcher shim.
///
/// Resolves the launcher's own directory, probes the CPU's SIMD
/// features via `x86_simd`, and tries to exec the best-matching
/// `bwa-mem2.{avx512bw|avx2|avx|sse42|sse41}` binary. Returns 1 if
/// the path resolution fails and 2 if no SIMD-suffixed binary is
/// found.
#[doc = "Original function: main:176"]
pub fn main(argc: i32, argv: &[String]) -> i32 {
    let _ = argc;
    let Some(argv0) = argv.first() else {
        eprintln!("ERROR: prefix is too long!");
        return 1;
    };
    let mut buf = String::new();
    let mut base_st = 0_i32;
    let ret = exe_path(argv0, 4096, &mut buf, Some(&mut base_st));
    if ret != 0 {
        eprintln!("ERROR: prefix is too long!");
        return 1;
    }
    let mut prefix = buf;
    prefix.push_str(&argv0[usize::try_from(base_st).expect("base_st")..]);
    let simd = x86_simd();
    if simd & SIMD_AVX512BW != 0 {
        test_and_launch(argv, &mut prefix, ".avx512bw");
    }
    if simd & SIMD_AVX2 != 0 {
        test_and_launch(argv, &mut prefix, ".avx2");
    }
    if simd & SIMD_AVX != 0 {
        test_and_launch(argv, &mut prefix, ".avx");
    }
    if simd & SIMD_SSE4_2 != 0 {
        test_and_launch(argv, &mut prefix, ".sse42");
    }
    if simd & SIMD_SSE4_1 != 0 {
        test_and_launch(argv, &mut prefix, ".sse41");
    }
    eprintln!("ERROR: fail to find the right executable");
    2
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let mut dir = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        dir.push(format!("bwa_mem2_rs_runsimd_{name}_{nanos}"));
        fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    #[test]
    fn cpuidex_returns_zero_leaf_data_consistently() {
        let mut cpuid = [0_i32; 4];
        cpuidex(&mut cpuid, 0, 0);
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        assert!(cpuid[0] >= 0);
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        assert_eq!(cpuid, [0; 4]);
    }

    #[test]
    fn exe_path_handles_absolute_and_relative_paths() {
        let dir = temp_dir("path");
        let exe = dir.join("prog");
        fs::write(&exe, b"#!/bin/sh\n").expect("write exe");
        #[cfg(unix)]
        {
            let mut perms = fs::metadata(&exe).expect("meta").permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&exe, perms).expect("chmod");
        }

        let mut buf = String::new();
        let mut base = -1_i32;
        assert_eq!(
            exe_path(exe.to_str().expect("utf8"), 4096, &mut buf, Some(&mut base)),
            0
        );
        assert!(buf.ends_with('/'));
        assert_eq!(
            base,
            1 + i32::try_from(buf.rfind('/').expect("slash")).expect("idx")
        );

        let cwd = env::current_dir().expect("cwd");
        env::set_current_dir(&dir).expect("cd temp");
        let mut rel_buf = String::new();
        assert_eq!(exe_path("./prog", 4096, &mut rel_buf, None), 0);
        assert!(rel_buf.ends_with("./"));
        env::set_current_dir(cwd).expect("restore cwd");
        fs::remove_dir_all(dir).expect("cleanup");
    }

    #[test]
    fn exe_path_searches_path_for_bare_name() {
        let dir = temp_dir("search");
        let exe = dir.join("tool");
        fs::write(&exe, b"#!/bin/sh\n").expect("write exe");
        #[cfg(unix)]
        {
            let mut perms = fs::metadata(&exe).expect("meta").permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&exe, perms).expect("chmod");
        }

        let old_path = env::var("PATH").unwrap_or_default();
        env::set_var("PATH", dir.to_string_lossy().to_string());
        let mut buf = String::new();
        let rc = exe_path("tool", 4096, &mut buf, None);
        env::set_var("PATH", old_path);
        assert_eq!(rc, 0);
        assert!(buf.ends_with('/'));
        fs::remove_dir_all(dir).expect("cleanup");
    }

    #[test]
    fn main_returns_error_when_no_simd_binary_exists() {
        let dir = temp_dir("main");
        let exe = dir.join("bwa-mem2");
        fs::write(&exe, b"#!/bin/sh\n").expect("write exe");
        #[cfg(unix)]
        {
            let mut perms = fs::metadata(&exe).expect("meta").permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&exe, perms).expect("chmod");
        }

        let args = vec![exe.to_string_lossy().to_string()];
        assert_eq!(main(1, &args), 2);
        fs::remove_dir_all(dir).expect("cleanup");
    }
}
