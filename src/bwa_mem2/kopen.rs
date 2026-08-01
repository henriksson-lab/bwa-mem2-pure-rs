#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kopen.cpp`.

use std::fs::File;
use std::io::Read;
use std::process::Child;

#[cfg(unix)]
use std::os::fd::AsRawFd;

// --- kopen.cpp ---

const KO_STDIN: i32 = 1;
const KO_FILE: i32 = 2;
const KO_PIPE: i32 = 3;
const KO_HTTP: i32 = 4;
const KO_FTP: i32 = 5;

#[doc = "Original struct: ftpaux_t (bwa-mem2/src/kopen.cpp)"]
#[derive(Debug, Default, Clone)]
pub struct ftpaux_t {
    pub fd: i32,
}

#[doc = "Original struct: koaux_t (bwa-mem2/src/kopen.cpp)"]
pub struct koaux_t {
    pub type_: i32,
    pub fd: i32,
    pub pid: i32,
    pub file: Option<File>,
    pub child: Option<Child>,
    pub pipe: Option<Box<dyn Read + Send>>,
}

impl std::fmt::Debug for koaux_t {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("koaux_t")
            .field("type_", &self.type_)
            .field("fd", &self.fd)
            .field("pid", &self.pid)
            .field("has_file", &self.file.is_some())
            .field("has_child", &self.child.is_some())
            .field("has_pipe", &self.pipe.is_some())
            .finish()
    }
}

impl Default for koaux_t {
    fn default() -> Self {
        Self {
            type_: 0,
            fd: -1,
            pid: 0,
            file: None,
            child: None,
            pipe: None,
        }
    }
}

fn open_remote_reader(url: &str) -> Option<Box<dyn Read + Send>> {
    let response = ureq::get(url)
        .header("User-Agent", "bwa-mem2-rs")
        .call()
        .ok()?;
    Some(Box::new(response.into_body().into_reader()))
}

fn remote_aux(url: &str, type_: i32) -> Option<koaux_t> {
    let pipe = open_remote_reader(url)?;
    Some(koaux_t {
        type_,
        fd: -1,
        pid: 0,
        file: None,
        child: None,
        pipe: Some(pipe),
    })
}

#[doc = "Original function: socket_wait:65"]
pub(crate) fn socket_wait(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("socket_wait")
}

#[doc = "Original function: socket_connect:80"]
pub(crate) fn socket_connect(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("socket_connect")
}

#[doc = "Original function: write_bytes:102"]
pub(crate) fn write_bytes(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("write_bytes")
}

#[doc = "Original function: http_open:117"]
pub(crate) fn http_open(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("http_open")
}

#[doc = "Original function: kftp_get_response:191"]
pub(crate) fn kftp_get_response(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("kftp_get_response")
}

#[doc = "Original function: kftp_send_cmd:215"]
pub(crate) fn kftp_send_cmd(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("kftp_send_cmd")
}

#[doc = "Original function: ftp_open:222"]
pub(crate) fn ftp_open(_arg0: crate::support::Opaque) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("ftp_open")
}

/// Split a shell-style command string into an argv vector.
///
/// Returns `None` if `cmd` is empty after trimming; otherwise splits on
/// whitespace. Used by the pipe-open path of `kopen` (currently inert
/// in the Rust scaffold).
#[doc = "Original function: cmd2argv:278"]
pub fn cmd2argv(cmd: &str) -> Option<Vec<String>> {
    let trimmed = cmd.trim();
    if trimmed.is_empty() {
        return None;
    }
    Some(
        trimmed
            .split_whitespace()
            .map(ToString::to_string)
            .collect(),
    )
}

/// Open a file/stdin/http/ftp/pipe spec and return a `koaux_t` handle.
///
/// Mirrors C's `void *kopen(const char *fn, int *_fd)`. `"-"` selects
/// stdin; `"http://"` / `"https://"` are streamed with a Rust HTTP client.
/// `"<"`-prefixed shell command specs are intentionally unsupported; pipe
/// command output to stdin and pass `"-"` instead.
///
/// # Arguments
/// * `fn_` - filename or URL to open
/// * `fd_out` - out parameter for the raw fd, `-1` if unknown
#[doc = "Original function: kopen:312"]
pub fn kopen(fn_: &str, fd_out: &mut i32) -> Option<koaux_t> {
    *fd_out = -1;
    if fn_.starts_with("http://") || fn_.starts_with("https://") {
        let aux = remote_aux(fn_, KO_HTTP)?;
        *fd_out = aux.fd;
        return Some(aux);
    }
    if fn_.starts_with("ftp://") {
        eprintln!("[E::kopen] ftp:// inputs are not supported");
        return None;
    }
    if fn_ == "-" {
        *fd_out = 0;
        return Some(koaux_t {
            type_: KO_STDIN,
            fd: 0,
            pid: 0,
            file: None,
            child: None,
            pipe: None,
        });
    }

    let trimmed = fn_.trim_start();
    if trimmed.starts_with('<') {
        eprintln!(
            "[E::kopen] shell command inputs are not supported; pipe command output to stdin and use '-'"
        );
        return None;
    }

    let file = File::open(fn_).ok()?;
    #[cfg(unix)]
    let fd = file.as_raw_fd();
    #[cfg(not(unix))]
    let fd = -1;
    *fd_out = fd;
    Some(koaux_t {
        type_: KO_FILE,
        fd,
        pid: 0,
        file: Some(file),
        child: None,
        pipe: None,
    })
}

/// Release a handle returned by `kopen`.
///
/// For `KO_FILE` aux objects the inner `File` is dropped, closing the
/// fd via the standard library rather than leaking it as the upstream
/// C code does. Always returns 0.
#[doc = "Original function: kclose:386"]
pub fn kclose(mut aux: koaux_t) -> i32 {
    if let Some(file) = aux.file.take() {
        drop(file);
    }
    if let Some(pipe) = aux.pipe.take() {
        drop(pipe);
    }
    if let Some(mut child) = aux.child.take() {
        let _ = child.wait();
    }
    0
}

#[doc = "Original function: main:401"]
pub(crate) fn main(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("main")
}

#[cfg(test)]
mod tests {
    use super::{cmd2argv, kclose, kopen, KO_FILE, KO_STDIN};
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before epoch")
            .as_nanos();
        p.push(format!("bwa_mem2_rs_kopen_{name}_{nanos}.txt"));
        p
    }

    #[test]
    fn cmd2argv_splits_trimmed_whitespace_separated_words() {
        let argv = cmd2argv("  bwa-mem2   mem  ref.fa  ").expect("argv");
        assert_eq!(argv, vec!["bwa-mem2", "mem", "ref.fa"]);
        assert!(cmd2argv("   ").is_none());
    }

    #[test]
    fn kopen_opens_plain_file_and_reports_fd() {
        let path = temp_path("plain");
        fs::write(&path, b"hello").expect("write file");
        let mut fd = -1;
        let aux = kopen(path.to_str().expect("utf8"), &mut fd).expect("kopen");
        assert_eq!(aux.type_, KO_FILE);
        #[cfg(unix)]
        assert!(fd >= 0);
        #[cfg(not(unix))]
        assert_eq!(fd, -1);
        assert_eq!(aux.fd, fd);
        assert_eq!(kclose(aux), 0);
        fs::remove_file(&path).expect("cleanup");
    }

    #[test]
    fn kopen_dash_returns_stdin_handle() {
        let mut fd = -1;
        let aux = kopen("-", &mut fd).expect("stdin aux");
        assert_eq!(aux.type_, KO_STDIN);
        assert_eq!(fd, 0);
        assert_eq!(kclose(aux), 0);
    }

    #[test]
    fn kopen_rejects_shell_command_inputs() {
        let mut fd = -1;
        assert!(kopen("< printf ACGT", &mut fd).is_none());
        assert_eq!(fd, -1);
    }
}
