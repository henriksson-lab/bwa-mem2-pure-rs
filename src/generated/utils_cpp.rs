#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/utils.cpp`.

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fmt::Arguments;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug)]
pub enum ErrFile {
    File(File),
    Stdin,
    Stdout,
}

impl ErrFile {
    fn open_path(path: &str, mode: &str) -> std::io::Result<Self> {
        let mut opts = OpenOptions::new();
        if mode.contains('r') {
            opts.read(true);
        }
        if mode.contains('w') {
            opts.write(true).create(true).truncate(true);
        } else if mode.contains('a') {
            opts.write(true).create(true).append(true);
        }
        if mode.contains('+') {
            opts.read(true).write(true);
        }
        opts.open(path).map(ErrFile::File)
    }

    fn write_all_checked(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match self {
            ErrFile::File(file) => file.write_all(buf),
            ErrFile::Stdout => std::io::stdout().write_all(buf),
            ErrFile::Stdin => Err(std::io::Error::new(std::io::ErrorKind::BrokenPipe, "stdin is not writable")),
        }
    }

    fn flush_checked(&mut self) -> std::io::Result<()> {
        match self {
            ErrFile::File(file) => file.flush(),
            ErrFile::Stdout => std::io::stdout().flush(),
            ErrFile::Stdin => Ok(()),
        }
    }

    fn read_exact_checked(&mut self, buf: &mut [u8]) -> std::io::Result<()> {
        match self {
            ErrFile::File(file) => file.read_exact(buf),
            ErrFile::Stdin => std::io::stdin().read_exact(buf),
            ErrFile::Stdout => Err(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "stdout is not readable")),
        }
    }

    fn read_byte_checked(&mut self) -> std::io::Result<Option<u8>> {
        let mut byte = [0_u8; 1];
        let n = match self {
            ErrFile::File(file) => file.read(&mut byte)?,
            ErrFile::Stdin => std::io::stdin().read(&mut byte)?,
            ErrFile::Stdout => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    "stdout is not readable",
                ))
            }
        };
        if n == 0 { Ok(None) } else { Ok(Some(byte[0])) }
    }

    fn seek_checked(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            ErrFile::File(file) => file.seek(pos),
            ErrFile::Stdin | ErrFile::Stdout => Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "stream does not support seek",
            )),
        }
    }

    fn sync_if_regular_file(&self) -> std::io::Result<()> {
        if let ErrFile::File(file) = self {
            let meta = file.metadata()?;
            if meta.is_file() {
                file.sync_all()?;
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum ErrGzFile {
    ReaderFile(MultiGzDecoder<BufReader<File>>),
    ReaderStdin(MultiGzDecoder<BufReader<std::io::Stdin>>),
    WriterFile(GzEncoder<File>),
    WriterStdout(GzEncoder<std::io::Stdout>),
}

impl ErrGzFile {
    fn open_path(path: &str, mode: &str) -> std::io::Result<Self> {
        if mode.contains('r') {
            let file = File::open(path)?;
            Ok(ErrGzFile::ReaderFile(MultiGzDecoder::new(BufReader::new(file))))
        } else {
            let file = OpenOptions::new().write(true).create(true).truncate(true).open(path)?;
            Ok(ErrGzFile::WriterFile(GzEncoder::new(file, Compression::default())))
        }
    }
}

fn format_error_line(header: &str, message: &str, abort_suffix: bool) -> String {
    let suffix = if abort_suffix { " Abort!" } else { "" };
    format!("[{}] {}{}", header, message, suffix)
}

#[doc = "Original function: err_xopen_core:56"]
pub fn err_xopen_core(func: &str, fn_: &str, mode: &str) -> ErrFile {
    if fn_ == "-" {
        return if mode.contains('r') { ErrFile::Stdin } else { ErrFile::Stdout };
    }
    match ErrFile::open_path(fn_, mode) {
        Ok(fp) => fp,
        Err(err) => err_fatal(func, format_args!("fail to open file '{}' : {}", fn_, err)),
    }
}

#[doc = "Original function: err_xreopen_core:67"]
pub fn err_xreopen_core(func: &str, fn_: &str, mode: &str, _fp: ErrFile) -> ErrFile {
    match ErrFile::open_path(fn_, mode) {
        Ok(fp) => fp,
        Err(err) => err_fatal(func, format_args!("fail to open file '{}' : {}", fn_, err)),
    }
}

#[doc = "Original function: err_xzopen_core:75"]
pub fn err_xzopen_core(func: &str, fn_: &str, mode: &str) -> ErrGzFile {
    if fn_ == "-" {
        if mode.contains('r') {
            return ErrGzFile::ReaderStdin(MultiGzDecoder::new(BufReader::new(std::io::stdin())));
        }
        return ErrGzFile::WriterStdout(GzEncoder::new(std::io::stdout(), Compression::default()));
    }
    match ErrGzFile::open_path(fn_, mode) {
        Ok(fp) => fp,
        Err(err) => err_fatal(func, format_args!("fail to open file '{}' : {}", fn_, err)),
    }
}

#[doc = "Original function: err_fatal:90"]
pub fn err_fatal(header: &str, args: Arguments<'_>) -> ! {
    eprintln!("{}", format_error_line(header, &args.to_string(), false));
    std::process::exit(libc::EXIT_FAILURE);
}

#[doc = "Original function: err_fatal_core:101"]
pub fn err_fatal_core(header: &str, args: Arguments<'_>) -> ! {
    eprintln!("{}", format_error_line(header, &args.to_string(), true));
    std::process::abort();
}

#[doc = "Original function: _err_fatal_simple:112"]
pub fn err_fatal_simple(func: &str, msg: &str) -> ! {
    eprintln!("{}", format_error_line(func, msg, false));
    std::process::exit(libc::EXIT_FAILURE);
}

#[doc = "Original function: _err_fatal_simple_core:118"]
pub fn err_fatal_simple_core(func: &str, msg: &str) -> ! {
    eprintln!("{}", format_error_line(func, msg, true));
    std::process::abort();
}

#[doc = "Original function: err_fwrite:124"]
pub fn err_fwrite(ptr: &[u8], size: usize, nmemb: usize, stream: &mut ErrFile) -> usize {
    let want = size.saturating_mul(nmemb);
    let slice = &ptr[..want];
    if let Err(err) = stream.write_all_checked(slice) {
        err_fatal_simple("fwrite", &err.to_string());
    }
    nmemb
}

#[doc = "Original function: err_fread_noeof:132"]
pub fn err_fread_noeof(ptr: &mut [u8], size: usize, nmemb: usize, stream: &mut ErrFile) -> usize {
    let want = size.saturating_mul(nmemb);
    let slice = &mut ptr[..want];
    if let Err(err) = stream.read_exact_checked(slice) {
        let msg = if err.kind() == std::io::ErrorKind::UnexpectedEof {
            "Unexpected end of file".to_string()
        } else {
            err.to_string()
        };
        err_fatal_simple("fread", &msg);
    }
    nmemb
}

#[doc = "Original function: err_gzread:142"]
pub fn err_gzread(file: &mut ErrGzFile, ptr: &mut [u8], len: u32) -> i32 {
    let want = usize::try_from(len).expect("len");
    let slice = &mut ptr[..want];
    let ret = match file {
        ErrGzFile::ReaderFile(reader) => reader.read(slice),
        ErrGzFile::ReaderStdin(reader) => reader.read(slice),
        ErrGzFile::WriterFile(_) | ErrGzFile::WriterStdout(_) => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "gz stream is not readable",
        )),
    };
    match ret {
        Ok(n) => i32::try_from(n).expect("gzread length overflow"),
        Err(err) => err_fatal_simple("gzread", &err.to_string()),
    }
}

#[doc = "Original function: err_fseek:156"]
pub fn err_fseek(stream: &mut ErrFile, offset: i64, whence: i32) -> i32 {
    let pos = match whence {
        x if x == libc::SEEK_SET => SeekFrom::Start(u64::try_from(offset).expect("offset")),
        x if x == libc::SEEK_CUR => SeekFrom::Current(offset),
        x if x == libc::SEEK_END => SeekFrom::End(offset),
        _ => err_fatal_simple("fseek", "invalid seek mode"),
    };
    if let Err(err) = stream.seek_checked(pos) {
        err_fatal_simple("fseek", &err.to_string());
    }
    0
}

#[doc = "Original function: err_ftell:166"]
pub fn err_ftell(stream: &mut ErrFile) -> i64 {
    match stream.seek_checked(SeekFrom::Current(0)) {
        Ok(pos) => i64::try_from(pos).expect("position overflow"),
        Err(err) => err_fatal_simple("ftell", &err.to_string()),
    }
}

#[doc = "Original function: err_printf:176"]
pub fn err_printf(args: Arguments<'_>) -> i32 {
    let rendered = args.to_string();
    match write!(std::io::stdout(), "{}", rendered) {
        Ok(()) => i32::try_from(rendered.len()).expect("printf length overflow"),
        Err(err) => err_fatal_simple("vfprintf(stdout)", &err.to_string()),
    }
}

#[doc = "Original function: err_fprintf:188"]
pub fn err_fprintf(stream: &mut ErrFile, args: Arguments<'_>) -> i32 {
    let rendered = args.to_string();
    if let Err(err) = stream.write_all_checked(rendered.as_bytes()) {
        err_fatal_simple("vfprintf", &err.to_string());
    }
    i32::try_from(rendered.len()).expect("fprintf length overflow")
}

#[doc = "Original function: err_fputc:200"]
pub fn err_fputc(c: i32, stream: &mut ErrFile) -> i32 {
    let byte = [u8::try_from(c).expect("c")];
    if let Err(err) = stream.write_all_checked(&byte) {
        err_fatal_simple("fputc", &err.to_string());
    }
    c
}

#[doc = "Original function: err_fgets:211"]
pub fn err_fgets(buf: &mut String, size: i32, stream: &mut ErrFile) -> String {
    buf.clear();
    let want = usize::try_from(size.max(0)).expect("size");
    while buf.len() + 1 < want {
        match stream.read_byte_checked() {
            Ok(Some(byte)) => {
                buf.push(char::from(byte));
                if byte == b'\n' {
                    break;
                }
            }
            Ok(None) => err_fatal_simple("fgets", "Unexpected end of file"),
            Err(err) => err_fatal_simple("fgets", &err.to_string()),
        }
    }
    buf.clone()
}

#[doc = "Original function: err_fputs:222"]
pub fn err_fputs(s: &str, stream: &mut ErrFile) -> i32 {
    if let Err(err) = stream.write_all_checked(s.as_bytes()) {
        err_fatal_simple("fputs", &err.to_string());
    }
    0
}

#[doc = "Original function: err_puts:233"]
pub fn err_puts(s: &str) -> i32 {
    let line = format!("{s}\n");
    match write!(std::io::stdout(), "{}", line) {
        Ok(()) => 0,
        Err(err) => err_fatal_simple("puts", &err.to_string()),
    }
}

#[doc = "Original function: err_fflush:244"]
pub fn err_fflush(stream: &mut ErrFile) -> i32 {
    if let Err(err) = stream.flush_checked() {
        err_fatal_simple("fflush", &err.to_string());
    }
    if let Err(err) = stream.sync_if_regular_file() {
        err_fatal_simple("fsync", &err.to_string());
    }
    0
}

#[doc = "Original function: err_fclose:271"]
pub fn err_fclose(mut stream: ErrFile) -> i32 {
    if let Err(err) = stream.flush_checked() {
        err_fatal_simple("fclose", &err.to_string());
    }
    0
}

#[doc = "Original function: err_gzclose:278"]
pub fn err_gzclose(file: ErrGzFile) -> i32 {
    let result = match file {
        ErrGzFile::ReaderFile(_) | ErrGzFile::ReaderStdin(_) => Ok(()),
        ErrGzFile::WriterFile(writer) => writer.finish().map(|_| ()),
        ErrGzFile::WriterStdout(writer) => writer.finish().map(|_| ()),
    };
    if let Err(err) = result {
        err_fatal_simple("gzclose", &err.to_string());
    }
    0
}

#[doc = "Original function: cputime:293"]
pub fn cputime() -> f64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if rc != 0 {
        err_fatal_simple("getrusage", &std::io::Error::last_os_error().to_string());
    }
    let usage = unsafe { usage.assume_init() };
    let user = usage.ru_utime.tv_sec as f64 + usage.ru_utime.tv_usec as f64 * 1e-6;
    let system = usage.ru_stime.tv_sec as f64 + usage.ru_stime.tv_usec as f64 * 1e-6;
    user + system
}

#[doc = "Original function: realtime:300"]
pub fn realtime() -> f64 {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock is before Unix epoch");
    now.as_secs() as f64 + now.subsec_micros() as f64 * 1e-6
}

#[cfg(test)]
mod tests {
    use super::{
        cputime, err_fclose, err_fflush, err_fgets, err_fprintf, err_fputc, err_fread_noeof, err_fseek, err_ftell, err_fwrite,
        err_gzclose, err_gzread, err_xopen_core, err_xreopen_core, err_xzopen_core, format_error_line, realtime, ErrFile,
    };
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::{self, File};
    use std::io::Write;
    use std::path::{Path, PathBuf};

    fn temp_path(label: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        p.push(format!("bwa_mem2_rs_utils_{label}_{stamp}"));
        p
    }

    #[test]
    fn error_line_matches_original_shape() {
        assert_eq!(format_error_line("hdr", "msg", false), "[hdr] msg");
        assert_eq!(format_error_line("hdr", "msg", true), "[hdr] msg Abort!");
    }

    #[test]
    fn cputime_is_non_negative() {
        assert!(cputime() >= 0.0);
    }

    #[test]
    fn realtime_tracks_unix_time() {
        let now = realtime();
        assert!(now > 1_000_000_000.0);
    }

    #[test]
    fn err_xopen_reopen_and_plain_io_round_trip() {
        let path1 = temp_path("plain1");
        let path2 = temp_path("plain2");
        let mut stream = err_xopen_core("test", path1.to_str().expect("utf8"), "w+");
        assert_eq!(err_fwrite(b"hello world", 1, 11, &mut stream), 11);
        assert_eq!(err_fputc(b'\n' as i32, &mut stream), b'\n' as i32);
        assert_eq!(err_fflush(&mut stream), 0);
        assert_eq!(err_fseek(&mut stream, 0, libc::SEEK_SET), 0);
        assert_eq!(err_ftell(&mut stream), 0);
        let mut buf = [0_u8; 11];
        assert_eq!(err_fread_noeof(&mut buf, 1, 11, &mut stream), 11);
        assert_eq!(&buf, b"hello world");
        let mut line = String::new();
        assert_eq!(err_fgets(&mut line, 8, &mut stream), "\n");
        let stream = err_xreopen_core("test", path2.to_str().expect("utf8"), "w", stream);
        assert_eq!(err_fclose(stream), 0);
        assert!(Path::new(&path2).exists());
        let _ = fs::remove_file(path1);
        let _ = fs::remove_file(path2);
    }

    #[test]
    fn err_fprintf_writes_text_to_file() {
        let path = temp_path("fprintf");
        let mut stream = err_xopen_core("test", path.to_str().expect("utf8"), "w");
        assert_eq!(err_fprintf(&mut stream, format_args!("value={}", 42)), 8);
        assert_eq!(err_fflush(&mut stream), 0);
        assert_eq!(err_fclose(stream), 0);
        assert_eq!(fs::read_to_string(&path).expect("read"), "value=42");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn err_xzopen_gzread_and_close_round_trip() {
        let path = temp_path("gzip").with_extension("gz");
        let mut encoder = GzEncoder::new(File::create(&path).expect("create"), Compression::default());
        encoder.write_all(b"gzip-data").expect("write");
        encoder.finish().expect("finish");

        let mut gz = err_xzopen_core("test", path.to_str().expect("utf8"), "r");
        let mut buf = [0_u8; 9];
        assert_eq!(err_gzread(&mut gz, &mut buf, 9), 9);
        assert_eq!(&buf, b"gzip-data");
        assert_eq!(err_gzclose(gz), 0);
        let _ = fs::remove_file(path);
    }

    #[test]
    fn stdin_stdout_handles_are_selected_for_dash() {
        assert!(matches!(err_xopen_core("test", "-", "r"), ErrFile::Stdin));
        assert!(matches!(err_xopen_core("test", "-", "w"), ErrFile::Stdout));
    }
}
