#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kseq.h`.

use crate::bwa_mem2::kstring::{ks_resize, kstring_t};
use flate2::read::MultiGzDecoder;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::process::Child;

// --- kseq.h ---

pub const KS_SEP_SPACE: i32 = 0;
pub const KS_SEP_TAB: i32 = 1;
pub const KS_SEP_LINE: i32 = 2;
pub const KS_SEP_MAX: i32 = 2;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct KseqRecord {
    name: String,
    comment: Option<String>,
    seq: String,
    qual: Option<String>,
}

fn trim_cr(line: &str) -> &str {
    line.strip_suffix('\r').unwrap_or(line)
}

fn split_header(header: &str) -> (String, Option<String>) {
    let (name, comment) = split_header_borrow(header);
    (name.to_string(), comment.map(|c| c.to_string()))
}

#[inline]
fn split_header_borrow(header: &str) -> (&str, Option<&str>) {
    let bytes = header.as_bytes();
    let Some(split_at) = bytes.iter().position(|b| b.is_ascii_whitespace()) else {
        return (header, None);
    };
    let name = &header[..split_at];
    let comment = &header[split_at + 1..];
    if comment.is_empty() {
        (name, None)
    } else {
        (name, Some(comment))
    }
}

fn parse_records(text: &str) -> Vec<KseqRecord> {
    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0_usize;
    let mut records = Vec::new();
    while i < lines.len() {
        let line = trim_cr(lines[i]);
        if line.is_empty() {
            i += 1;
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            let (name, comment) = split_header(header);
            i += 1;
            let mut seq = String::new();
            while i < lines.len() {
                let line = trim_cr(lines[i]);
                if line.starts_with('>') || line.starts_with('@') {
                    break;
                }
                if !line.is_empty() {
                    seq.push_str(line);
                }
                i += 1;
            }
            records.push(KseqRecord {
                name,
                comment,
                seq,
                qual: None,
            });
            continue;
        }
        if let Some(header) = line.strip_prefix('@') {
            let (name, comment) = split_header(header);
            i += 1;
            let mut seq = String::new();
            while i < lines.len() {
                let line = trim_cr(lines[i]);
                if line.starts_with('+') {
                    i += 1;
                    break;
                }
                if !line.is_empty() {
                    seq.push_str(line);
                }
                i += 1;
            }
            let mut qual = String::new();
            while i < lines.len() && qual.len() < seq.len() {
                qual.push_str(trim_cr(lines[i]));
                i += 1;
            }
            records.push(KseqRecord {
                name,
                comment,
                seq,
                qual: Some(qual),
            });
            continue;
        }
        i += 1;
    }
    records
}

fn set_kstring(dst: &mut kstring_t, text: &str) {
    ks_resize(dst, text.len() + 1);
    if dst.s.len() < text.len() + 1 {
        dst.s.resize(text.len() + 1, 0);
        dst.m = dst.s.len();
    }
    dst.s[..text.len()].copy_from_slice(text.as_bytes());
    dst.l = text.len();
    dst.s[dst.l] = 0;
}

#[doc = "Original struct: __kstream_t (bwa-mem2/src/kseq.h)"]
pub struct __kstream_t {
    pub begin: usize,
    pub end: usize,
    pub is_eof: i32,
    input: Vec<u8>,
    cursor: usize,
    records: Vec<KseqRecord>,
    record_index: usize,
    reader: Option<Box<dyn BufRead + Send>>,
    child: Option<Child>,
    pending_header: Option<String>,
    line_buf: String,
}

impl fmt::Debug for __kstream_t {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("__kstream_t")
            .field("begin", &self.begin)
            .field("end", &self.end)
            .field("is_eof", &self.is_eof)
            .field("input_len", &self.input.len())
            .field("cursor", &self.cursor)
            .field("records_len", &self.records.len())
            .field("record_index", &self.record_index)
            .field("has_reader", &self.reader.is_some())
            .field("has_child", &self.child.is_some())
            .field("pending_header", &self.pending_header)
            .finish()
    }
}

impl Default for __kstream_t {
    fn default() -> Self {
        Self {
            begin: 0,
            end: 0,
            is_eof: 0,
            input: Vec::new(),
            cursor: 0,
            records: Vec::new(),
            record_index: 0,
            reader: None,
            child: None,
            pending_header: None,
            line_buf: String::new(),
        }
    }
}

impl __kstream_t {
    pub fn from_text(text: &str) -> Self {
        Self {
            begin: 0,
            end: text.len(),
            is_eof: 0,
            input: text.as_bytes().to_vec(),
            cursor: 0,
            records: parse_records(text),
            record_index: 0,
            reader: None,
            child: None,
            pending_header: None,
            line_buf: String::new(),
        }
    }

    pub fn from_reader(reader: Box<dyn BufRead + Send>) -> Self {
        let mut stream = Self::default();
        stream.reader = Some(reader);
        stream
    }

    pub fn from_reader_with_child(reader: Box<dyn BufRead + Send>, child: Child) -> Self {
        let mut stream = Self::default();
        stream.reader = Some(reader);
        stream.child = Some(child);
        stream
    }

    pub fn from_path(path: &str) -> std::io::Result<Self> {
        let mut probe = File::open(path)?;
        let mut magic = [0_u8; 2];
        let n = probe.read(&mut magic)?;
        drop(probe);
        let file = File::open(path)?;
        let gzipped = path.ends_with(".gz") || (n == 2 && magic == [0x1f, 0x8b]);
        if gzipped {
            let gz = MultiGzDecoder::new(BufReader::new(file));
            Ok(Self::from_reader(Box::new(BufReader::new(gz))))
        } else {
            Ok(Self::from_reader(Box::new(BufReader::new(file))))
        }
    }
}

impl Drop for __kstream_t {
    fn drop(&mut self) {
        let _ = self.reader.take();
        if let Some(mut child) = self.child.take() {
            let _ = child.wait();
        }
    }
}

pub type kstream_t = __kstream_t;
pub type __kstring_t = kstring_t;

#[doc = "Original struct: kseq_t (bwa-mem2/src/kseq.h)"]
#[derive(Debug, Default)]
pub struct kseq_t {
    pub name: kstring_t,
    pub comment: kstring_t,
    pub seq: kstring_t,
    pub qual: kstring_t,
    pub last_char: i32,
    pub f: kstream_t,
}

impl kseq_t {
    pub fn from_text(text: &str) -> Self {
        Self {
            f: kstream_t::from_text(text),
            ..Default::default()
        }
    }

    pub fn from_reader(reader: Box<dyn BufRead + Send>) -> Self {
        Self {
            f: kstream_t::from_reader(reader),
            ..Default::default()
        }
    }

    pub fn from_reader_with_child(reader: Box<dyn BufRead + Send>, child: Child) -> Self {
        Self {
            f: kstream_t::from_reader_with_child(reader, child),
            ..Default::default()
        }
    }

    pub fn from_path(path: &str) -> std::io::Result<Self> {
        Ok(Self {
            f: kstream_t::from_path(path)?,
            ..Default::default()
        })
    }
}

fn read_trimmed_line(
    reader: &mut dyn BufRead,
    buf: &mut String,
) -> std::io::Result<Option<String>> {
    buf.clear();
    let n = reader.read_line(buf)?;
    if n == 0 {
        return Ok(None);
    }
    while matches!(buf.as_bytes().last(), Some(b'\n' | b'\r')) {
        buf.pop();
    }
    // Take buf's allocation rather than cloning. Caller owns the line; next read_line() will
    // allocate fresh into buf. Saves the memcpy of the line bytes (~600 bytes per FASTQ record
    // at 4 lines/record × 1.4M records).
    Ok(Some(std::mem::take(buf)))
}

fn next_stream_header(ks: &mut kstream_t) -> Option<String> {
    if let Some(header) = ks.pending_header.take() {
        return Some(header);
    }
    let reader = ks.reader.as_mut()?;
    loop {
        let line = read_trimmed_line(reader.as_mut(), &mut ks.line_buf).ok()??;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') || line.starts_with('@') {
            return Some(line);
        }
    }
}

// Read next header line into ks.line_buf (preserves capacity across calls instead of
// mem::take'ing the buffer). Returns false on EOF / error. On success, ks.line_buf holds the
// header line (without trailing \r\n). Caller can borrow ks.line_buf as &str directly.
fn next_stream_header_in_place(ks: &mut kstream_t) -> bool {
    if let Some(header) = ks.pending_header.take() {
        ks.line_buf.clear();
        ks.line_buf.push_str(&header);
        return true;
    }
    let Some(reader) = ks.reader.as_mut() else {
        return false;
    };
    loop {
        ks.line_buf.clear();
        let n = match reader.read_line(&mut ks.line_buf) {
            Ok(n) => n,
            Err(_) => return false,
        };
        if n == 0 {
            return false;
        }
        while matches!(ks.line_buf.as_bytes().last(), Some(b'\n' | b'\r')) {
            ks.line_buf.pop();
        }
        if ks.line_buf.is_empty() {
            continue;
        }
        if ks.line_buf.starts_with('>') || ks.line_buf.starts_with('@') {
            return true;
        }
    }
}

fn streaming_kseq_read(seq: &mut kseq_t) -> i64 {
    if !next_stream_header_in_place(&mut seq.f) {
        seq.f.is_eof = 1;
        return -1;
    }

    // Determine record type and emit name/comment directly into the kstring buffers,
    // borrowing into seq.f.line_buf — saves the per-record String allocation that
    // `next_stream_header` previously produced via mem::take.
    let is_fastq = {
        let header_line: &str = seq.f.line_buf.as_str();
        let (header, is_fastq) = if let Some(h) = header_line.strip_prefix('>') {
            (h, false)
        } else if let Some(h) = header_line.strip_prefix('@') {
            (h, true)
        } else {
            seq.f.is_eof = 1;
            return -1;
        };
        let (name, comment) = split_header_borrow(header);
        // set_kstring copies bytes out of the &str borrow into seq.name/seq.comment so the
        // borrow into seq.f.line_buf can end before the body reader needs &mut seq.f.
        set_kstring(&mut seq.name, name);
        set_kstring(&mut seq.comment, comment.unwrap_or(""));
        is_fastq
    };
    // Build body and qual directly in the destination kstring buffers — saves a per-record
    // copy round-trip through intermediate Strings + the .clone()/mem::take dance in
    // read_trimmed_line. read_until appends to a Vec<u8> directly without UTF-8 validation
    // (FASTQ payload is ASCII by spec).
    let kseq_t {
        seq: seq_buf,
        qual: qual_buf,
        f: kstream,
        ..
    } = seq;
    seq_buf.s.clear();
    seq_buf.l = 0;
    qual_buf.s.clear();
    qual_buf.l = 0;
    let reader = kstream.reader.as_mut().expect("streaming reader");

    let body_len: usize;
    if !is_fastq {
        loop {
            let start = seq_buf.s.len();
            let n = reader.read_until(b'\n', &mut seq_buf.s).unwrap_or(0);
            if n == 0 {
                kstream.is_eof = 1;
                break;
            }
            // Trim trailing \r\n / \n / \r from this line.
            while seq_buf.s.len() > start && matches!(seq_buf.s.last(), Some(b'\n' | b'\r')) {
                seq_buf.s.pop();
            }
            // Header signals end of body — save line into pending_header (as String) and back out.
            if matches!(seq_buf.s.get(start), Some(b'>' | b'@')) {
                let header_bytes = seq_buf.s.split_off(start);
                kstream.pending_header = String::from_utf8(header_bytes).ok();
                break;
            }
        }
        body_len = seq_buf.s.len();
    } else {
        loop {
            let start = seq_buf.s.len();
            let n = reader.read_until(b'\n', &mut seq_buf.s).unwrap_or(0);
            if n == 0 {
                kstream.is_eof = 1;
                break;
            }
            while seq_buf.s.len() > start && matches!(seq_buf.s.last(), Some(b'\n' | b'\r')) {
                seq_buf.s.pop();
            }
            if matches!(seq_buf.s.get(start), Some(b'+')) {
                seq_buf.s.truncate(start);
                break;
            }
        }
        body_len = seq_buf.s.len();
        // Read qual lines until we have enough bytes to match body_len.
        while qual_buf.s.len() < body_len {
            let n = reader.read_until(b'\n', &mut qual_buf.s).unwrap_or(0);
            if n == 0 {
                kstream.is_eof = 1;
                break;
            }
            while matches!(qual_buf.s.last(), Some(b'\n' | b'\r')) {
                qual_buf.s.pop();
            }
        }
        // qual length mismatch → return -2 (matches old contract).
        if qual_buf.s.len() != body_len {
            return -2;
        }
    }

    seq_buf.l = body_len;
    if !seq_buf.s.is_empty() {
        seq_buf.s.push(0);
        seq_buf.m = seq_buf.s.len();
    }
    if is_fastq {
        qual_buf.l = qual_buf.s.len();
        if !qual_buf.s.is_empty() {
            qual_buf.s.push(0);
            qual_buf.m = qual_buf.s.len();
        }
    }
    body_len as i64
}

pub fn kseq_read(seq: &mut kseq_t) -> i64 {
    if seq.f.reader.is_some() {
        return streaming_kseq_read(seq);
    }
    if seq.f.record_index >= seq.f.records.len() {
        seq.f.is_eof = 1;
        seq.f.begin = seq.f.end;
        return -1;
    }
    let record = seq.f.records[seq.f.record_index].clone();
    seq.f.record_index += 1;
    set_kstring(&mut seq.name, &record.name);
    set_kstring(&mut seq.comment, record.comment.as_deref().unwrap_or(""));
    set_kstring(&mut seq.seq, &record.seq);
    set_kstring(&mut seq.qual, record.qual.as_deref().unwrap_or(""));
    if let Some(qual) = record.qual.as_ref() {
        if qual.len() != record.seq.len() {
            return -2;
        }
    }
    i64::try_from(record.seq.len()).expect("sequence length")
}

pub fn kseq_destroy(ks: kseq_t) {
    drop(ks);
}

#[doc = "Original function: ks_getuntil:152"]
pub fn ks_getuntil(
    ks: &mut kstream_t,
    delimiter: i32,
    str_: &mut kstring_t,
    dret: Option<&mut i32>,
) -> i32 {
    if ks.input.is_empty() && ks.reader.is_some() {
        let mut bytes = Vec::new();
        if let Some(reader) = ks.reader.as_mut() {
            if reader.read_to_end(&mut bytes).is_err() {
                ks.is_eof = 1;
                return -1;
            }
        }
        ks.input = bytes;
        ks.cursor = 0;
        ks.begin = 0;
        ks.end = ks.input.len();
    }
    let bytes = &ks.input;
    if ks.cursor >= bytes.len() {
        ks.is_eof = 1;
        if let Some(dret) = dret {
            *dret = 0;
        }
        str_.l = 0;
        return -1;
    }
    let start = ks.cursor;
    let mut end = start;
    while end < bytes.len() {
        let b = bytes[end];
        let stop = match delimiter {
            KS_SEP_LINE => b == b'\n',
            KS_SEP_SPACE => b.is_ascii_whitespace(),
            KS_SEP_TAB => b.is_ascii_whitespace() && b != b' ',
            d if d > KS_SEP_MAX => b == u8::try_from(d).unwrap_or_default(),
            _ => false,
        };
        if stop {
            break;
        }
        end += 1;
    }
    let mut out = &bytes[start..end];
    if delimiter == KS_SEP_LINE && out.last() == Some(&b'\r') {
        out = &out[..out.len() - 1];
    }
    ks.cursor = if end < bytes.len() { end + 1 } else { end };
    ks.begin = ks.cursor;
    ks.end = bytes.len();
    if ks.cursor >= bytes.len() {
        ks.is_eof = 1;
    }
    set_kstring(str_, std::str::from_utf8(out).expect("kseq text utf8"));
    if let Some(dret) = dret {
        *dret = if end < bytes.len() {
            i32::from(bytes[end])
        } else {
            0
        };
    }
    i32::try_from(str_.l).expect("ks_getuntil length")
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::process::{Command, Stdio};

    use super::{
        ks_getuntil, kseq_destroy, kseq_read, kseq_t, kstream_t, kstring_t, KS_SEP_LINE,
        KS_SEP_SPACE,
    };

    #[test]
    fn kseq_read_parses_fasta_and_fastq_records() {
        let mut fasta = kseq_t::from_text(">r0 comment\nAC\nGT\n>r1\nTTA\n");
        assert_eq!(kseq_read(&mut fasta), 4);
        assert_eq!(fasta.name.as_str(), "r0");
        assert_eq!(fasta.comment.as_str(), "comment");
        assert_eq!(fasta.seq.as_str(), "ACGT");
        assert_eq!(fasta.qual.as_str(), "");
        assert_eq!(kseq_read(&mut fasta), 3);
        assert_eq!(fasta.name.as_str(), "r1");
        assert_eq!(fasta.seq.as_str(), "TTA");
        assert_eq!(kseq_read(&mut fasta), -1);

        let mut fastq = kseq_t::from_text("@q0 note\nAC\nGT\n+\nII\nJJ\n");
        assert_eq!(kseq_read(&mut fastq), 4);
        assert_eq!(fastq.name.as_str(), "q0");
        assert_eq!(fastq.comment.as_str(), "note");
        assert_eq!(fastq.seq.as_str(), "ACGT");
        assert_eq!(fastq.qual.as_str(), "IIJJ");
    }

    #[test]
    fn ks_getuntil_honors_space_and_line_delimiters() {
        let mut ks = kstream_t::from_text("name comment here\nnext");
        let mut out = kstring_t::default();
        let mut delim = 0;
        assert_eq!(
            ks_getuntil(&mut ks, KS_SEP_SPACE, &mut out, Some(&mut delim)),
            4
        );
        assert_eq!(out.as_str(), "name");
        assert_eq!(delim, i32::from(b' '));

        assert_eq!(
            ks_getuntil(&mut ks, KS_SEP_LINE, &mut out, Some(&mut delim)),
            12
        );
        assert_eq!(out.as_str(), "comment here");
        assert_eq!(delim, i32::from(b'\n'));
    }

    #[test]
    fn ks_getuntil_reads_from_streaming_reader() {
        let reader = Box::new(BufReader::new("name comment\n".as_bytes()));
        let mut ks = kstream_t::from_reader(reader);
        let mut out = kstring_t::default();
        let mut delim = 0;
        assert_eq!(
            ks_getuntil(&mut ks, KS_SEP_SPACE, &mut out, Some(&mut delim)),
            4
        );
        assert_eq!(out.as_str(), "name");
        assert_eq!(delim, i32::from(b' '));
    }

    #[test]
    fn kseq_destroy_waits_for_child_backed_stream() {
        let mut child = Command::new("rustc")
            .arg("--version")
            .stdout(Stdio::piped())
            .spawn()
            .expect("spawn rustc");
        let stdout = child.stdout.take().expect("stdout");
        let ks = kseq_t::from_reader_with_child(Box::new(BufReader::new(stdout)), child);
        kseq_destroy(ks);
    }

    #[test]
    fn kseq_read_streams_fastq_records_without_preparse() {
        let reader = Box::new(BufReader::new(
            "@r0 note\nAC\nGT\n+\nII\nJJ\n@r1\nTT\n+\nHH\n".as_bytes(),
        ));
        let mut ks = kseq_t::from_reader(reader);
        assert_eq!(kseq_read(&mut ks), 4);
        assert_eq!(ks.name.as_str(), "r0");
        assert_eq!(ks.comment.as_str(), "note");
        assert_eq!(ks.seq.as_str(), "ACGT");
        assert_eq!(ks.qual.as_str(), "IIJJ");
        assert_eq!(kseq_read(&mut ks), 2);
        assert_eq!(ks.name.as_str(), "r1");
        assert_eq!(ks.seq.as_str(), "TT");
        assert_eq!(kseq_read(&mut ks), -1);
    }
}
