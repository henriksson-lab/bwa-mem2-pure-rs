#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/kstring.h` + `bwa-mem2/src/kstring.cpp`.

use std::fmt::Arguments;

// --- kstring.h ---

#[inline]
fn kroundup32(mut x: usize) -> usize {
    if x == 0 {
        return 0;
    }
    x -= 1;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    if usize::BITS > 32 {
        x |= x >> 32;
    }
    x + 1
}

#[doc = "Original struct: kstring_t (bwa-mem2/src/kstring.h)"]
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct kstring_t {
    pub l: usize,
    pub m: usize,
    pub s: Vec<u8>,
}

#[doc = "Original struct: __kstring_t (bwa-mem2/src/kstring.h)"]
pub type __kstring_t = kstring_t;

impl kstring_t {
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.s[..self.l]
    }

    #[inline]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).expect("kstring_t contains invalid UTF-8")
    }

    #[inline]
    fn ensure_allocated(&mut self, size: usize) {
        if self.m < size {
            self.m = kroundup32(size);
            self.s.resize(self.m, 0);
        }
    }
}

#[doc = "Original function: ks_resize:53"]
#[inline]
pub fn ks_resize(s: &mut kstring_t, size: usize) {
    if s.m < size {
        s.ensure_allocated(size);
    }
}

#[doc = "Original function: kputsn:62"]
#[inline]
pub fn kputsn(p: &[u8], l: i32, s: &mut kstring_t) -> i32 {
    let l = l as usize;
    if s.l + l + 1 >= s.m {
        s.ensure_allocated(s.l + l + 2);
    }
    s.s[s.l..s.l + l].copy_from_slice(&p[..l]);
    s.l += l;
    s.s[s.l] = 0;
    l as i32
}

#[doc = "Original function: kputs:75"]
#[inline]
pub fn kputs(p: &str, s: &mut kstring_t) -> i32 {
    kputsn(p.as_bytes(), p.len() as i32, s)
}

#[doc = "Original function: kputc:80"]
#[inline]
pub fn kputc(c: i32, s: &mut kstring_t) -> i32 {
    if s.l + 1 >= s.m {
        s.ensure_allocated(s.l + 2);
    }
    s.s[s.l] = c as u8;
    s.l += 1;
    s.s[s.l] = 0;
    c
}

#[inline]
fn append_u64(mut n: u64, neg: bool, s: &mut kstring_t) {
    // Up to 20 digits for u64 + 1 sign byte.
    let mut buf = [0u8; 21];
    let mut i = buf.len();
    while n >= 10 {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
    }
    i -= 1;
    buf[i] = b'0' + n as u8;
    if neg {
        i -= 1;
        buf[i] = b'-';
    }
    let len = buf.len() - i;
    if s.l + len + 1 >= s.m {
        s.ensure_allocated(s.l + len + 2);
    }
    s.s[s.l..s.l + len].copy_from_slice(&buf[i..]);
    s.l += len;
    s.s[s.l] = 0;
}

#[doc = "Original function: kputw:92"]
#[inline]
pub fn kputw(c: i32, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    let (n, neg) = if c < 0 {
        ((c as i64).unsigned_abs(), true)
    } else {
        (c as u64, false)
    };
    append_u64(n, neg, s);
    0
}

#[doc = "Original function: kputuw:109"]
#[inline]
pub fn kputuw(c: u32, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    append_u64(c as u64, false, s);
    0
}

#[doc = "Original function: kputl:126"]
#[inline]
pub fn kputl(c: i64, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    let (n, neg) = if c < 0 {
        (c.unsigned_abs(), true)
    } else {
        (c as u64, false)
    };
    append_u64(n, neg, s);
    0
}

#[cfg(test)]
mod tests_h {
    use super::{kputc, kputl, kputs, kputsn, kputuw, kputw, ks_resize, kstring_t};

    #[test]
    fn ks_resize_rounds_up_and_zero_fills() {
        let mut s = kstring_t::default();
        ks_resize(&mut s, 5);
        assert_eq!(s.m, 8);
        assert_eq!(s.s.len(), 8);
        assert!(s.s.iter().all(|&b| b == 0));
    }

    #[test]
    fn kputsn_appends_and_nul_terminates() {
        let mut s = kstring_t::default();
        assert_eq!(kputsn(b"abc", 3, &mut s), 3);
        assert_eq!(s.as_bytes(), b"abc");
        assert_eq!(s.s[s.l], 0);
    }

    #[test]
    fn put_helpers_match_expected_text() {
        let mut s = kstring_t::default();
        assert_eq!(kputs("ab", &mut s), 2);
        assert_eq!(kputc(i32::from(b'!'), &mut s), i32::from(b'!'));
        assert_eq!(kputw(-12, &mut s), 0);
        assert_eq!(kputuw(34, &mut s), 0);
        assert_eq!(kputl(-56, &mut s), 0);
        assert_eq!(s.as_str(), "ab!-1234-56");
        assert_eq!(s.s[s.l], 0);
    }

    #[test]
    fn zero_numeric_helpers_return_ascii_zero() {
        let mut s = kstring_t::default();
        assert_eq!(kputw(0, &mut s), i32::from(b'0'));
        assert_eq!(kputuw(0, &mut s), i32::from(b'0'));
        assert_eq!(kputl(0, &mut s), i32::from(b'0'));
        assert_eq!(s.as_str(), "000");
    }
}

// --- kstring.cpp ---

#[doc = "Original function: ksprintf:34"]
pub fn ksprintf(s: &mut kstring_t, args: Arguments<'_>) -> i32 {
    let rendered = args.to_string();
    let l = rendered.len();
    let remaining = s.m.saturating_sub(s.l);
    if l + 1 > remaining {
        let needed = s.l + l + 2;
        let mut rounded = needed;
        rounded -= 1;
        rounded |= rounded >> 1;
        rounded |= rounded >> 2;
        rounded |= rounded >> 4;
        rounded |= rounded >> 8;
        rounded |= rounded >> 16;
        if usize::BITS > 32 {
            rounded |= rounded >> 32;
        }
        s.m = rounded + 1;
        s.s.resize(s.m, 0);
    }
    s.s[s.l..s.l + l].copy_from_slice(rendered.as_bytes());
    s.l += l;
    if s.l == s.s.len() {
        s.s.push(0);
        s.m = s.s.len();
    } else {
        s.s[s.l] = 0;
    }
    i32::try_from(l).expect("ksprintf length overflow")
}

#[doc = "Original function: main:55"]
pub fn main() -> i32 {
    let mut s = kstring_t::default();
    let _ = ksprintf(&mut s, format_args!("abcdefg: {}", 100));
    println!("{}", s.as_str());
    0
}

#[cfg(test)]
mod tests {
    use super::ksprintf;
    use crate::bwa_mem2::kstring::kstring_t;

    #[test]
    fn ksprintf_appends_and_grows_like_c_version() {
        let mut s = kstring_t::default();
        assert_eq!(ksprintf(&mut s, format_args!("abcdefg: {}", 100)), 12);
        assert_eq!(s.as_str(), "abcdefg: 100");
        assert_eq!(s.m, 16);
        assert_eq!(s.s[s.l], 0);
    }

    #[test]
    fn ksprintf_reuses_remaining_space_when_exactly_enough() {
        let mut s = kstring_t {
            l: 2,
            m: 8,
            s: vec![0; 8],
        };
        s.s[..2].copy_from_slice(b"ab");
        assert_eq!(ksprintf(&mut s, format_args!("cdefg")), 5);
        assert_eq!(s.as_str(), "abcdefg");
        assert_eq!(s.m, 8);
        assert_eq!(s.s[s.l], 0);
    }
}
