#![allow(dead_code, non_snake_case, non_camel_case_types, non_upper_case_globals)]

//! Generated scaffold for `bwa-mem2/src/kstring.h`.

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
    pub fn as_bytes(&self) -> &[u8] {
        &self.s[..self.l]
    }

    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).expect("kstring_t contains invalid UTF-8")
    }

    fn ensure_allocated(&mut self, size: usize) {
        if self.m < size {
            self.m = kroundup32(size);
            self.s.resize(self.m, 0);
        }
    }
}

#[doc = "Original function: ks_resize:53"]
pub fn ks_resize(s: &mut kstring_t, size: usize) {
    if s.m < size {
        s.ensure_allocated(size);
    }
}

#[doc = "Original function: kputsn:62"]
pub fn kputsn(p: &[u8], l: i32, s: &mut kstring_t) -> i32 {
    let l = usize::try_from(l).expect("kputsn length must be non-negative");
    if s.l + l + 1 >= s.m {
        s.ensure_allocated(s.l + l + 2);
    }
    s.s[s.l..s.l + l].copy_from_slice(&p[..l]);
    s.l += l;
    s.s[s.l] = 0;
    i32::try_from(l).expect("kputsn length overflow")
}

#[doc = "Original function: kputs:75"]
pub fn kputs(p: &str, s: &mut kstring_t) -> i32 {
    kputsn(p.as_bytes(), i32::try_from(p.len()).expect("kputs length overflow"), s)
}

#[doc = "Original function: kputc:80"]
pub fn kputc(c: i32, s: &mut kstring_t) -> i32 {
    if s.l + 1 >= s.m {
        s.ensure_allocated(s.l + 2);
    }
    s.s[s.l] = c as u8;
    s.l += 1;
    s.s[s.l] = 0;
    c
}

fn append_decimal(text: &str, s: &mut kstring_t) {
    let len = text.len();
    if s.l + len + 1 >= s.m {
        s.ensure_allocated(s.l + len + 2);
    }
    s.s[s.l..s.l + len].copy_from_slice(text.as_bytes());
    s.l += len;
    s.s[s.l] = 0;
}

#[doc = "Original function: kputw:92"]
pub fn kputw(c: i32, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    append_decimal(&c.to_string(), s);
    0
}

#[doc = "Original function: kputuw:109"]
pub fn kputuw(c: u32, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    append_decimal(&c.to_string(), s);
    0
}

#[doc = "Original function: kputl:126"]
pub fn kputl(c: i64, s: &mut kstring_t) -> i32 {
    if c == 0 {
        return kputc(i32::from(b'0'), s);
    }
    append_decimal(&c.to_string(), s);
    0
}

#[cfg(test)]
mod tests {
    use super::{kputc, kputs, kputl, kputsn, kputuw, kputw, ks_resize, kstring_t};

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
