#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Generated scaffold for `bwa-mem2/src/kstring.cpp`.

use std::fmt::Arguments;

use crate::generated::kstring_h::kstring_t;

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
    use crate::generated::kstring_h::kstring_t;

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
