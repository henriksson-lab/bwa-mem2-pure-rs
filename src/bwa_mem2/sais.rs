#![allow(
    dead_code,
    non_snake_case,
    non_camel_case_types,
    non_upper_case_globals
)]

//! Port of `bwa-mem2/src/sais.h`.

// --- sais.h ---

#[doc = "Original function: getCounts:52"]
pub(crate) fn get_counts(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
) {
    crate::support::stub::<()>("get_counts")
}

#[doc = "Original function: getBuckets:59"]
pub(crate) fn get_buckets(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
) {
    crate::support::stub::<()>("get_buckets")
}

#[doc = "Original function: LMSsort1:68"]
pub(crate) fn lms_sort1(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) {
    crate::support::stub::<()>("lms_sort1")
}

#[doc = "Original function: LMSpostproc1:111"]
pub(crate) fn lms_postproc1(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("lms_postproc1")
}

#[doc = "Original function: LMSsort2:159"]
pub(crate) fn lms_sort2(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) {
    crate::support::stub::<()>("lms_sort2")
}

#[doc = "Original function: LMSpostproc2:219"]
pub(crate) fn lms_postproc2(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("lms_postproc2")
}

#[doc = "Original function: induceSA:261"]
pub(crate) fn induce_sa(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) {
    crate::support::stub::<()>("induce_sa")
}

#[doc = "Original function: computeBWT:294"]
pub(crate) fn compute_bwt(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("compute_bwt")
}

#[doc = "Original function: stage1sort:335"]
pub(crate) fn stage1sort(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("stage1sort")
}

#[doc = "Original function: stage3sort:395"]
pub(crate) fn stage3sort(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
    _arg6: crate::support::Opaque,
    _arg7: crate::support::Opaque,
    _arg8: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("stage3sort")
}

#[doc = "Original function: suffixsort:426"]
pub(crate) fn suffixsort(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
    _arg4: crate::support::Opaque,
    _arg5: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("suffixsort")
}

#[doc = "Original function: saisxx:556"]
pub(crate) fn saisxx(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("saisxx")
}

trait SaisSymbol: Copy + Ord {
    fn bucket(self) -> usize;
}

impl SaisSymbol for u8 {
    #[inline]
    fn bucket(self) -> usize {
        usize::from(self)
    }
}

impl SaisSymbol for usize {
    #[inline]
    fn bucket(self) -> usize {
        self
    }
}

impl SaisSymbol for u32 {
    #[inline]
    fn bucket(self) -> usize {
        usize::try_from(self).expect("u32 bucket")
    }
}

impl SaisSymbol for i64 {
    #[inline]
    fn bucket(self) -> usize {
        usize::try_from(self).expect("i64 bucket")
    }
}

#[inline]
fn bucket_sizes<T: SaisSymbol>(text: &[T], alphabet: usize) -> Vec<usize> {
    let mut sizes = vec![0_usize; alphabet];
    for &c in text {
        sizes[c.bucket()] += 1;
    }
    sizes
}

#[inline]
fn bucket_heads(sizes: &[usize]) -> Vec<usize> {
    let mut heads = vec![0_usize; sizes.len()];
    let mut sum = 0_usize;
    for (idx, &size) in sizes.iter().enumerate() {
        heads[idx] = sum;
        sum += size;
    }
    heads
}

#[inline]
fn bucket_tails(sizes: &[usize]) -> Vec<usize> {
    let mut tails = vec![0_usize; sizes.len()];
    let mut sum = 0_usize;
    for (idx, &size) in sizes.iter().enumerate() {
        sum += size;
        tails[idx] = sum;
    }
    tails
}

#[inline]
fn classify_s_types<T: SaisSymbol>(text: &[T]) -> Vec<u8> {
    let n = text.len();
    let mut is_s = vec![0_u8; n];
    is_s[n - 1] = 1;
    for i in (0..n - 1).rev() {
        let c0 = unsafe { *text.get_unchecked(i) };
        let c1 = unsafe { *text.get_unchecked(i + 1) };
        let next_is_s = unsafe { *is_s.get_unchecked(i + 1) } != 0;
        unsafe {
            *is_s.get_unchecked_mut(i) = u8::from(c0 < c1 || (c0 == c1 && next_is_s));
        }
    }
    is_s
}

#[inline]
fn is_lms_pos(is_s: &[u8], pos: usize) -> bool {
    pos > 0 && is_s[pos] != 0 && is_s[pos - 1] == 0
}

fn induce_sa_impl<T: SaisSymbol>(
    text: &[T],
    alphabet: usize,
    is_s: &[u8],
    lms_order: &[u32],
) -> Vec<i32> {
    let n = text.len();
    let sizes = bucket_sizes(text, alphabet);
    let mut sa = vec![-1_i32; n];

    let mut tails = bucket_tails(&sizes);
    for &pos in lms_order.iter().rev() {
        let pos = pos as usize;
        let c = unsafe { text.get_unchecked(pos).bucket() };
        tails[c] -= 1;
        unsafe {
            *sa.get_unchecked_mut(tails[c]) = pos as i32;
        }
    }

    let mut heads = bucket_heads(&sizes);
    for i in 0..n {
        let j = unsafe { *sa.get_unchecked(i) };
        if j > 0 {
            let prev = (j as usize) - 1;
            if unsafe { *is_s.get_unchecked(prev) } == 0 {
                let c = unsafe { text.get_unchecked(prev).bucket() };
                unsafe {
                    *sa.get_unchecked_mut(heads[c]) = prev as i32;
                }
                heads[c] += 1;
            }
        }
    }

    let mut tails = bucket_tails(&sizes);
    for i in (0..n).rev() {
        let j = unsafe { *sa.get_unchecked(i) };
        if j > 0 {
            let prev = (j as usize) - 1;
            if unsafe { *is_s.get_unchecked(prev) } != 0 {
                let c = unsafe { text.get_unchecked(prev).bucket() };
                tails[c] -= 1;
                unsafe {
                    *sa.get_unchecked_mut(tails[c]) = prev as i32;
                }
            }
        }
    }

    sa
}

fn same_lms_substring<T: SaisSymbol>(text: &[T], is_s: &[u8], a: usize, b: usize) -> bool {
    if a == b {
        return true;
    }
    let n = text.len();
    let mut offset = 0_usize;
    loop {
        let a_pos = a + offset;
        let b_pos = b + offset;
        if a_pos >= n || b_pos >= n {
            return a_pos == b_pos;
        }
        if text[a_pos] != text[b_pos] {
            return false;
        }
        let a_lms = offset > 0 && is_lms_pos(is_s, a_pos);
        let b_lms = offset > 0 && is_lms_pos(is_s, b_pos);
        if a_lms || b_lms {
            return a_lms && b_lms;
        }
        offset += 1;
    }
}

fn sais_with_sentinel<T: SaisSymbol>(text: &[T], alphabet: usize) -> Vec<i32> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![0];
    }

    let is_s = classify_s_types(text);
    let mut lms_positions = Vec::new();
    for pos in 1..n {
        if unsafe { *is_s.get_unchecked(pos) != 0 && *is_s.get_unchecked(pos - 1) == 0 } {
            lms_positions.push(u32::try_from(pos).expect("lms position"));
        }
    }
    let sa = induce_sa_impl(text, alphabet, &is_s, &lms_positions);
    let mut sorted_lms = Vec::with_capacity(lms_positions.len());
    for &pos in &sa {
        if pos >= 0 {
            let pos_us = pos as usize;
            if is_lms_pos(&is_s, pos_us) {
                sorted_lms.push(u32::try_from(pos_us).expect("lms pos"));
            }
        }
    }
    drop(sa);

    let mut name = -1_i32;
    let mut prev = None;
    let mut named_lms = Vec::with_capacity(sorted_lms.len());
    for &pos in &sorted_lms {
        let pos_us = usize::try_from(pos).expect("lms pos");
        if prev.map_or(true, |p| !same_lms_substring(text, &is_s, p, pos_us)) {
            name += 1;
        }
        named_lms.push((pos, u32::try_from(name).expect("lms name")));
        prev = Some(pos_us);
    }
    drop(sorted_lms);
    let name_count = usize::try_from(name + 1).expect("name count");

    named_lms.sort_unstable_by_key(|&(pos, _)| pos);
    let mut summary = Vec::with_capacity(lms_positions.len());
    let mut unique_order = (name_count == lms_positions.len()).then(|| vec![0_u32; name_count]);
    for (idx, (&pos, &(named_pos, lms_name))) in
        lms_positions.iter().zip(named_lms.iter()).enumerate()
    {
        debug_assert_eq!(pos, named_pos);
        summary.push(lms_name);
        if let Some(order) = unique_order.as_mut() {
            order[usize::try_from(lms_name).expect("lms name")] =
                u32::try_from(idx).expect("lms idx");
        }
    }
    drop(named_lms);

    let (ordered_lms, is_s) = if name_count == lms_positions.len() {
        drop(summary);
        (
            unique_order
                .take()
                .expect("unique order")
                .into_iter()
                .map(|idx| lms_positions[usize::try_from(idx).expect("lms idx")])
                .collect::<Vec<_>>(),
            is_s,
        )
    } else {
        drop(unique_order);
        drop(is_s);
        let summary_sa = sais_with_sentinel(&summary, name_count);
        drop(summary);
        let ordered_lms = summary_sa
            .into_iter()
            .map(|idx| lms_positions[usize::try_from(idx).expect("summary sa idx")])
            .collect::<Vec<_>>();
        (ordered_lms, classify_s_types(text))
    };

    induce_sa_impl(text, alphabet, &is_s, &ordered_lms)
}

pub fn sais_suffixes_u32(binary_ref_seq: &[u8]) -> Vec<u32> {
    let n = binary_ref_seq.len();
    assert!(
        n < i32::MAX as usize,
        "SA-IS u32 path requires n < i32::MAX"
    );
    let mut text = Vec::with_capacity(n + 1);
    text.extend(binary_ref_seq.iter().map(|&base| base + 1));
    text.push(0);
    let sa = sais_with_sentinel(&text, 5);
    debug_assert_eq!(usize::try_from(sa[0]).expect("sentinel"), n);
    sa.into_iter()
        .skip(1)
        .map(|pos| u32::try_from(pos).expect("suffix coordinate"))
        .collect()
}

fn upstream_counts<T: SaisSymbol>(text: &[T], k: usize) -> Vec<i64> {
    let mut counts = vec![0_i64; k];
    for &c in text {
        counts[c.bucket()] += 1;
    }
    counts
}

fn upstream_buckets(counts: &[i64], end: bool) -> Vec<i64> {
    let mut buckets = vec![0_i64; counts.len()];
    upstream_buckets_into(counts, &mut buckets, end);
    buckets
}

fn upstream_buckets_into(counts: &[i64], buckets: &mut [i64], end: bool) {
    let mut sum = 0_i64;
    if end {
        for (idx, &count) in counts.iter().enumerate() {
            sum += count;
            buckets[idx] = sum;
        }
    } else {
        for (idx, &count) in counts.iter().enumerate() {
            sum += count;
            buckets[idx] = sum - count;
        }
    }
}

fn upstream_lms_sort1<T: SaisSymbol>(
    text: &[T],
    sa: &mut [i64],
    counts: &[i64],
    buckets: &mut [i64],
) {
    let n = i64::try_from(text.len()).expect("len");
    upstream_buckets_into(counts, buckets, false);
    let mut j = n - 1;
    let mut c1 = text[j as usize].bucket();
    let mut b = buckets[c1] as usize;
    j -= 1;
    sa[b] = if text[j as usize].bucket() < c1 {
        !j
    } else {
        j
    };
    b += 1;

    for i in 0..n as usize {
        j = sa[i];
        if j > 0 {
            let c0 = text[j as usize].bucket();
            if c0 != c1 {
                buckets[c1] = b as i64;
                c1 = c0;
                b = buckets[c1] as usize;
            }
            j -= 1;
            sa[b] = if text[j as usize].bucket() < c1 {
                !j
            } else {
                j
            };
            b += 1;
            sa[i] = 0;
        } else if j < 0 {
            sa[i] = !j;
        }
    }

    upstream_buckets_into(counts, buckets, true);
    c1 = 0;
    b = buckets[c1] as usize;
    let mut i = n - 1;
    loop {
        j = sa[i as usize];
        if j > 0 {
            let c0 = text[j as usize].bucket();
            if c0 != c1 {
                buckets[c1] = b as i64;
                c1 = c0;
                b = buckets[c1] as usize;
            }
            j -= 1;
            b -= 1;
            sa[b] = if text[j as usize].bucket() > c1 {
                !(j + 1)
            } else {
                j
            };
            sa[i as usize] = 0;
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }
}

fn upstream_lms_postproc1<T: SaisSymbol>(text: &[T], sa: &mut [i64], n: i64, m: i64) -> i64 {
    let mut i = 0_i64;
    while sa[i as usize] < 0 {
        sa[i as usize] = !sa[i as usize];
        i += 1;
    }
    if i < m {
        let mut d = i;
        i += 1;
        loop {
            let p = sa[i as usize];
            if p < 0 {
                sa[d as usize] = !p;
                sa[i as usize] = 0;
                d += 1;
                if d == m {
                    break;
                }
            }
            i += 1;
        }
    }

    let mut i = n - 1;
    let mut j = n - 1;
    let mut c0 = text[(n - 1) as usize].bucket();
    loop {
        i -= 1;
        if i < 0 {
            break;
        }
        let c1 = c0;
        c0 = text[i as usize].bucket();
        if c0 < c1 {
            break;
        }
    }
    while i >= 0 {
        loop {
            i -= 1;
            if i < 0 {
                break;
            }
            let c1 = c0;
            c0 = text[i as usize].bucket();
            if c0 > c1 {
                break;
            }
        }
        if i >= 0 {
            sa[(m + ((i + 1) >> 1)) as usize] = j - i;
            j = i + 1;
            loop {
                i -= 1;
                if i < 0 {
                    break;
                }
                let c1 = c0;
                c0 = text[i as usize].bucket();
                if c0 < c1 {
                    break;
                }
            }
        }
    }

    let mut name = 0_i64;
    let mut q = n;
    let mut qlen = 0_i64;
    for i in 0..m as usize {
        let p = sa[i];
        let plen = sa[(m + (p >> 1)) as usize];
        let mut diff = true;
        if plen == qlen && q + plen < n {
            let mut offset = 0_i64;
            while offset < plen
                && text[(p + offset) as usize].bucket() == text[(q + offset) as usize].bucket()
            {
                offset += 1;
            }
            if offset == plen {
                diff = false;
            }
        }
        if diff {
            name += 1;
            q = p;
            qlen = plen;
        }
        sa[(m + (p >> 1)) as usize] = name;
    }

    name
}

fn upstream_lms_sort2<T: SaisSymbol>(
    text: &[T],
    sa: &mut [i64],
    counts: &[i64],
    buckets: &mut [i64],
    d_seen: &mut [i64],
) {
    let n = i64::try_from(text.len()).expect("len");
    upstream_buckets_into(counts, buckets, false);
    let mut j = n - 1;
    let mut c1 = text[j as usize].bucket();
    let mut b = buckets[c1] as usize;
    j -= 1;
    let mut t = i64::from(text[j as usize].bucket() < c1);
    j += n;
    sa[b] = if (t & 1) != 0 { !j } else { j };
    b += 1;

    let mut d = 0_i64;
    for i in 0..n as usize {
        j = sa[i];
        if j > 0 {
            if n <= j {
                d += 1;
                j -= n;
            }
            let c0 = text[j as usize].bucket();
            if c0 != c1 {
                buckets[c1] = b as i64;
                c1 = c0;
                b = buckets[c1] as usize;
            }
            j -= 1;
            t = ((c0 as i64) << 1) | i64::from(text[j as usize].bucket() < c1);
            if d_seen[t as usize] != d {
                j += n;
                d_seen[t as usize] = d;
            }
            sa[b] = if (t & 1) != 0 { !j } else { j };
            b += 1;
            sa[i] = 0;
        } else if j < 0 {
            sa[i] = !j;
        }
    }

    let mut i = n - 1;
    while i >= 0 {
        let iu = i as usize;
        if sa[iu] > 0 && sa[iu] < n {
            sa[iu] += n;
            let mut jj = i - 1;
            while sa[jj as usize] < n {
                jj -= 1;
            }
            sa[jj as usize] -= n;
            i = jj;
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }

    upstream_buckets_into(counts, buckets, true);
    d += 1;
    c1 = 0;
    b = buckets[c1] as usize;
    let mut i = n - 1;
    loop {
        j = sa[i as usize];
        if j > 0 {
            if n <= j {
                d += 1;
                j -= n;
            }
            let c0 = text[j as usize].bucket();
            if c0 != c1 {
                buckets[c1] = b as i64;
                c1 = c0;
                b = buckets[c1] as usize;
            }
            j -= 1;
            t = ((c0 as i64) << 1) | i64::from(text[j as usize].bucket() > c1);
            if d_seen[t as usize] != d {
                j += n;
                d_seen[t as usize] = d;
            }
            b -= 1;
            sa[b] = if (t & 1) != 0 { !(j + 1) } else { j };
            sa[i as usize] = 0;
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }
}

fn upstream_lms_postproc2(sa: &mut [i64], n: i64, m: i64) -> i64 {
    let mut i = 0_i64;
    let mut name = 0_i64;
    while sa[i as usize] < 0 {
        let j = !sa[i as usize];
        if n <= j {
            name += 1;
        }
        sa[i as usize] = j;
        i += 1;
    }
    if i < m {
        let mut d = i;
        i += 1;
        loop {
            let j0 = sa[i as usize];
            if j0 < 0 {
                let j = !j0;
                if n <= j {
                    name += 1;
                }
                sa[d as usize] = j;
                sa[i as usize] = 0;
                d += 1;
                if d == m {
                    break;
                }
            }
            i += 1;
        }
    }
    if name < m {
        let mut i = m - 1;
        let mut d = name + 1;
        loop {
            let mut j = sa[i as usize];
            if n <= j {
                j -= n;
                d -= 1;
            }
            sa[(m + (j >> 1)) as usize] = d;
            if i == 0 {
                break;
            }
            i -= 1;
        }
    } else {
        for i in 0..m as usize {
            let j = sa[i];
            if n <= j {
                sa[i] = j - n;
            }
        }
    }
    name
}

fn upstream_induce_sa<T: SaisSymbol>(
    text: &[T],
    sa: &mut [i64],
    counts: &[i64],
    buckets: &mut [i64],
) {
    let n = i64::try_from(text.len()).expect("len");

    upstream_buckets_into(counts, buckets, false);
    let mut j = n - 1;
    let mut c1 = unsafe { text.get_unchecked(j as usize).bucket() };
    let mut b = unsafe { *buckets.get_unchecked(c1) as usize };
    unsafe {
        *sa.get_unchecked_mut(b) = if j > 0 && text.get_unchecked((j - 1) as usize).bucket() < c1 {
            !j
        } else {
            j
        };
    }
    b += 1;

    for i in 0..n as usize {
        j = unsafe { *sa.get_unchecked(i) };
        unsafe {
            *sa.get_unchecked_mut(i) = !j;
        }
        if j > 0 {
            j -= 1;
            let c0 = unsafe { text.get_unchecked(j as usize).bucket() };
            if c0 != c1 {
                unsafe {
                    *buckets.get_unchecked_mut(c1) = b as i64;
                }
                c1 = c0;
                b = unsafe { *buckets.get_unchecked(c1) as usize };
            }
            unsafe {
                *sa.get_unchecked_mut(b) =
                    if j > 0 && text.get_unchecked((j - 1) as usize).bucket() < c1 {
                        !j
                    } else {
                        j
                    };
            }
            b += 1;
        }
    }

    upstream_buckets_into(counts, buckets, true);
    c1 = 0;
    b = unsafe { *buckets.get_unchecked(c1) as usize };
    let mut i = n - 1;
    loop {
        j = unsafe { *sa.get_unchecked(i as usize) };
        if j > 0 {
            j -= 1;
            let c0 = unsafe { text.get_unchecked(j as usize).bucket() };
            if c0 != c1 {
                unsafe {
                    *buckets.get_unchecked_mut(c1) = b as i64;
                }
                c1 = c0;
                b = unsafe { *buckets.get_unchecked(c1) as usize };
            }
            b -= 1;
            unsafe {
                *sa.get_unchecked_mut(b) =
                    if j == 0 || text.get_unchecked((j - 1) as usize).bucket() > c1 {
                        !j
                    } else {
                        j
                    };
            }
        } else {
            unsafe {
                *sa.get_unchecked_mut(i as usize) = !j;
            }
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }
}

fn upstream_stage3sort<T: SaisSymbol>(
    text: &[T],
    sa: &mut [i64],
    m: i64,
    k: usize,
    flags: u32,
    counts: &mut Vec<i64>,
    buckets: &mut Vec<i64>,
) {
    let n = i64::try_from(text.len()).expect("len");
    if (flags & 8) != 0 {
        counts.clear();
        counts.resize(k, 0);
        for &c in text {
            counts[c.bucket()] += 1;
        }
    }
    if buckets.len() != k {
        buckets.resize(k, 0);
    }
    if m > 1 {
        upstream_buckets_into(counts, buckets, true);
        let mut i = m - 1;
        let mut j = n;
        let mut p = unsafe { *sa.get_unchecked((m - 1) as usize) };
        let mut c1 = unsafe { text.get_unchecked(p as usize).bucket() };
        while i >= 0 {
            let c0 = c1;
            let q = unsafe { *buckets.get_unchecked(c0) };
            while q < j {
                j -= 1;
                unsafe {
                    *sa.get_unchecked_mut(j as usize) = 0;
                }
            }
            loop {
                j -= 1;
                unsafe {
                    *sa.get_unchecked_mut(j as usize) = p;
                }
                i -= 1;
                if i < 0 {
                    break;
                }
                p = unsafe { *sa.get_unchecked(i as usize) };
                c1 = unsafe { text.get_unchecked(p as usize).bucket() };
                if c1 != c0 {
                    break;
                }
            }
        }
        while j > 0 {
            j -= 1;
            unsafe {
                *sa.get_unchecked_mut(j as usize) = 0;
            }
        }
    }
    upstream_induce_sa(text, sa, counts, buckets);
}

fn upstream_stage1sort<T: SaisSymbol>(
    text: &[T],
    sa: &mut [i64],
    k: usize,
    flags: u32,
) -> (i64, i64, Vec<i64>, Vec<i64>) {
    let n = i64::try_from(text.len()).expect("len");
    let counts = upstream_counts(text, k);
    let mut buckets = upstream_buckets(&counts, true);
    sa.fill(0);

    let mut b = n - 1;
    let mut i = n - 1;
    let mut j = n;
    let mut m = 0_i64;
    let mut c0 = text[(n - 1) as usize].bucket();
    loop {
        i -= 1;
        if i < 0 {
            break;
        }
        let c1 = c0;
        c0 = text[i as usize].bucket();
        if c0 < c1 {
            break;
        }
    }
    while i >= 0 {
        let mut lms_bucket = 0_usize;
        loop {
            i -= 1;
            if i < 0 {
                break;
            }
            let c1 = c0;
            c0 = text[i as usize].bucket();
            if c0 > c1 {
                lms_bucket = c1;
                break;
            }
        }
        if i >= 0 {
            sa[b as usize] = j;
            buckets[lms_bucket] -= 1;
            b = buckets[lms_bucket];
            j = i;
            m += 1;
            loop {
                i -= 1;
                if i < 0 {
                    break;
                }
                let c1 = c0;
                c0 = text[i as usize].bucket();
                if c0 < c1 {
                    break;
                }
            }
        }
    }
    sa[(n - 1) as usize] = 0;

    let name = if m > 1 {
        if (flags & (16 | 32)) != 0 {
            if j + 1 < n {
                buckets[text[(j + 1) as usize].bucket()] += 1;
            }
            let mut d_seen = vec![0_i64; k * 2];
            let mut end = 0_i64;
            for idx in 0..k {
                end += counts[idx];
                let bucket_pos = buckets[idx] as usize;
                if buckets[idx] != end && bucket_pos < sa.len() && sa[bucket_pos] != 0 {
                    sa[bucket_pos] += n;
                }
                d_seen[idx] = 0;
                d_seen[idx + k] = 0;
            }
            upstream_lms_sort2(text, sa, &counts, &mut buckets, &mut d_seen);
            upstream_lms_postproc2(sa, n, m)
        } else {
            upstream_lms_sort1(text, sa, &counts, &mut buckets);
            upstream_lms_postproc1(text, sa, n, m)
        }
    } else if m == 1 {
        sa[b as usize] = j + 1;
        1
    } else {
        0
    };

    (m, name, counts, buckets)
}

fn upstream_suffixsort_into<T: SaisSymbol>(text: &[T], sa_work: &mut [i64], k: usize, fs: i64) {
    let n = i64::try_from(text.len()).expect("len");
    if n == 0 {
        return;
    }
    if n == 1 {
        sa_work[0] = 0;
        return;
    }

    debug_assert!(sa_work.len() >= usize::try_from(n + fs).expect("workspace len"));
    let sa = &mut sa_work[..n as usize];
    let mut flags = if k <= 256 {
        if i64::try_from(k).expect("k") <= fs {
            1_u32
        } else {
            3_u32
        }
    } else if i64::try_from(k).expect("k") <= fs {
        if i64::try_from(k * 2).expect("2k") <= fs {
            0_u32
        } else if k <= 1024 {
            2_u32
        } else {
            64 | 8
        }
    } else {
        4 | 8
    };
    if n <= i64::MAX / 2 && 2 <= n / i64::try_from(k).expect("k") {
        if (flags & 1) != 0 {
            flags |= if i64::try_from(k * 2).expect("2k") <= fs - i64::try_from(k).expect("k") {
                32
            } else {
                16
            };
        } else if flags == 0
            && i64::try_from(k * 2).expect("2k") <= fs - i64::try_from(k * 2).expect("2k")
        {
            flags |= 32;
        }
    }
    let (m, name, mut counts, mut buckets) = upstream_stage1sort(text, sa, k, flags);

    if name < m {
        let mut newfs = (n + fs) - (m * 2);
        if (flags & (1 | 4 | 64)) == 0 {
            if i64::try_from(k).expect("k") + name <= newfs {
                newfs -= i64::try_from(k).expect("k");
            } else {
                flags |= 8;
            }
        }
        let ra_start = usize::try_from(m + newfs).expect("ra start");
        let m_usize = usize::try_from(m).expect("m");
        let mut j = m - 1;
        let mut i = m + (n >> 1) - 1;
        while i >= m {
            let v = sa_work[i as usize];
            if v != 0 {
                sa_work[ra_start + j as usize] = v - 1;
                if j == 0 {
                    break;
                }
                j -= 1;
            }
            if i == m {
                break;
            }
            i -= 1;
        }

        {
            let (recursive_workspace, reduced_text_region) = sa_work.split_at_mut(ra_start);
            let reduced_text = &reduced_text_region[..m_usize];
            upstream_suffixsort_into(reduced_text, recursive_workspace, name as usize, newfs);
        }

        let mut j = m - 1;
        let mut i = n - 1;
        let mut c0 = text[(n - 1) as usize].bucket();
        loop {
            i -= 1;
            if i < 0 {
                break;
            }
            let c1 = c0;
            c0 = text[i as usize].bucket();
            if c0 < c1 {
                break;
            }
        }
        while i >= 0 {
            loop {
                i -= 1;
                if i < 0 {
                    break;
                }
                let c1 = c0;
                c0 = text[i as usize].bucket();
                if c0 > c1 {
                    break;
                }
            }
            if i >= 0 {
                sa_work[ra_start + j as usize] = i + 1;
                if j == 0 {
                    break;
                }
                j -= 1;
                loop {
                    i -= 1;
                    if i < 0 {
                        break;
                    }
                    let c1 = c0;
                    c0 = text[i as usize].bucket();
                    if c0 < c1 {
                        break;
                    }
                }
            }
        }

        for i in 0..m_usize {
            sa_work[i] = sa_work[ra_start + sa_work[i] as usize];
        }

        upstream_stage3sort(
            text,
            &mut sa_work[..n as usize],
            m,
            k,
            flags,
            &mut counts,
            &mut buckets,
        );
        return;
    }

    upstream_stage3sort(
        text,
        &mut sa_work[..n as usize],
        m,
        k,
        flags,
        &mut counts,
        &mut buckets,
    );
}

fn upstream_suffixsort_with_fs<T: SaisSymbol>(text: &[T], k: usize, fs: i64) -> Vec<i64> {
    let n = i64::try_from(text.len()).expect("len");
    let mut sa = vec![0_i64; usize::try_from(n + fs).expect("workspace len")];
    upstream_suffixsort_into(text, &mut sa, k, fs);
    sa.truncate(usize::try_from(n).expect("len"));
    sa
}

fn upstream_suffixsort2<T: SaisSymbol>(text: &[T], k: usize) -> Vec<i64> {
    upstream_suffixsort_with_fs(text, k, 0)
}

pub fn sais_suffixes_u32_upstream_stage1_probe(binary_ref_seq: &[u8]) -> Vec<u32> {
    let n = binary_ref_seq.len();
    if n == 0 {
        return Vec::new();
    }
    let text: Vec<u8> = binary_ref_seq.iter().map(|&base| base + 1).collect();
    let mut sa = vec![0_i64; n];
    let (_m, _name, _counts, _buckets) = upstream_stage1sort(&text, &mut sa, 256, 0);
    sa.into_iter()
        .filter(|&x| x >= 0 && (x as usize) < n)
        .map(|x| x as u32)
        .collect()
}

pub fn sais_suffixes_i64_upstream_port(binary_ref_seq: &[u8]) -> Vec<i64> {
    let text: Vec<u8> = binary_ref_seq.iter().map(|&base| base + 1).collect();
    sais_suffixes_i64_upstream_port_mapped(&text)
}

pub fn sais_suffixes_i64_upstream_port_mapped(mapped_ref_seq: &[u8]) -> Vec<i64> {
    debug_assert!(
        mapped_ref_seq.iter().all(|&base| base > 0),
        "mapped SA-IS text must reserve zero below the alphabet"
    );
    upstream_suffixsort2(mapped_ref_seq, 256)
}

pub fn sais_suffixes_u32_upstream_port(binary_ref_seq: &[u8]) -> Vec<u32> {
    sais_suffixes_i64_upstream_port(binary_ref_seq)
        .into_iter()
        .map(|pos| u32::try_from(pos).expect("suffix coordinate"))
        .collect()
}

#[cfg(test)]
mod real_sais_tests {
    use super::{
        sais_suffixes_i64_upstream_port_mapped, sais_suffixes_u32, sais_suffixes_u32_upstream_port,
    };

    fn naive_suffixes(text: &[u8]) -> Vec<u32> {
        let mut suffixes: Vec<u32> = (0..u32::try_from(text.len()).expect("len")).collect();
        suffixes.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
        suffixes
    }

    #[test]
    fn sais_suffixes_match_naive_small_cases() {
        let cases: &[&[u8]] = &[
            b"",
            b"A",
            b"AAAA",
            b"ACGT",
            b"GATAGACA",
            b"banana",
            b"mississippi",
            b"ACGTCGTCGAAAAACGT",
        ];
        for case in cases {
            let mapped: Vec<u8> = case.iter().map(|b| b % 4).collect();
            assert_eq!(
                sais_suffixes_u32(&mapped),
                naive_suffixes(&mapped),
                "{case:?}"
            );
        }
    }

    #[test]
    fn sais_suffixes_match_naive_exhaustive_short_dna() {
        let mut text = Vec::new();
        for len in 0..=8 {
            let total = 4_usize.pow(len);
            for mut word in 0..total {
                text.clear();
                for _ in 0..len {
                    text.push((word & 3) as u8);
                    word >>= 2;
                }
                assert_eq!(sais_suffixes_u32(&text), naive_suffixes(&text), "{text:?}");
            }
        }
    }

    #[test]
    fn sais_suffixes_match_naive_repetitive_and_lcg_cases() {
        let mut cases = vec![
            vec![0; 128],
            (0..128).map(|i| (i & 3) as u8).collect::<Vec<_>>(),
            (0..128)
                .map(|i| ((i * 3 + 1) & 3) as u8)
                .collect::<Vec<_>>(),
        ];

        let mut state = 0x1234_5678_u32;
        for len in 1..=96 {
            let mut text = Vec::with_capacity(len);
            for _ in 0..len {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                text.push(((state >> 29) & 3) as u8);
            }
            cases.push(text);
        }

        for text in cases {
            assert_eq!(sais_suffixes_u32(&text), naive_suffixes(&text), "{text:?}");
        }
    }

    #[test]
    fn upstream_port_suffixes_match_naive_small_cases() {
        let cases: &[&[u8]] = &[
            b"A",
            b"AAAA",
            b"ACGT",
            b"GATAGACA",
            b"banana",
            b"mississippi",
            b"ACGTCGTCGAAAAACGT",
        ];
        for case in cases {
            let mapped: Vec<u8> = case.iter().map(|b| b % 4).collect();
            assert_eq!(
                sais_suffixes_u32_upstream_port(&mapped),
                naive_suffixes(&mapped),
                "{case:?}"
            );
        }
    }

    #[test]
    fn upstream_port_sorts_original_ascii_reference_order() {
        let cases: &[&[u8]] = &[
            b"ATN",
            b"ANTT",
            b"ACGTNNNT",
            b"TNGCATNA",
            b"ACGTMRWSYKVHDBN",
        ];
        for case in cases {
            let actual = sais_suffixes_i64_upstream_port_mapped(case)
                .into_iter()
                .map(|pos| u32::try_from(pos).expect("suffix coordinate"))
                .collect::<Vec<_>>();
            assert_eq!(actual, naive_suffixes(case), "{case:?}");
        }
    }

    #[test]
    fn upstream_port_suffixes_match_naive_exhaustive_short_dna() {
        let mut text = Vec::new();
        for len in 0..=8 {
            let total = 4_usize.pow(len);
            for mut word in 0..total {
                text.clear();
                for _ in 0..len {
                    text.push((word & 3) as u8);
                    word >>= 2;
                }
                assert_eq!(
                    sais_suffixes_u32_upstream_port(&text),
                    naive_suffixes(&text),
                    "{text:?}"
                );
            }
        }
    }

    #[test]
    fn upstream_port_suffixes_match_naive_repetitive_and_lcg_cases() {
        let mut cases = vec![
            vec![0; 128],
            (0..128).map(|i| (i & 3) as u8).collect::<Vec<_>>(),
            (0..128)
                .map(|i| ((i * 3 + 1) & 3) as u8)
                .collect::<Vec<_>>(),
        ];

        let mut state = 0x1234_5678_u32;
        for len in 1..=96 {
            let mut text = Vec::with_capacity(len);
            for _ in 0..len {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                text.push(((state >> 29) & 3) as u8);
            }
            cases.push(text);
        }

        for text in cases {
            assert_eq!(
                sais_suffixes_u32_upstream_port(&text),
                naive_suffixes(&text),
                "{text:?}"
            );
        }
    }
}

#[doc = "Original function: saisxx_bwt:578"]
pub(crate) fn saisxx_bwt(
    _arg0: crate::support::Opaque,
    _arg1: crate::support::Opaque,
    _arg2: crate::support::Opaque,
    _arg3: crate::support::Opaque,
) -> crate::support::Opaque {
    crate::support::stub::<crate::support::Opaque>("saisxx_bwt")
}
