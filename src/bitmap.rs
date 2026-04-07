use std::os::raw::c_int;

#[inline(always)]
pub unsafe fn test(bm: *const u8, ndx: c_int) -> u8 {
    let byte = *bm.add((ndx >> 3) as usize);
    if byte & (1 << (ndx & 0x07)) != 0 {
        1
    } else {
        0
    }
}

#[inline(always)]
pub unsafe fn clear(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p &= !(1u8 << (ndx & 0x07));
}

#[inline(always)]
pub unsafe fn set(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p |= 1u8 << (ndx & 0x07);
}

#[inline(always)]
pub unsafe fn toggle(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p ^= 1u8 << (ndx & 0x07);
}
