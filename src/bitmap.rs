use std::os::raw::c_int;

#[no_mangle]
pub unsafe extern "C" fn test(bm: *const u8, ndx: c_int) -> u8 {
    let byte = *bm.add((ndx >> 3) as usize);
    if byte & (1 << (ndx & 0x07)) != 0 {
        1
    } else {
        0
    }
}

#[no_mangle]
pub unsafe extern "C" fn clear(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p &= !(1u8 << (ndx & 0x07));
}

#[no_mangle]
pub unsafe extern "C" fn set(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p |= 1u8 << (ndx & 0x07);
}

#[no_mangle]
pub unsafe extern "C" fn toggle(bm: *mut u8, ndx: c_int) {
    let p = bm.add((ndx >> 3) as usize);
    *p ^= 1u8 << (ndx & 0x07);
}
