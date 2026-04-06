use crate::training_data;
use crate::types::Training;
use std::os::raw::{c_char, c_int};

#[no_mangle]
pub unsafe extern "C" fn read_training_file(
    path: *const c_char,
    tinf: *mut Training,
) -> c_int {
    let fh = libc::fopen(path, b"rb\0".as_ptr() as *const c_char);
    if fh.is_null() {
        return 1;
    }
    let rv = libc::fread(
        tinf as *mut libc::c_void,
        std::mem::size_of::<Training>(),
        1,
        fh,
    );
    libc::fclose(fh);
    if rv != 1 {
        return -1;
    }
    0
}

#[no_mangle]
pub unsafe extern "C" fn write_training_file(
    path: *const c_char,
    tinf: *mut Training,
) -> c_int {
    let fh = libc::fopen(path, b"wb\0".as_ptr() as *const c_char);
    if fh.is_null() {
        return -1;
    }
    let rv = libc::fwrite(
        tinf as *const libc::c_void,
        std::mem::size_of::<Training>(),
        1,
        fh,
    );
    libc::fclose(fh);
    if rv != 1 {
        return -1;
    }
    0
}

// Generate all 50 initialize_metagenome_N functions
macro_rules! init_metagenome_fn {
    ($n:literal, $name:ident) => {
        #[no_mangle]
        pub unsafe extern "C" fn $name(tptr: *mut Training) {
            training_data::load_metagenome($n, tptr);
        }
    };
}

init_metagenome_fn!(0, initialize_metagenome_0);
init_metagenome_fn!(1, initialize_metagenome_1);
init_metagenome_fn!(2, initialize_metagenome_2);
init_metagenome_fn!(3, initialize_metagenome_3);
init_metagenome_fn!(4, initialize_metagenome_4);
init_metagenome_fn!(5, initialize_metagenome_5);
init_metagenome_fn!(6, initialize_metagenome_6);
init_metagenome_fn!(7, initialize_metagenome_7);
init_metagenome_fn!(8, initialize_metagenome_8);
init_metagenome_fn!(9, initialize_metagenome_9);
init_metagenome_fn!(10, initialize_metagenome_10);
init_metagenome_fn!(11, initialize_metagenome_11);
init_metagenome_fn!(12, initialize_metagenome_12);
init_metagenome_fn!(13, initialize_metagenome_13);
init_metagenome_fn!(14, initialize_metagenome_14);
init_metagenome_fn!(15, initialize_metagenome_15);
init_metagenome_fn!(16, initialize_metagenome_16);
init_metagenome_fn!(17, initialize_metagenome_17);
init_metagenome_fn!(18, initialize_metagenome_18);
init_metagenome_fn!(19, initialize_metagenome_19);
init_metagenome_fn!(20, initialize_metagenome_20);
init_metagenome_fn!(21, initialize_metagenome_21);
init_metagenome_fn!(22, initialize_metagenome_22);
init_metagenome_fn!(23, initialize_metagenome_23);
init_metagenome_fn!(24, initialize_metagenome_24);
init_metagenome_fn!(25, initialize_metagenome_25);
init_metagenome_fn!(26, initialize_metagenome_26);
init_metagenome_fn!(27, initialize_metagenome_27);
init_metagenome_fn!(28, initialize_metagenome_28);
init_metagenome_fn!(29, initialize_metagenome_29);
init_metagenome_fn!(30, initialize_metagenome_30);
init_metagenome_fn!(31, initialize_metagenome_31);
init_metagenome_fn!(32, initialize_metagenome_32);
init_metagenome_fn!(33, initialize_metagenome_33);
init_metagenome_fn!(34, initialize_metagenome_34);
init_metagenome_fn!(35, initialize_metagenome_35);
init_metagenome_fn!(36, initialize_metagenome_36);
init_metagenome_fn!(37, initialize_metagenome_37);
init_metagenome_fn!(38, initialize_metagenome_38);
init_metagenome_fn!(39, initialize_metagenome_39);
init_metagenome_fn!(40, initialize_metagenome_40);
init_metagenome_fn!(41, initialize_metagenome_41);
init_metagenome_fn!(42, initialize_metagenome_42);
init_metagenome_fn!(43, initialize_metagenome_43);
init_metagenome_fn!(44, initialize_metagenome_44);
init_metagenome_fn!(45, initialize_metagenome_45);
init_metagenome_fn!(46, initialize_metagenome_46);
init_metagenome_fn!(47, initialize_metagenome_47);
init_metagenome_fn!(48, initialize_metagenome_48);
init_metagenome_fn!(49, initialize_metagenome_49);
