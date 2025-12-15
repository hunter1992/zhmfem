use std::os::raw::c_void;

// 链接动态库
#[link(name = "pardiso", kind = "dylib")]
unsafe extern "C" {
    pub fn pardisoinit(
        pt: *mut *mut c_void,
        mtype: *mut i32,
        solvr: *mut i32,
        iparm: *mut i32,
        dparm: *mut f64,
        error: *mut i32,
    );

    pub fn pardiso(
        pt: *mut *mut c_void,
        maxfct: *mut i32,
        mnum: *mut i32,
        mtype: *mut i32,
        phase: *mut i32,
        n: *mut i32,
        a: *mut f64,
        ia: *mut i32,
        ja: *mut i32,
        perm: *mut i32,
        nrhs: *mut i32,
        iparm: *mut i32,
        msglvl: *mut i32,
        b: *mut f64,
        x: *mut f64,
        error: *mut i32,
        dparm: *mut f64,
    );
}
