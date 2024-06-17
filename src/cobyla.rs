#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut,
    clippy::needless_return,
    clippy::zero_ptr,
    clippy::toplevel_ref_arg,
    clippy::nonminimal_bool,
    clippy::assign_op_pattern,
    clippy::collapsible_if,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::single_match,
    clippy::unnecessary_cast
)]

use std::convert::TryFrom;

#[repr(C)]
pub(crate) enum CobylaStatus {
    COBYLA_INITIAL_ITERATE = 2,
    COBYLA_ITERATE = 1,
    COBYLA_SUCCESS = 0,
    COBYLA_ROUNDING_ERRORS = -1,
    COBYLA_TOO_MANY_EVALUATIONS = -2,
    COBYLA_BAD_ADDRESS = -3,
    COBYLA_CORRUPTED = -4,
}

//  {
//     // pub type _IO_wide_data;
//     // pub type _IO_codecvt;
//     // pub type _IO_marker;
//     // fn malloc(_: u64) -> *mut std::ffi::c_void;
//     // fn free(__ptr: *mut std::ffi::c_void);
//     // fn memset(_: *mut std::ffi::c_void, _: i32, _: u64) -> *mut std::ffi::c_void;
//     // fn __errno_location() -> *mut i32;
//     // fn sqrt(_: f64) -> f64;
//     // fn fabs(_: f64) -> f64;
//     // static mut stdout: *mut FILE;
//     // static mut stderr: *mut FILE;
//     // fn fprintf(_: *mut FILE, _: *const u8, _: ...) -> i32;
// }
pub(crate) type size_t = u64;
pub(crate) type __off_t = i64;
pub(crate) type __off64_t = i64;
// #[derive(Copy, Clone)]
// #[repr(C)]
// pub struct _IO_FILE {
//     pub _flags: i32,
//     pub _IO_read_ptr: *mut u8,
//     pub _IO_read_end: *mut u8,
//     pub _IO_read_base: *mut u8,
//     pub _IO_write_base: *mut u8,
//     pub _IO_write_ptr: *mut u8,
//     pub _IO_write_end: *mut u8,
//     pub _IO_buf_base: *mut u8,
//     pub _IO_buf_end: *mut u8,
//     pub _IO_save_base: *mut u8,
//     pub _IO_backup_base: *mut u8,
//     pub _IO_save_end: *mut u8,
//     pub _markers: *mut _IO_marker,
//     pub _chain: *mut _IO_FILE,
//     pub _fileno: i32,
//     pub _flags2: i32,
//     pub _old_offset: __off_t,
//     pub _cur_column: libc::c_ushort,
//     pub _vtable_offset: libc::c_schar,
//     pub _shortbuf: [u8; 1],
//     pub _lock: *mut std::ffi::c_void,
//     pub _offset: __off64_t,
//     pub _codecvt: *mut _IO_codecvt,
//     pub _wide_data: *mut _IO_wide_data,
//     pub _freeres_list: *mut _IO_FILE,
//     pub _freeres_buf: *mut std::ffi::c_void,
//     pub __pad5: size_t,
//     pub _mode: i32,
//     pub _unused2: [u8; 20],
// }
// pub type _IO_lock_t = ();
// pub type FILE = _IO_FILE;
pub(crate) type cobyla_calcfc = unsafe fn(
    i64,
    i64,
    *const f64,
    *mut f64,
    *mut std::ffi::c_void,
) -> f64;
#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct _cobyla_context {
    pub n: i64,
    pub m: i64,
    pub iprint: i64,
    pub maxfun: i64,
    pub nfvals: i64,
    pub rhobeg: f64,
    pub rhoend: f64,
    pub iact: *mut i64,
    pub con: *mut f64,
    pub sim: *mut f64,
    pub simi: *mut f64,
    pub datmat: *mut f64,
    pub a: *mut f64,
    pub vsig: *mut f64,
    pub veta: *mut f64,
    pub sigbar: *mut f64,
    pub dx: *mut f64,
    pub w: *mut f64,
    pub parmu: f64,
    pub parsig: f64,
    pub prerec: f64,
    pub prerem: f64,
    pub rho: f64,
    pub f: f64,
    pub ibrnch: i64,
    pub iflag: i64,
    pub ifull: i64,
    pub jdrop: i64,
    pub status: i32,
}

pub type cobyla_context_t = _cobyla_context;

impl Default for cobyla_context_t {
    fn default() -> Self {
        cobyla_context_t {
            n: 0,
            m: 0,
            iprint: 0,
            maxfun: 0,
            nfvals: 0,
            rhobeg: 0.,
            rhoend: 0.,
            iact: 0 as *mut i64,
            con: 0 as *mut f64,
            sim: 0 as *mut f64,
            simi: 0 as *mut f64,
            datmat: 0 as *mut f64,
            a: 0 as *mut f64,
            vsig: 0 as *mut f64,
            veta: 0 as *mut f64,
            sigbar: 0 as *mut f64,
            dx: 0 as *mut f64,
            w: 0 as *mut f64,
            parmu: 0.,
            parsig: 0.,
            prerec: 0.,
            prerem: 0.,
            rho: 0.,
            f: 0.,
            ibrnch: 0,
            iflag: 0,
            ifull: 0,
            jdrop: 0,
            status: 0,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn raw_cobyla(
    mut n: i64,
    mut m: i64,
    mut calcfc: Option<cobyla_calcfc>,
    mut calcfc_data: *mut std::ffi::c_void,
    mut x: *mut f64,
    mut rhobeg: f64,
    mut rhoend: f64,
    mut iprint: i64,
    mut maxfun: *mut i64,
    mut w: *mut f64,
    mut iact: *mut i64,
) -> i32 {
    let mut mpp: i64 = m + 2 as i32 as i64;
    let mut con: *mut f64 = w;
    let mut sim: *mut f64 = con.offset(mpp as isize);
    let mut simi: *mut f64 = sim.offset((n * n) as isize).offset(n as isize);
    let mut datmat: *mut f64 = simi.offset((n * n) as isize);
    let mut a: *mut f64 = datmat.offset((n * mpp) as isize).offset(mpp as isize);
    let mut vsig: *mut f64 = a.offset((m * n) as isize).offset(n as isize);
    let mut veta: *mut f64 = vsig.offset(n as isize);
    let mut sigbar: *mut f64 = veta.offset(n as isize);
    let mut dx: *mut f64 = sigbar.offset(n as isize);
    let mut work: *mut f64 = dx.offset(n as isize);
    return cobylb(
        n,
        m,
        calcfc,
        calcfc_data,
        x,
        rhobeg,
        rhoend,
        iprint,
        maxfun,
        con,
        sim,
        simi,
        datmat,
        a,
        vsig,
        veta,
        sigbar,
        dx,
        work,
        iact,
    );
}

pub(crate) unsafe fn cobyla_create(
    mut n: i64,
    mut m: i64,
    mut rhobeg: f64,
    mut rhoend: f64,
    mut iprint: i64,
    mut maxfun: i64,
) -> *mut cobyla_context_t {
    let mut ctx: *mut cobyla_context_t = 0 as *mut cobyla_context_t;
    let mut size: i64 = 0;
    let mut offset1: i64 = 0;
    let mut offset2: i64 = 0;
    let mut mpp: i64 = 0;
    if n < 1 as i32 as i64
        || m < 0 as i32 as i64
        || rhobeg < rhoend
        || rhoend <= 0 as i32 as f64
        || maxfun < 1 as i32 as i64
    {
        // *__errno_location() = 22 as i32;
        return 0 as *mut cobyla_context_t;
    }
    size = ::std::mem::size_of::<cobyla_context_t>() as u64 as i64;
    offset1 = (::std::mem::size_of::<i64>() as u64)
        .wrapping_sub(1 as i32 as u64)
        .wrapping_add(size as u64)
        .wrapping_div(::std::mem::size_of::<i64>() as u64)
        .wrapping_mul(::std::mem::size_of::<i64>() as u64)
        as i64;
    size = (offset1 as u64).wrapping_add(
        ((m + 1 as i32 as i64) as u64)
            .wrapping_mul(::std::mem::size_of::<i64>() as u64),
    ) as i64;
    offset2 = (::std::mem::size_of::<f64>() as u64)
        .wrapping_sub(1 as i32 as u64)
        .wrapping_add(size as u64)
        .wrapping_div(::std::mem::size_of::<f64>() as u64)
        .wrapping_mul(::std::mem::size_of::<f64>() as u64)
        as i64;
    size = (offset2 as u64).wrapping_add(
        ((n * (3 as i32 as i64 * n
            + 2 as i32 as i64 * m
            + 11 as i32 as i64)
            + 4 as i32 as i64 * m
            + 6 as i32 as i64) as u64)
            .wrapping_mul(::std::mem::size_of::<f64>() as u64),
    ) as i64;
    // ctx = malloc(size as u64) as *mut cobyla_context_t;
    let mut vec: Box<Vec<i32>> = Box::new(vec![0; usize::try_from(size).unwrap()]);
    ctx = vec.as_mut_ptr() as *mut cobyla_context_t;
    std::mem::forget(vec);

    if ctx.is_null() {
        return 0 as *mut cobyla_context_t;
    }
    // memset(
    //     ctx as *mut std::ffi::c_void,
    //     0 as i32,
    //     size as u64,
    // );
    (*ctx).n = n;
    (*ctx).m = m;
    (*ctx).nfvals = 0 as i32 as i64;
    (*ctx).status = 1 as i32;
    (*ctx).iprint = iprint;
    (*ctx).maxfun = maxfun;
    (*ctx).rhobeg = rhobeg;
    (*ctx).rhoend = rhoend;
    mpp = m + 2 as i32 as i64;
    let ref mut fresh0 = (*ctx).iact;
    *fresh0 = (ctx as *mut u8).offset(offset1 as isize) as *mut i64;
    let ref mut fresh1 = (*ctx).con;
    *fresh1 = (ctx as *mut u8).offset(offset2 as isize) as *mut f64;
    let ref mut fresh2 = (*ctx).sim;
    *fresh2 = ((*ctx).con).offset(mpp as isize);
    let ref mut fresh3 = (*ctx).simi;
    *fresh3 = ((*ctx).sim).offset((n * n) as isize).offset(n as isize);
    let ref mut fresh4 = (*ctx).datmat;
    *fresh4 = ((*ctx).simi).offset((n * n) as isize);
    let ref mut fresh5 = (*ctx).a;
    *fresh5 = ((*ctx).datmat)
        .offset((n * mpp) as isize)
        .offset(mpp as isize);
    let ref mut fresh6 = (*ctx).vsig;
    *fresh6 = ((*ctx).a).offset((m * n) as isize).offset(n as isize);
    let ref mut fresh7 = (*ctx).veta;
    *fresh7 = ((*ctx).vsig).offset(n as isize);
    let ref mut fresh8 = (*ctx).sigbar;
    *fresh8 = ((*ctx).veta).offset(n as isize);
    let ref mut fresh9 = (*ctx).dx;
    *fresh9 = ((*ctx).sigbar).offset(n as isize);
    let ref mut fresh10 = (*ctx).w;
    *fresh10 = ((*ctx).dx).offset(n as isize);
    return ctx;
}

pub(crate) unsafe fn cobyla_delete(mut ctx: *mut cobyla_context_t) {
    if !ctx.is_null() {
        let _ = Box::from_raw(ctx);
    }
}

pub(crate) unsafe fn cobyla_restart(mut ctx: *mut cobyla_context_t) -> i32 {
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(3 as i32);
    }
    (*ctx).nfvals = 0 as i32 as i64;
    (*ctx).status = 1 as i32;
    return (*ctx).status;
}

pub(crate) unsafe fn cobyla_get_status(mut ctx: *const cobyla_context_t) -> i32 {
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(3 as i32);
    }
    return (*ctx).status;
}

pub(crate) unsafe fn cobyla_get_nevals(mut ctx: *const cobyla_context_t) -> i64 {
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(1 as i32) as i64;
    }
    return (*ctx).nfvals;
}

pub(crate) unsafe fn cobyla_get_rho(mut ctx: *const cobyla_context_t) -> f64 {
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(1 as i32) as f64;
    }
    return (*ctx).rho;
}

pub(crate) unsafe fn cobyla_get_last_f(mut ctx: *const cobyla_context_t) -> f64 {
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(1 as i32) as f64;
    }
    return (*ctx).f;
}

pub(crate) unsafe fn cobyla_iterate(
    mut ctx: *mut cobyla_context_t,
    mut f: f64,
    mut x: *mut f64,
    mut c: *mut f64,
) -> i32 {
    let mut current_block: u64;
    let zero: f64 = 0.0f64;
    let one: f64 = 1.0f64;
    let alpha: f64 = 0.25f64;
    let beta: f64 = 2.1f64;
    let gamma: f64 = 0.5f64;
    let delta: f64 = 1.1f64;
    let mut parmu: f64 = 0.;
    let mut parsig: f64 = 0.;
    let mut prerec: f64 = 0.;
    let mut prerem: f64 = 0.;
    let mut rho: f64 = 0.;
    let mut barmu: f64 = 0.;
    let mut cvmaxm: f64 = 0.;
    let mut cvmaxp: f64 = 0.;
    let mut dxsign: f64 = 0.;
    let mut edgmax: f64 = 0.;
    let mut error: f64 = 0.;
    let mut pareta: f64 = 0.;
    let mut phi: f64 = 0.;
    let mut phimin: f64 = 0.;
    let mut ratio: f64 = 0.;
    let mut resmax: f64 = 0.;
    let mut resnew: f64 = 0.;
    let mut sum: f64 = 0.;
    let mut temp: f64 = 0.;
    let mut tempa: f64 = 0.;
    let mut trured: f64 = 0.;
    let mut vmnew: f64 = 0.;
    let mut vmold: f64 = 0.;
    let mut ibrnch: i64 = 0;
    let mut iflag: i64 = 0;
    let mut ifull: i64 = 0;
    let mut jdrop: i64 = 0;
    let mut nfvals: i64 = 0;
    let mut maxfun: i64 = 0;
    let mut i: i64 = 0;
    let mut j: i64 = 0;
    let mut k: i64 = 0;
    let mut l: i64 = 0;
    let mut mp: i64 = 0;
    let mut mpp: i64 = 0;
    let mut np: i64 = 0;
    let mut nbest: i64 = 0;
    let mut n: i64 = 0;
    let mut m: i64 = 0;
    let mut iprint: i64 = 0;
    let mut rhobeg: f64 = 0.;
    let mut rhoend: f64 = 0.;
    let mut con: *mut f64 = 0 as *mut f64;
    let mut sim: *mut f64 = 0 as *mut f64;
    let mut simi: *mut f64 = 0 as *mut f64;
    let mut datmat: *mut f64 = 0 as *mut f64;
    let mut a: *mut f64 = 0 as *mut f64;
    let mut vsig: *mut f64 = 0 as *mut f64;
    let mut veta: *mut f64 = 0 as *mut f64;
    let mut sigbar: *mut f64 = 0 as *mut f64;
    let mut dx: *mut f64 = 0 as *mut f64;
    let mut w: *mut f64 = 0 as *mut f64;
    let mut iact: *mut i64 = 0 as *mut i64;
    let mut status: i32 = 0;
    if ctx.is_null() {
        // *__errno_location() = 14 as i32;
        return -(3 as i32);
    }
    n = (*ctx).n;
    m = (*ctx).m;
    iprint = (*ctx).iprint;
    maxfun = (*ctx).maxfun;
    nfvals = (*ctx).nfvals;
    rhobeg = (*ctx).rhobeg;
    rhoend = (*ctx).rhoend;
    iact = (*ctx).iact;
    con = (*ctx).con;
    sim = (*ctx).sim;
    simi = (*ctx).simi;
    datmat = (*ctx).datmat;
    a = (*ctx).a;
    vsig = (*ctx).vsig;
    veta = (*ctx).veta;
    sigbar = (*ctx).sigbar;
    dx = (*ctx).dx;
    w = (*ctx).w;
    status = (*ctx).status;
    if x.is_null() || c.is_null() && m > 0 as i32 as i64 {
        // *__errno_location() = 14 as i32;
        (*ctx).status = -(3 as i32);
        return -(3 as i32);
    }
    if status != 1 as i32 || nfvals < 0 as i32 as i64 {
        // *__errno_location() = 22 as i32;
        (*ctx).status = -(4 as i32);
        return (*ctx).status;
    }
    np = n + 1 as i32 as i64;
    mp = m + 1 as i32 as i64;
    mpp = m + 2 as i32 as i64;
    if nfvals == 0 as i32 as i64 {
        status = 2 as i32;
        prerec = zero;
        prerem = zero;
        parsig = zero;
        iflag = 0 as i32 as i64;
        rho = rhobeg;
        parmu = zero;
        jdrop = np;
        ibrnch = 0 as i32 as i64;
        temp = one / rho;
        i = 1 as i32 as i64;
        while i <= n {
            *sim.offset(
                (i - 1 as i32 as i64 + n * (np - 1 as i32 as i64))
                    as isize,
            ) = *x.offset((i - 1 as i32 as i64) as isize);
            j = 1 as i32 as i64;
            while j <= n {
                *sim.offset(
                    (i - 1 as i32 as i64
                        + n * (j - 1 as i32 as i64)) as isize,
                ) = zero;
                *simi.offset(
                    (i - 1 as i32 as i64
                        + n * (j - 1 as i32 as i64)) as isize,
                ) = zero;
                j += 1;
            }
            *sim.offset(
                (i - 1 as i32 as i64 + n * (i - 1 as i32 as i64))
                    as isize,
            ) = rho;
            *simi.offset(
                (i - 1 as i32 as i64 + n * (i - 1 as i32 as i64))
                    as isize,
            ) = temp;
            i += 1;
        }
        if iprint >= 2 as i32 as i64 {
            println!(
                "\n   The initial value of RHO is {}  and PARMU is set to zero.\n",
                rho,
            );
        }
        current_block = 14453151562619017203;
    } else {
        parmu = (*ctx).parmu;
        parsig = (*ctx).parsig;
        prerec = (*ctx).prerec;
        prerem = (*ctx).prerem;
        rho = (*ctx).rho;
        ibrnch = (*ctx).ibrnch;
        iflag = (*ctx).iflag;
        ifull = (*ctx).ifull;
        jdrop = (*ctx).jdrop;
        current_block = 2561692241959298350;
    }
    'c_12387: loop {
        match current_block {
            14453151562619017203 => {
                if nfvals >= maxfun && nfvals > 0 as i32 as i64 {
                    status = -(2 as i32);
                    if iprint > 0 as i32 as i64 {
                        println!(
                            "Return from subroutine COBYLA because {:?}.\n",
                            cobyla_reason(status),
                        );
                    }
                    current_block = 2880604979707412946;
                    break;
                } else if status == 2 as i32 {
                    status = 1 as i32;
                    current_block = 2561692241959298350;
                } else {
                    status = 1 as i32;
                    current_block = 14015749458106396694;
                    break;
                }
            }
            _ => {
                nfvals += 1;
                resmax = zero;
                k = 0 as i32 as i64;
                while k < m {
                    temp = *c.offset(k as isize);
                    *con.offset(k as isize) = temp;
                    if resmax < -temp {
                        resmax = -temp;
                    }
                    k += 1;
                }
                if nfvals == iprint - 1 as i32 as i64
                    || iprint == 3 as i32 as i64
                {
                    print_calcfc(
                        n,
                        nfvals,
                        f,
                        resmax,
                        &mut *x.offset(0) as *mut f64 as *const f64,
                    );
                }
                *con.offset((mp - 1 as i32 as i64) as isize) = f;
                *con.offset((mpp - 1 as i32 as i64) as isize) = resmax;
                if ibrnch == 1 as i32 as i64 {
                    vmold = *datmat.offset(
                        (mp - 1 as i32 as i64
                            + mpp * (np - 1 as i32 as i64))
                            as isize,
                    ) + parmu
                        * *datmat.offset(
                            (mpp - 1 as i32 as i64
                                + mpp * (np - 1 as i32 as i64))
                                as isize,
                        );
                    vmnew = f + parmu * resmax;
                    trured = vmold - vmnew;
                    if parmu == zero
                        && f == *datmat.offset(
                            (mp - 1 as i32 as i64
                                + mpp * (np - 1 as i32 as i64))
                                as isize,
                        )
                    {
                        prerem = prerec;
                        trured = *datmat.offset(
                            (mpp - 1 as i32 as i64
                                + mpp * (np - 1 as i32 as i64))
                                as isize,
                        ) - resmax;
                    }
                    ratio = if trured <= zero { one } else { zero };
                    jdrop = 0 as i32 as i64;
                    j = 1 as i32 as i64;
                    while j <= n {
                        temp = zero;
                        i = 1 as i32 as i64;
                        while i <= n {
                            temp += *simi.offset(
                                (j - 1 as i32 as i64
                                    + n * (i - 1 as i32 as i64))
                                    as isize,
                            ) * *dx.offset((i - 1 as i32 as i64) as isize);
                            i += 1;
                        }
                        temp = temp.abs();
                        if temp > ratio {
                            jdrop = j;
                            ratio = temp;
                        }
                        *sigbar.offset((j - 1 as i32 as i64) as isize) =
                            temp * *vsig.offset((j - 1 as i32 as i64) as isize);
                        j += 1;
                    }
                    edgmax = delta * rho;
                    l = 0 as i32 as i64;
                    j = 1 as i32 as i64;
                    while j <= n {
                        if *sigbar.offset((j - 1 as i32 as i64) as isize) >= parsig
                            || *sigbar.offset((j - 1 as i32 as i64) as isize)
                                >= *vsig.offset((j - 1 as i32 as i64) as isize)
                        {
                            temp = *veta.offset((j - 1 as i32 as i64) as isize);
                            if trured > zero {
                                temp = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    let mut tempb: f64 = *dx
                                        .offset((i - 1 as i32 as i64) as isize)
                                        - *sim.offset(
                                            (i - 1 as i32 as i64
                                                + n * (j - 1 as i32 as i64))
                                                as isize,
                                        );
                                    temp += tempb * tempb;
                                    i += 1;
                                }
                                temp = temp.sqrt();
                            }
                            if temp > edgmax {
                                l = j;
                                edgmax = temp;
                            }
                        }
                        j += 1;
                    }
                    if l > 0 as i32 as i64 {
                        jdrop = l;
                    }
                    if jdrop == 0 as i32 as i64 {
                        current_block = 12414752556692412193;
                    } else {
                        temp = zero;
                        i = 1 as i32 as i64;
                        while i <= n {
                            *sim.offset(
                                (i - 1 as i32 as i64
                                    + n * (jdrop - 1 as i32 as i64))
                                    as isize,
                            ) = *dx.offset((i - 1 as i32 as i64) as isize);
                            temp += *simi.offset(
                                (jdrop - 1 as i32 as i64
                                    + n * (i - 1 as i32 as i64))
                                    as isize,
                            ) * *dx.offset((i - 1 as i32 as i64) as isize);
                            i += 1;
                        }
                        i = 1 as i32 as i64;
                        while i <= n {
                            *simi.offset(
                                (jdrop - 1 as i32 as i64
                                    + n * (i - 1 as i32 as i64))
                                    as isize,
                            ) = *simi.offset(
                                (jdrop - 1 as i32 as i64
                                    + n * (i - 1 as i32 as i64))
                                    as isize,
                            ) / temp;
                            i += 1;
                        }
                        j = 1 as i32 as i64;
                        while j <= n {
                            if j != jdrop {
                                temp = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp += *simi.offset(
                                        (j - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) * *dx
                                        .offset((i - 1 as i32 as i64) as isize);
                                    i += 1;
                                }
                                i = 1 as i32 as i64;
                                while i <= n {
                                    *simi.offset(
                                        (j - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) -= temp
                                        * *simi.offset(
                                            (jdrop - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                    i += 1;
                                }
                            }
                            j += 1;
                        }
                        k = 1 as i32 as i64;
                        while k <= mpp {
                            *datmat.offset(
                                (k - 1 as i32 as i64
                                    + mpp * (jdrop - 1 as i32 as i64))
                                    as isize,
                            ) = *con.offset((k - 1 as i32 as i64) as isize);
                            k += 1;
                        }
                        if trured > zero && trured >= 0.1f64 * prerem {
                            current_block = 10213206415157612706;
                        } else {
                            current_block = 12414752556692412193;
                        }
                    }
                } else {
                    k = 1 as i32 as i64;
                    while k <= mpp {
                        *datmat.offset(
                            (k - 1 as i32 as i64
                                + mpp * (jdrop - 1 as i32 as i64))
                                as isize,
                        ) = *con.offset((k - 1 as i32 as i64) as isize);
                        k += 1;
                    }
                    if !(nfvals > np) {
                        if jdrop <= n {
                            if *datmat.offset(
                                (mp - 1 as i32 as i64
                                    + mpp * (np - 1 as i32 as i64))
                                    as isize,
                            ) <= f
                            {
                                *x.offset((jdrop - 1 as i32 as i64) as isize) =
                                    *sim.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                            } else {
                                *sim.offset(
                                    (jdrop - 1 as i32 as i64
                                        + n * (np - 1 as i32 as i64))
                                        as isize,
                                ) = *x.offset((jdrop - 1 as i32 as i64) as isize);
                                k = 1 as i32 as i64;
                                while k <= mpp {
                                    *datmat.offset(
                                        (k - 1 as i32 as i64
                                            + mpp * (jdrop - 1 as i32 as i64))
                                            as isize,
                                    ) = *datmat.offset(
                                        (k - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (k - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    ) = *con
                                        .offset((k - 1 as i32 as i64) as isize);
                                    k += 1;
                                }
                                k = 1 as i32 as i64;
                                while k <= jdrop {
                                    temp = zero;
                                    *sim.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) = -rho;
                                    i = k;
                                    while i <= jdrop {
                                        temp -= *simi.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    *simi.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) = temp;
                                    k += 1;
                                }
                            }
                        }
                        if nfvals <= n {
                            jdrop = nfvals;
                            *x.offset((jdrop - 1 as i32 as i64) as isize) += rho;
                            current_block = 14453151562619017203;
                            continue;
                        }
                    }
                    ibrnch = 1 as i32 as i64;
                    current_block = 10213206415157612706;
                }
                'c_12399: loop {
                    match current_block {
                        12414752556692412193 => {
                            if iflag == 0 as i32 as i64 {
                                ibrnch = 0 as i32 as i64;
                                current_block = 10213206415157612706;
                            } else if rho > rhoend {
                                rho = 0.5f64 * rho;
                                if rho <= 1.5f64 * rhoend {
                                    rho = rhoend;
                                }
                                if parmu > zero {
                                    let mut cmin: f64 = 0.;
                                    let mut cmax: f64 = 0.;
                                    let mut denom: f64 = 0.;
                                    denom = zero;
                                    k = 1 as i32 as i64;
                                    while k <= mp {
                                        cmin = *datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        );
                                        cmax = cmin;
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            temp = *datmat.offset(
                                                (k - 1 as i32 as i64
                                                    + mpp * (i - 1 as i32 as i64))
                                                    as isize,
                                            );
                                            cmin = if cmin <= temp { cmin } else { temp };
                                            cmax = if cmax >= temp { cmax } else { temp };
                                            i += 1;
                                        }
                                        if k <= m && cmin < 0.5f64 * cmax {
                                            temp = (if cmax >= zero { cmax } else { zero }) - cmin;
                                            if denom <= zero {
                                                denom = temp;
                                            } else {
                                                denom = if denom <= temp { denom } else { temp };
                                            }
                                        }
                                        k += 1;
                                    }
                                    if denom == zero {
                                        parmu = zero;
                                    } else if cmax - cmin < parmu * denom {
                                        parmu = (cmax - cmin) / denom;
                                    }
                                }
                                if iprint >= 2 as i32 as i64 {
                                    println!(
                                        "\n   Reduction in RHO to {}  and PARMU ={}\n",
                                        rho, parmu,
                                    );
                                    if iprint == 2 as i32 as i64 {
                                        print_calcfc(
                                            n,
                                            nfvals,
                                            *datmat.offset(
                                                (mp - 1 as i32 as i64
                                                    + mpp * (np - 1 as i32 as i64))
                                                    as isize,
                                            ),
                                            *datmat.offset(
                                                (mpp - 1 as i32 as i64
                                                    + mpp * (np - 1 as i32 as i64))
                                                    as isize,
                                            ),
                                            &mut *sim.offset(
                                                (n * (np - 1 as i32 as i64))
                                                    as isize,
                                            )
                                                as *mut f64
                                                as *const f64,
                                        );
                                    }
                                }
                                current_block = 10213206415157612706;
                            } else {
                                if iprint >= 1 as i32 as i64 {
                                    println!("\n   Normal return from subroutine COBYLA\n");
                                }
                                status = 0 as i32;
                                if ifull == 1 as i32 as i64 {
                                    current_block = 18000475564953935197;
                                    break 'c_12387;
                                } else {
                                    current_block = 2880604979707412946;
                                    break 'c_12387;
                                }
                            }
                        }
                        _ => {
                            phimin = *datmat.offset(
                                (mp - 1 as i32 as i64
                                    + mpp * (np - 1 as i32 as i64))
                                    as isize,
                            ) + parmu
                                * *datmat.offset(
                                    (mpp - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                );
                            nbest = np;
                            j = 1 as i32 as i64;
                            while j <= n {
                                temp = *datmat.offset(
                                    (mp - 1 as i32 as i64
                                        + mpp * (j - 1 as i32 as i64))
                                        as isize,
                                ) + parmu
                                    * *datmat.offset(
                                        (mpp - 1 as i32 as i64
                                            + mpp * (j - 1 as i32 as i64))
                                            as isize,
                                    );
                                if temp < phimin {
                                    nbest = j;
                                    phimin = temp;
                                } else if temp == phimin && parmu == zero {
                                    if *datmat.offset(
                                        (mpp - 1 as i32 as i64
                                            + mpp * (j - 1 as i32 as i64))
                                            as isize,
                                    ) < *datmat.offset(
                                        (mpp - 1 as i32 as i64
                                            + mpp * (nbest - 1 as i32 as i64))
                                            as isize,
                                    ) {
                                        nbest = j;
                                    }
                                }
                                j += 1;
                            }
                            if nbest <= n {
                                i = 1 as i32 as i64;
                                while i <= mpp {
                                    temp = *datmat.offset(
                                        (i - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (i - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    ) = *datmat.offset(
                                        (i - 1 as i32 as i64
                                            + mpp * (nbest - 1 as i32 as i64))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (i - 1 as i32 as i64
                                            + mpp * (nbest - 1 as i32 as i64))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp = *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nbest - 1 as i32 as i64))
                                            as isize,
                                    );
                                    *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nbest - 1 as i32 as i64))
                                            as isize,
                                    ) = zero;
                                    *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (np - 1 as i32 as i64))
                                            as isize,
                                    ) += temp;
                                    tempa = zero;
                                    k = 1 as i32 as i64;
                                    while k <= n {
                                        *sim.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) -= temp;
                                        tempa -= *simi.offset(
                                            (k - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                        k += 1;
                                    }
                                    *simi.offset(
                                        (nbest - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) = tempa;
                                    i += 1;
                                }
                            }
                            error = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                j = 1 as i32 as i64;
                                while j <= n {
                                    temp = if i == j { -one } else { zero };
                                    k = 1 as i32 as i64;
                                    while k <= n {
                                        temp += *simi.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) * *sim.offset(
                                            (k - 1 as i32 as i64
                                                + n * (j - 1 as i32 as i64))
                                                as isize,
                                        );
                                        k += 1;
                                    }
                                    temp = temp.abs();
                                    error = if error >= temp { error } else { temp };
                                    j += 1;
                                }
                                i += 1;
                            }
                            if error > 0.1f64 {
                                status = -(1 as i32);
                                if iprint >= 1 as i32 as i64 {
                                    println!(
                                        "Return from subroutine COBYLA because {:?}.\n",
                                        cobyla_reason(status),
                                    );
                                }
                                current_block = 2880604979707412946;
                                break 'c_12387;
                            } else {
                                k = 1 as i32 as i64;
                                while k <= mp {
                                    *con.offset((k - 1 as i32 as i64) as isize) =
                                        -*datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        );
                                    j = 1 as i32 as i64;
                                    while j <= n {
                                        *w.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        ) = *datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (j - 1 as i32 as i64))
                                                as isize,
                                        ) + *con.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                        j += 1;
                                    }
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        temp = zero;
                                        j = 1 as i32 as i64;
                                        while j <= n {
                                            temp += *w.offset(
                                                (j - 1 as i32 as i64) as isize,
                                            ) * *simi.offset(
                                                (j - 1 as i32 as i64
                                                    + n * (i - 1 as i32 as i64))
                                                    as isize,
                                            );
                                            j += 1;
                                        }
                                        *a.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) = if k == mp { -temp } else { temp };
                                        i += 1;
                                    }
                                    k += 1;
                                }
                                iflag = 1 as i32 as i64;
                                parsig = alpha * rho;
                                pareta = beta * rho;
                                j = 1 as i32 as i64;
                                while j <= n {
                                    let mut wsig: f64 = zero;
                                    let mut weta: f64 = zero;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        wsig += *simi.offset(
                                            (j - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        ) * *simi.offset(
                                            (j - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                        weta += *sim.offset(
                                            (i - 1 as i32 as i64
                                                + n * (j - 1 as i32 as i64))
                                                as isize,
                                        ) * *sim.offset(
                                            (i - 1 as i32 as i64
                                                + n * (j - 1 as i32 as i64))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    *vsig.offset((j - 1 as i32 as i64) as isize) =
                                        one / wsig.sqrt();
                                    *veta.offset((j - 1 as i32 as i64) as isize) =
                                        weta.sqrt();
                                    if *vsig.offset((j - 1 as i32 as i64) as isize)
                                        < parsig
                                        || *veta
                                            .offset((j - 1 as i32 as i64) as isize)
                                            > pareta
                                    {
                                        iflag = 0 as i32 as i64;
                                    }
                                    j += 1;
                                }
                                if ibrnch == 1 as i32 as i64
                                    || iflag == 1 as i32 as i64
                                {
                                    let mut z: *mut f64 = w;
                                    let mut zdota: *mut f64 = z.offset((n * n) as isize);
                                    let mut vmc: *mut f64 = zdota.offset(n as isize);
                                    let mut sdirn: *mut f64 = vmc.offset(mp as isize);
                                    let mut dxnew: *mut f64 = sdirn.offset(n as isize);
                                    let mut vmd: *mut f64 = dxnew.offset(n as isize);
                                    trstlp(
                                        n,
                                        m,
                                        a as *const f64,
                                        con as *const f64,
                                        rho,
                                        dx,
                                        &mut ifull,
                                        iact,
                                        z,
                                        zdota,
                                        vmc,
                                        sdirn,
                                        dxnew,
                                        vmd,
                                    );
                                    if ifull == 0 as i32 as i64 {
                                        temp = zero;
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            temp += *dx.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            ) * *dx.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            );
                                            i += 1;
                                        }
                                        if temp < 0.25f64 * rho * rho {
                                            ibrnch = 1 as i32 as i64;
                                            current_block = 12414752556692412193;
                                            continue;
                                        }
                                    }
                                    resnew = zero;
                                    *con.offset((mp - 1 as i32 as i64) as isize) =
                                        zero;
                                    sum = zero;
                                    k = 1 as i32 as i64;
                                    while k <= mp {
                                        sum = *con.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            sum -= *a.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (k - 1 as i32 as i64))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            );
                                            i += 1;
                                        }
                                        if k < mp {
                                            resnew = if resnew >= sum { resnew } else { sum };
                                        }
                                        k += 1;
                                    }
                                    prerec = *datmat.offset(
                                        (mpp - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    ) - resnew;
                                    barmu = if prerec > zero { sum / prerec } else { zero };
                                    if !(parmu < 1.5f64 * barmu) {
                                        break;
                                    }
                                    parmu = 2.0f64 * barmu;
                                    if iprint >= 2 as i32 as i64 {
                                        println!("\n   Increase in PARMU to {}\n", parmu);
                                    }
                                    phi = *datmat.offset(
                                        (mp - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    ) + parmu
                                        * *datmat.offset(
                                            (mpp - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        );
                                    j = 1 as i32 as i64;
                                    loop {
                                        if !(j <= n) {
                                            break 'c_12399;
                                        }
                                        temp = *datmat.offset(
                                            (mp - 1 as i32 as i64
                                                + mpp * (j - 1 as i32 as i64))
                                                as isize,
                                        ) + parmu
                                            * *datmat.offset(
                                                (mpp - 1 as i32 as i64
                                                    + mpp * (j - 1 as i32 as i64))
                                                    as isize,
                                            );
                                        if temp < phi {
                                            current_block = 10213206415157612706;
                                            break;
                                        }
                                        if temp == phi && parmu == zero {
                                            if *datmat.offset(
                                                (mpp - 1 as i32 as i64
                                                    + mpp * (j - 1 as i32 as i64))
                                                    as isize,
                                            ) < *datmat.offset(
                                                (mpp - 1 as i32 as i64
                                                    + mpp * (np - 1 as i32 as i64))
                                                    as isize,
                                            ) {
                                                current_block = 10213206415157612706;
                                                break;
                                            }
                                        }
                                        j += 1;
                                    }
                                } else {
                                    jdrop = 0 as i32 as i64;
                                    temp = pareta;
                                    j = 1 as i32 as i64;
                                    while j <= n {
                                        if *veta
                                            .offset((j - 1 as i32 as i64) as isize)
                                            > temp
                                        {
                                            jdrop = j;
                                            temp = *veta.offset(
                                                (j - 1 as i32 as i64) as isize,
                                            );
                                        }
                                        j += 1;
                                    }
                                    if jdrop == 0 as i32 as i64 {
                                        j = 1 as i32 as i64;
                                        while j <= n {
                                            if *vsig.offset(
                                                (j - 1 as i32 as i64) as isize,
                                            ) < temp
                                            {
                                                jdrop = j;
                                                temp = *vsig.offset(
                                                    (j - 1 as i32 as i64) as isize,
                                                );
                                            }
                                            j += 1;
                                        }
                                    }
                                    temp = gamma
                                        * rho
                                        * *vsig.offset(
                                            (jdrop - 1 as i32 as i64) as isize,
                                        );
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        ) = temp
                                            * *simi.offset(
                                                (jdrop - 1 as i32 as i64
                                                    + n * (i - 1 as i32 as i64))
                                                    as isize,
                                            );
                                        i += 1;
                                    }
                                    cvmaxp = zero;
                                    cvmaxm = zero;
                                    sum = zero;
                                    k = 1 as i32 as i64;
                                    while k <= mp {
                                        sum = zero;
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            sum = sum
                                                + *a.offset(
                                                    (i - 1 as i32 as i64
                                                        + n * (k - 1 as i32
                                                            as i64))
                                                        as isize,
                                                ) * *dx.offset(
                                                    (i - 1 as i32 as i64) as isize,
                                                );
                                            i += 1;
                                        }
                                        if k < mp {
                                            temp = *datmat.offset(
                                                (k - 1 as i32 as i64
                                                    + mpp * (np - 1 as i32 as i64))
                                                    as isize,
                                            );
                                            cvmaxp = if cvmaxp >= -sum - temp {
                                                cvmaxp
                                            } else {
                                                -sum - temp
                                            };
                                            cvmaxm = if cvmaxm >= sum - temp {
                                                cvmaxm
                                            } else {
                                                sum - temp
                                            };
                                        }
                                        k += 1;
                                    }
                                    if parmu * (cvmaxp - cvmaxm) > sum + sum {
                                        dxsign = -one;
                                    } else {
                                        dxsign = one;
                                    }
                                    temp = zero;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        ) *= dxsign;
                                        *sim.offset(
                                            (i - 1 as i32 as i64
                                                + n * (jdrop - 1 as i32 as i64))
                                                as isize,
                                        ) = *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        );
                                        temp += *simi.offset(
                                            (jdrop - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        ) * *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        );
                                        i += 1;
                                    }
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        *simi.offset(
                                            (jdrop - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        ) /= temp;
                                        i += 1;
                                    }
                                    j = 1 as i32 as i64;
                                    while j <= n {
                                        if j != jdrop {
                                            temp = zero;
                                            i = 1 as i32 as i64;
                                            while i <= n {
                                                temp += *simi.offset(
                                                    (j - 1 as i32 as i64
                                                        + n * (i - 1 as i32
                                                            as i64))
                                                        as isize,
                                                ) * *dx.offset(
                                                    (i - 1 as i32 as i64) as isize,
                                                );
                                                i += 1;
                                            }
                                            i = 1 as i32 as i64;
                                            while i <= n {
                                                *simi.offset(
                                                    (j - 1 as i32 as i64
                                                        + n * (i - 1 as i32
                                                            as i64))
                                                        as isize,
                                                ) -= temp
                                                    * *simi.offset(
                                                        (jdrop - 1 as i32 as i64
                                                            + n * (i - 1 as i32
                                                                as i64))
                                                            as isize,
                                                    );
                                                i += 1;
                                            }
                                        }
                                        *x.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        ) = *sim.offset(
                                            (j - 1 as i32 as i64
                                                + n * (np - 1 as i32 as i64))
                                                as isize,
                                        ) + *dx.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        );
                                        j += 1;
                                    }
                                    current_block = 14453151562619017203;
                                    continue 'c_12387;
                                }
                            }
                        }
                    }
                }
                prerem = parmu * prerec - sum;
                i = 1 as i32 as i64;
                while i <= n {
                    *x.offset((i - 1 as i32 as i64) as isize) = *sim.offset(
                        (i - 1 as i32 as i64
                            + n * (np - 1 as i32 as i64))
                            as isize,
                    ) + *dx
                        .offset((i - 1 as i32 as i64) as isize);
                    i += 1;
                }
                ibrnch = 1 as i32 as i64;
                current_block = 14453151562619017203;
            }
        }
    }
    match current_block {
        2880604979707412946 => {
            i = 1 as i32 as i64;
            while i <= n {
                *x.offset((i - 1 as i32 as i64) as isize) = *sim.offset(
                    (i - 1 as i32 as i64
                        + n * (np - 1 as i32 as i64)) as isize,
                );
                i += 1;
            }
            f = *datmat.offset(
                (mp - 1 as i32 as i64
                    + mpp * (np - 1 as i32 as i64)) as isize,
            );
            resmax = *datmat.offset(
                (mpp - 1 as i32 as i64
                    + mpp * (np - 1 as i32 as i64)) as isize,
            );
            current_block = 18000475564953935197;
        }
        _ => {}
    }
    match current_block {
        18000475564953935197 => {
            if iprint >= 1 as i32 as i64 {
                print_calcfc(
                    n,
                    nfvals,
                    f,
                    resmax,
                    &mut *x.offset(0 as isize) as *mut f64 as *const f64,
                );
            }
        }
        _ => {}
    }
    (*ctx).nfvals = nfvals;
    (*ctx).parmu = parmu;
    (*ctx).parsig = parsig;
    (*ctx).prerec = prerec;
    (*ctx).prerem = prerem;
    (*ctx).rho = rho;
    (*ctx).f = f;
    (*ctx).ibrnch = ibrnch;
    (*ctx).iflag = iflag;
    (*ctx).ifull = ifull;
    (*ctx).jdrop = jdrop;
    (*ctx).status = status;
    return status;
}

#[allow(clippy::too_many_arguments)]
unsafe fn cobylb(
    n: i64,
    m: i64,
    mut calcfc: Option<cobyla_calcfc>,
    mut calcfc_data: *mut std::ffi::c_void,
    mut x: *mut f64,
    rhobeg: f64,
    rhoend: f64,
    iprint: i64,
    mut _maxfun: *mut i64,
    mut con: *mut f64,
    mut sim: *mut f64,
    mut simi: *mut f64,
    mut datmat: *mut f64,
    mut a: *mut f64,
    mut vsig: *mut f64,
    mut veta: *mut f64,
    mut sigbar: *mut f64,
    mut dx: *mut f64,
    mut w: *mut f64,
    mut iact: *mut i64,
) -> i32 {
    let mut current_block: u64;
    let zero: f64 = 0.0f64;
    let one: f64 = 1.0f64;
    let alpha: f64 = 0.25f64;
    let beta: f64 = 2.1f64;
    let gamma: f64 = 0.5f64;
    let delta: f64 = 1.1f64;
    let mut parmu: f64 = 0.;
    let mut parsig: f64 = 0.;
    let mut prerec: f64 = 0.;
    let mut prerem: f64 = 0.;
    let mut rho: f64 = 0.;
    let mut barmu: f64 = 0.;
    let mut cvmaxm: f64 = 0.;
    let mut cvmaxp: f64 = 0.;
    let mut dxsign: f64 = 0.;
    let mut edgmax: f64 = 0.;
    let mut error: f64 = 0.;
    let mut pareta: f64 = 0.;
    let mut phi: f64 = 0.;
    let mut phimin: f64 = 0.;
    let mut ratio: f64 = 0.;
    let mut resmax: f64 = 0.;
    let mut resnew: f64 = 0.;
    let mut sum: f64 = 0.;
    let mut temp: f64 = 0.;
    let mut tempa: f64 = 0.;
    let mut trured: f64 = 0.;
    let mut vmnew: f64 = 0.;
    let mut vmold: f64 = 0.;
    let mut ibrnch: i64 = 0;
    let mut iflag: i64 = 0;
    let mut ifull: i64 = 0;
    let mut jdrop: i64 = 0;
    let mut nfvals: i64 = 0;
    let mut maxfun: i64 = 0;
    let mut i: i64 = 0;
    let mut j: i64 = 0;
    let mut k: i64 = 0;
    let mut l: i64 = 0;
    let mut mp: i64 = 0;
    let mut mpp: i64 = 0;
    let mut np: i64 = 0;
    let mut nbest: i64 = 0;
    let mut f: f64 = 0.;
    let mut status: i32 = 0;
    np = n + 1 as i32 as i64;
    mp = m + 1 as i32 as i64;
    mpp = m + 2 as i32 as i64;
    maxfun = *_maxfun;
    nfvals = 0 as i32 as i64;
    rho = rhobeg;
    parmu = zero;
    jdrop = np;
    ibrnch = 0 as i32 as i64;
    temp = one / rho;
    i = 1 as i32 as i64;
    while i <= n {
        *sim.offset(
            (i - 1 as i32 as i64 + n * (np - 1 as i32 as i64))
                as isize,
        ) = *x.offset((i - 1 as i32 as i64) as isize);
        j = 1 as i32 as i64;
        while j <= n {
            *sim.offset(
                (i - 1 as i32 as i64 + n * (j - 1 as i32 as i64))
                    as isize,
            ) = zero;
            *simi.offset(
                (i - 1 as i32 as i64 + n * (j - 1 as i32 as i64))
                    as isize,
            ) = zero;
            j += 1;
        }
        *sim.offset(
            (i - 1 as i32 as i64 + n * (i - 1 as i32 as i64))
                as isize,
        ) = rho;
        *simi.offset(
            (i - 1 as i32 as i64 + n * (i - 1 as i32 as i64))
                as isize,
        ) = temp;
        i += 1;
    }
    if iprint >= 2 as i32 as i64 {
        println!(
            "\n   The initial value of RHO is {}  and PARMU is set to zero.\n",
            rho,
        );
    }
    'c_3166: loop {
        if nfvals >= maxfun && nfvals > 0 as i32 as i64 {
            status = -(2 as i32);
            if iprint > 0 as i32 as i64 {
                println!(
                    "Return from subroutine COBYLA because {:?}.\n",
                    cobyla_reason(status),
                );
            }
            current_block = 18034400475116118263;
            break;
        } else {
            f = calcfc.expect("non-null function pointer")(
                n,
                m,
                x as *const f64,
                con,
                calcfc_data,
            );
            nfvals += 1;
            resmax = zero;
            k = 0 as i32 as i64;
            while k < m {
                temp = -*con.offset(k as isize);
                if resmax < temp {
                    resmax = temp;
                }
                k += 1;
            }
            if nfvals == iprint - 1 as i32 as i64
                || iprint == 3 as i32 as i64
            {
                print_calcfc(
                    n,
                    nfvals,
                    f,
                    resmax,
                    &mut *x.offset(0) as *mut f64 as *const f64,
                );
            }
            *con.offset((mp - 1 as i32 as i64) as isize) = f;
            *con.offset((mpp - 1 as i32 as i64) as isize) = resmax;
            if ibrnch == 1 as i32 as i64 {
                vmold = *datmat.offset(
                    (mp - 1 as i32 as i64
                        + mpp * (np - 1 as i32 as i64))
                        as isize,
                ) + parmu
                    * *datmat.offset(
                        (mpp - 1 as i32 as i64
                            + mpp * (np - 1 as i32 as i64))
                            as isize,
                    );
                vmnew = f + parmu * resmax;
                trured = vmold - vmnew;
                if parmu == zero
                    && f == *datmat.offset(
                        (mp - 1 as i32 as i64
                            + mpp * (np - 1 as i32 as i64))
                            as isize,
                    )
                {
                    prerem = prerec;
                    trured = *datmat.offset(
                        (mpp - 1 as i32 as i64
                            + mpp * (np - 1 as i32 as i64))
                            as isize,
                    ) - resmax;
                }
                ratio = if trured <= zero { one } else { zero };
                jdrop = 0 as i32 as i64;
                j = 1 as i32 as i64;
                while j <= n {
                    temp = zero;
                    i = 1 as i32 as i64;
                    while i <= n {
                        temp += *simi.offset(
                            (j - 1 as i32 as i64
                                + n * (i - 1 as i32 as i64))
                                as isize,
                        ) * *dx.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    temp = temp.abs();
                    if temp > ratio {
                        jdrop = j;
                        ratio = temp;
                    }
                    *sigbar.offset((j - 1 as i32 as i64) as isize) =
                        temp * *vsig.offset((j - 1 as i32 as i64) as isize);
                    j += 1;
                }
                edgmax = delta * rho;
                l = 0 as i32 as i64;
                j = 1 as i32 as i64;
                while j <= n {
                    if *sigbar.offset((j - 1 as i32 as i64) as isize) >= parsig
                        || *sigbar.offset((j - 1 as i32 as i64) as isize)
                            >= *vsig.offset((j - 1 as i32 as i64) as isize)
                    {
                        temp = *veta.offset((j - 1 as i32 as i64) as isize);
                        if trured > zero {
                            temp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                let mut tempb: f64 = *dx
                                    .offset((i - 1 as i32 as i64) as isize)
                                    - *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (j - 1 as i32 as i64))
                                            as isize,
                                    );
                                temp += tempb * tempb;
                                i += 1;
                            }
                            temp = temp.sqrt();
                        }
                        if temp > edgmax {
                            l = j;
                            edgmax = temp;
                        }
                    }
                    j += 1;
                }
                if l > 0 as i32 as i64 {
                    jdrop = l;
                }
                if jdrop == 0 as i32 as i64 {
                    current_block = 12206176237338959610;
                } else {
                    temp = zero;
                    i = 1 as i32 as i64;
                    while i <= n {
                        *sim.offset(
                            (i - 1 as i32 as i64
                                + n * (jdrop - 1 as i32 as i64))
                                as isize,
                        ) = *dx.offset((i - 1 as i32 as i64) as isize);
                        temp += *simi.offset(
                            (jdrop - 1 as i32 as i64
                                + n * (i - 1 as i32 as i64))
                                as isize,
                        ) * *dx.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    i = 1 as i32 as i64;
                    while i <= n {
                        *simi.offset(
                            (jdrop - 1 as i32 as i64
                                + n * (i - 1 as i32 as i64))
                                as isize,
                        ) = *simi.offset(
                            (jdrop - 1 as i32 as i64
                                + n * (i - 1 as i32 as i64))
                                as isize,
                        ) / temp;
                        i += 1;
                    }
                    j = 1 as i32 as i64;
                    while j <= n {
                        if j != jdrop {
                            temp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp += *simi.offset(
                                    (j - 1 as i32 as i64
                                        + n * (i - 1 as i32 as i64))
                                        as isize,
                                ) * *dx
                                    .offset((i - 1 as i32 as i64) as isize);
                                i += 1;
                            }
                            i = 1 as i32 as i64;
                            while i <= n {
                                *simi.offset(
                                    (j - 1 as i32 as i64
                                        + n * (i - 1 as i32 as i64))
                                        as isize,
                                ) -= temp
                                    * *simi.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                        }
                        j += 1;
                    }
                    k = 1 as i32 as i64;
                    while k <= mpp {
                        *datmat.offset(
                            (k - 1 as i32 as i64
                                + mpp * (jdrop - 1 as i32 as i64))
                                as isize,
                        ) = *con.offset((k - 1 as i32 as i64) as isize);
                        k += 1;
                    }
                    if trured > zero && trured >= 0.1f64 * prerem {
                        current_block = 1854612888966593411;
                    } else {
                        current_block = 12206176237338959610;
                    }
                }
            } else {
                k = 1 as i32 as i64;
                while k <= mpp {
                    *datmat.offset(
                        (k - 1 as i32 as i64
                            + mpp * (jdrop - 1 as i32 as i64))
                            as isize,
                    ) = *con.offset((k - 1 as i32 as i64) as isize);
                    k += 1;
                }
                if !(nfvals > np) {
                    if jdrop <= n {
                        if *datmat.offset(
                            (mp - 1 as i32 as i64
                                + mpp * (np - 1 as i32 as i64))
                                as isize,
                        ) <= f
                        {
                            *x.offset((jdrop - 1 as i32 as i64) as isize) = *sim
                                .offset(
                                    (jdrop - 1 as i32 as i64
                                        + n * (np - 1 as i32 as i64))
                                        as isize,
                                );
                        } else {
                            *sim.offset(
                                (jdrop - 1 as i32 as i64
                                    + n * (np - 1 as i32 as i64))
                                    as isize,
                            ) = *x.offset((jdrop - 1 as i32 as i64) as isize);
                            k = 1 as i32 as i64;
                            while k <= mpp {
                                *datmat.offset(
                                    (k - 1 as i32 as i64
                                        + mpp * (jdrop - 1 as i32 as i64))
                                        as isize,
                                ) = *datmat.offset(
                                    (k - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                );
                                *datmat.offset(
                                    (k - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                ) = *con.offset((k - 1 as i32 as i64) as isize);
                                k += 1;
                            }
                            k = 1 as i32 as i64;
                            while k <= jdrop {
                                temp = zero;
                                *sim.offset(
                                    (jdrop - 1 as i32 as i64
                                        + n * (k - 1 as i32 as i64))
                                        as isize,
                                ) = -rho;
                                i = k;
                                while i <= jdrop {
                                    temp -= *simi.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                *simi.offset(
                                    (jdrop - 1 as i32 as i64
                                        + n * (k - 1 as i32 as i64))
                                        as isize,
                                ) = temp;
                                k += 1;
                            }
                        }
                    }
                    if nfvals <= n {
                        jdrop = nfvals;
                        *x.offset((jdrop - 1 as i32 as i64) as isize) += rho;
                        continue;
                    }
                }
                ibrnch = 1 as i32 as i64;
                current_block = 1854612888966593411;
            }
            'c_3180: loop {
                match current_block {
                    12206176237338959610 => {
                        if iflag == 0 as i32 as i64 {
                            ibrnch = 0 as i32 as i64;
                            current_block = 1854612888966593411;
                        } else if rho > rhoend {
                            rho = 0.5f64 * rho;
                            if rho <= 1.5f64 * rhoend {
                                rho = rhoend;
                            }
                            if parmu > zero {
                                let mut cmin: f64 = 0.;
                                let mut cmax: f64 = 0.;
                                let mut denom: f64 = 0.;
                                denom = zero;
                                k = 1 as i32 as i64;
                                while k <= mp {
                                    cmin = *datmat.offset(
                                        (k - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                                    cmax = cmin;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        temp = *datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                        cmin = if cmin <= temp { cmin } else { temp };
                                        cmax = if cmax >= temp { cmax } else { temp };
                                        i += 1;
                                    }
                                    if k <= m && cmin < 0.5f64 * cmax {
                                        temp = (if cmax >= zero { cmax } else { zero }) - cmin;
                                        if denom <= zero {
                                            denom = temp;
                                        } else {
                                            denom = if denom <= temp { denom } else { temp };
                                        }
                                    }
                                    k += 1;
                                }
                                if denom == zero {
                                    parmu = zero;
                                } else if cmax - cmin < parmu * denom {
                                    parmu = (cmax - cmin) / denom;
                                }
                            }
                            if iprint >= 2 as i32 as i64 {
                                println!(
                                    "\n   Reduction in RHO to {}  and PARMU = {}\n",
                                    rho, parmu,
                                );
                                if iprint == 2 as i32 as i64 {
                                    print_calcfc(
                                        n,
                                        nfvals,
                                        *datmat.offset(
                                            (mp - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        ),
                                        *datmat.offset(
                                            (mpp - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        ),
                                        &mut *sim.offset(
                                            (n * (np - 1 as i32 as i64)) as isize,
                                        )
                                            as *mut f64
                                            as *const f64,
                                    );
                                }
                            }
                            current_block = 1854612888966593411;
                        } else {
                            if iprint >= 1 as i32 as i64 {
                                println!("\n   Normal return from subroutine COBYLA\n");
                            }
                            status = 0 as i32;
                            if ifull == 1 as i32 as i64 {
                                current_block = 1854573226567202969;
                                break 'c_3166;
                            } else {
                                current_block = 18034400475116118263;
                                break 'c_3166;
                            }
                        }
                    }
                    _ => {
                        phimin = *datmat.offset(
                            (mp - 1 as i32 as i64
                                + mpp * (np - 1 as i32 as i64))
                                as isize,
                        ) + parmu
                            * *datmat.offset(
                                (mpp - 1 as i32 as i64
                                    + mpp * (np - 1 as i32 as i64))
                                    as isize,
                            );
                        nbest = np;
                        j = 1 as i32 as i64;
                        while j <= n {
                            temp = *datmat.offset(
                                (mp - 1 as i32 as i64
                                    + mpp * (j - 1 as i32 as i64))
                                    as isize,
                            ) + parmu
                                * *datmat.offset(
                                    (mpp - 1 as i32 as i64
                                        + mpp * (j - 1 as i32 as i64))
                                        as isize,
                                );
                            if temp < phimin {
                                nbest = j;
                                phimin = temp;
                            } else if temp == phimin && parmu == zero {
                                if *datmat.offset(
                                    (mpp - 1 as i32 as i64
                                        + mpp * (j - 1 as i32 as i64))
                                        as isize,
                                ) < *datmat.offset(
                                    (mpp - 1 as i32 as i64
                                        + mpp * (nbest - 1 as i32 as i64))
                                        as isize,
                                ) {
                                    nbest = j;
                                }
                            }
                            j += 1;
                        }
                        if nbest <= n {
                            i = 1 as i32 as i64;
                            while i <= mpp {
                                temp = *datmat.offset(
                                    (i - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                );
                                *datmat.offset(
                                    (i - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                ) = *datmat.offset(
                                    (i - 1 as i32 as i64
                                        + mpp * (nbest - 1 as i32 as i64))
                                        as isize,
                                );
                                *datmat.offset(
                                    (i - 1 as i32 as i64
                                        + mpp * (nbest - 1 as i32 as i64))
                                        as isize,
                                ) = temp;
                                i += 1;
                            }
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp = *sim.offset(
                                    (i - 1 as i32 as i64
                                        + n * (nbest - 1 as i32 as i64))
                                        as isize,
                                );
                                *sim.offset(
                                    (i - 1 as i32 as i64
                                        + n * (nbest - 1 as i32 as i64))
                                        as isize,
                                ) = zero;
                                *sim.offset(
                                    (i - 1 as i32 as i64
                                        + n * (np - 1 as i32 as i64))
                                        as isize,
                                ) += temp;
                                tempa = zero;
                                k = 1 as i32 as i64;
                                while k <= n {
                                    *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) -= temp;
                                    tempa -= *simi.offset(
                                        (k - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    );
                                    k += 1;
                                }
                                *simi.offset(
                                    (nbest - 1 as i32 as i64
                                        + n * (i - 1 as i32 as i64))
                                        as isize,
                                ) = tempa;
                                i += 1;
                            }
                        }
                        error = zero;
                        i = 1 as i32 as i64;
                        while i <= n {
                            j = 1 as i32 as i64;
                            while j <= n {
                                temp = if i == j { -one } else { zero };
                                k = 1 as i32 as i64;
                                while k <= n {
                                    temp += *simi.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) * *sim.offset(
                                        (k - 1 as i32 as i64
                                            + n * (j - 1 as i32 as i64))
                                            as isize,
                                    );
                                    k += 1;
                                }
                                temp = temp.abs();
                                error = if error >= temp { error } else { temp };
                                j += 1;
                            }
                            i += 1;
                        }
                        if error > 0.1f64 {
                            status = -(1 as i32);
                            if iprint >= 1 as i32 as i64 {
                                println!(
                                    "Return from subroutine COBYLA because {:?}.\n",
                                    cobyla_reason(status),
                                );
                            }
                            current_block = 18034400475116118263;
                            break 'c_3166;
                        } else {
                            k = 1 as i32 as i64;
                            while k <= mp {
                                *con.offset((k - 1 as i32 as i64) as isize) =
                                    -*datmat.offset(
                                        (k - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                                j = 1 as i32 as i64;
                                while j <= n {
                                    *w.offset((j - 1 as i32 as i64) as isize) =
                                        *datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (j - 1 as i32 as i64))
                                                as isize,
                                        ) + *con.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                    j += 1;
                                }
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp = zero;
                                    j = 1 as i32 as i64;
                                    while j <= n {
                                        temp += *w.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        ) * *simi.offset(
                                            (j - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                        j += 1;
                                    }
                                    *a.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) = if k == mp { -temp } else { temp };
                                    i += 1;
                                }
                                k += 1;
                            }
                            iflag = 1 as i32 as i64;
                            parsig = alpha * rho;
                            pareta = beta * rho;
                            j = 1 as i32 as i64;
                            while j <= n {
                                let mut wsig: f64 = zero;
                                let mut weta: f64 = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    wsig += *simi.offset(
                                        (j - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) * *simi.offset(
                                        (j - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    );
                                    weta += *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (j - 1 as i32 as i64))
                                            as isize,
                                    ) * *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (j - 1 as i32 as i64))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                *vsig.offset((j - 1 as i32 as i64) as isize) =
                                    one / wsig.sqrt();
                                *veta.offset((j - 1 as i32 as i64) as isize) =
                                    weta.sqrt();
                                if *vsig.offset((j - 1 as i32 as i64) as isize)
                                    < parsig
                                    || *veta.offset((j - 1 as i32 as i64) as isize)
                                        > pareta
                                {
                                    iflag = 0 as i32 as i64;
                                }
                                j += 1;
                            }
                            if ibrnch == 1 as i32 as i64
                                || iflag == 1 as i32 as i64
                            {
                                let mut z: *mut f64 = w;
                                let mut zdota: *mut f64 = z.offset((n * n) as isize);
                                let mut vmc: *mut f64 = zdota.offset(n as isize);
                                let mut sdirn: *mut f64 = vmc.offset(mp as isize);
                                let mut dxnew: *mut f64 = sdirn.offset(n as isize);
                                let mut vmd: *mut f64 = dxnew.offset(n as isize);
                                trstlp(
                                    n,
                                    m,
                                    a as *const f64,
                                    con as *const f64,
                                    rho,
                                    dx,
                                    &mut ifull,
                                    iact,
                                    z,
                                    zdota,
                                    vmc,
                                    sdirn,
                                    dxnew,
                                    vmd,
                                );
                                if ifull == 0 as i32 as i64 {
                                    temp = zero;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        temp += *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        ) * *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        );
                                        i += 1;
                                    }
                                    if temp < 0.25f64 * rho * rho {
                                        ibrnch = 1 as i32 as i64;
                                        current_block = 12206176237338959610;
                                        continue;
                                    }
                                }
                                resnew = zero;
                                *con.offset((mp - 1 as i32 as i64) as isize) =
                                    zero;
                                sum = zero;
                                k = 1 as i32 as i64;
                                while k <= mp {
                                    sum = *con
                                        .offset((k - 1 as i32 as i64) as isize);
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        sum -= *a.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) * *dx.offset(
                                            (i - 1 as i32 as i64) as isize,
                                        );
                                        i += 1;
                                    }
                                    if k < mp {
                                        resnew = if resnew >= sum { resnew } else { sum };
                                    }
                                    k += 1;
                                }
                                prerec = *datmat.offset(
                                    (mpp - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                ) - resnew;
                                barmu = if prerec > zero { sum / prerec } else { zero };
                                if !(parmu < 1.5f64 * barmu) {
                                    break;
                                }
                                parmu = 2.0f64 * barmu;
                                if iprint >= 2 as i32 as i64 {
                                    println!("\n   Increase in PARMU to {}\n", parmu);
                                }
                                phi = *datmat.offset(
                                    (mp - 1 as i32 as i64
                                        + mpp * (np - 1 as i32 as i64))
                                        as isize,
                                ) + parmu
                                    * *datmat.offset(
                                        (mpp - 1 as i32 as i64
                                            + mpp * (np - 1 as i32 as i64))
                                            as isize,
                                    );
                                j = 1 as i32 as i64;
                                loop {
                                    if !(j <= n) {
                                        break 'c_3180;
                                    }
                                    temp = *datmat.offset(
                                        (mp - 1 as i32 as i64
                                            + mpp * (j - 1 as i32 as i64))
                                            as isize,
                                    ) + parmu
                                        * *datmat.offset(
                                            (mpp - 1 as i32 as i64
                                                + mpp * (j - 1 as i32 as i64))
                                                as isize,
                                        );
                                    if temp < phi {
                                        current_block = 1854612888966593411;
                                        break;
                                    }
                                    if temp == phi && parmu == zero {
                                        if *datmat.offset(
                                            (mpp - 1 as i32 as i64
                                                + mpp * (j - 1 as i32 as i64))
                                                as isize,
                                        ) < *datmat.offset(
                                            (mpp - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        ) {
                                            current_block = 1854612888966593411;
                                            break;
                                        }
                                    }
                                    j += 1;
                                }
                            } else {
                                jdrop = 0 as i32 as i64;
                                temp = pareta;
                                j = 1 as i32 as i64;
                                while j <= n {
                                    if *veta.offset((j - 1 as i32 as i64) as isize)
                                        > temp
                                    {
                                        jdrop = j;
                                        temp = *veta.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        );
                                    }
                                    j += 1;
                                }
                                if jdrop == 0 as i32 as i64 {
                                    j = 1 as i32 as i64;
                                    while j <= n {
                                        if *vsig
                                            .offset((j - 1 as i32 as i64) as isize)
                                            < temp
                                        {
                                            jdrop = j;
                                            temp = *vsig.offset(
                                                (j - 1 as i32 as i64) as isize,
                                            );
                                        }
                                        j += 1;
                                    }
                                }
                                temp = gamma
                                    * rho
                                    * *vsig.offset(
                                        (jdrop - 1 as i32 as i64) as isize,
                                    );
                                i = 1 as i32 as i64;
                                while i <= n {
                                    *dx.offset((i - 1 as i32 as i64) as isize) =
                                        temp * *simi.offset(
                                            (jdrop - 1 as i32 as i64
                                                + n * (i - 1 as i32 as i64))
                                                as isize,
                                        );
                                    i += 1;
                                }
                                cvmaxp = zero;
                                cvmaxm = zero;
                                sum = zero;
                                k = 1 as i32 as i64;
                                while k <= mp {
                                    sum = zero;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        sum = sum
                                            + *a.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (k - 1 as i32 as i64))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            );
                                        i += 1;
                                    }
                                    if k < mp {
                                        temp = *datmat.offset(
                                            (k - 1 as i32 as i64
                                                + mpp * (np - 1 as i32 as i64))
                                                as isize,
                                        );
                                        cvmaxp = if cvmaxp >= -sum - temp {
                                            cvmaxp
                                        } else {
                                            -sum - temp
                                        };
                                        cvmaxm = if cvmaxm >= sum - temp {
                                            cvmaxm
                                        } else {
                                            sum - temp
                                        };
                                    }
                                    k += 1;
                                }
                                if parmu * (cvmaxp - cvmaxm) > sum + sum {
                                    dxsign = -one;
                                } else {
                                    dxsign = one;
                                }
                                temp = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    *dx.offset((i - 1 as i32 as i64) as isize) *=
                                        dxsign;
                                    *sim.offset(
                                        (i - 1 as i32 as i64
                                            + n * (jdrop - 1 as i32 as i64))
                                            as isize,
                                    ) = *dx.offset((i - 1 as i32 as i64) as isize);
                                    temp += *simi.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) * *dx
                                        .offset((i - 1 as i32 as i64) as isize);
                                    i += 1;
                                }
                                i = 1 as i32 as i64;
                                while i <= n {
                                    *simi.offset(
                                        (jdrop - 1 as i32 as i64
                                            + n * (i - 1 as i32 as i64))
                                            as isize,
                                    ) /= temp;
                                    i += 1;
                                }
                                j = 1 as i32 as i64;
                                while j <= n {
                                    if j != jdrop {
                                        temp = zero;
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            temp += *simi.offset(
                                                (j - 1 as i32 as i64
                                                    + n * (i - 1 as i32 as i64))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            );
                                            i += 1;
                                        }
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            *simi.offset(
                                                (j - 1 as i32 as i64
                                                    + n * (i - 1 as i32 as i64))
                                                    as isize,
                                            ) -= temp
                                                * *simi.offset(
                                                    (jdrop - 1 as i32 as i64
                                                        + n * (i - 1 as i32
                                                            as i64))
                                                        as isize,
                                                );
                                            i += 1;
                                        }
                                    }
                                    *x.offset((j - 1 as i32 as i64) as isize) =
                                        *sim.offset(
                                            (j - 1 as i32 as i64
                                                + n * (np - 1 as i32 as i64))
                                                as isize,
                                        ) + *dx.offset(
                                            (j - 1 as i32 as i64) as isize,
                                        );
                                    j += 1;
                                }
                                continue 'c_3166;
                            }
                        }
                    }
                }
            }
            prerem = parmu * prerec - sum;
            i = 1 as i32 as i64;
            while i <= n {
                *x.offset((i - 1 as i32 as i64) as isize) = *sim.offset(
                    (i - 1 as i32 as i64
                        + n * (np - 1 as i32 as i64)) as isize,
                ) + *dx
                    .offset((i - 1 as i32 as i64) as isize);
                i += 1;
            }
            ibrnch = 1 as i32 as i64;
        }
    }
    match current_block {
        18034400475116118263 => {
            i = 1 as i32 as i64;
            while i <= n {
                *x.offset((i - 1 as i32 as i64) as isize) = *sim.offset(
                    (i - 1 as i32 as i64
                        + n * (np - 1 as i32 as i64)) as isize,
                );
                i += 1;
            }
            f = *datmat.offset(
                (mp - 1 as i32 as i64
                    + mpp * (np - 1 as i32 as i64)) as isize,
            );
            resmax = *datmat.offset(
                (mpp - 1 as i32 as i64
                    + mpp * (np - 1 as i32 as i64)) as isize,
            );
        }
        _ => {}
    }
    if iprint >= 1 as i32 as i64 {
        print_calcfc(
            n,
            nfvals,
            f,
            resmax,
            &mut *x.offset(0) as *mut f64 as *const f64,
        );
    }
    *_maxfun = nfvals;
    return status;
}

#[allow(clippy::too_many_arguments)]
unsafe fn trstlp(
    n: i64,
    m: i64,
    mut a: *const f64,
    mut b: *const f64,
    rho: f64,
    mut dx: *mut f64,
    mut ifull: *mut i64,
    mut iact: *mut i64,
    mut z: *mut f64,
    mut zdota: *mut f64,
    mut vmultc: *mut f64,
    mut sdirn: *mut f64,
    mut dxnew: *mut f64,
    mut vmultd: *mut f64,
) {
    let mut current_block: u64;
    let zero: f64 = 0.0f64;
    let one: f64 = 1.0f64;
    let tiny: f64 = 1.0e-6f64;
    let Op1: f64 = 0.1f64;
    let Op2: f64 = 0.2f64;
    let mut acca: f64 = 0.;
    let mut accb: f64 = 0.;
    let mut alpha: f64 = 0.;
    let mut beta: f64 = 0.;
    let mut dd: f64 = 0.;
    let mut optnew: f64 = 0.;
    let mut optold: f64 = 0.;
    let mut ratio: f64 = 0.;
    let mut resmax: f64 = 0.;
    let mut resold: f64 = 0.;
    let mut sd: f64 = 0.;
    let mut sp: f64 = 0.;
    let mut spabs: f64 = 0.;
    let mut ss: f64 = 0.;
    let mut step: f64 = 0.;
    let mut stpful: f64 = 0.;
    let mut sum: f64 = 0.;
    let mut sumabs: f64 = 0.;
    let mut temp: f64 = 0.;
    let mut tempa: f64 = 0.;
    let mut tempb: f64 = 0.;
    let mut tot: f64 = 0.;
    let mut vsave: f64 = 0.;
    let mut zdotv: f64 = 0.;
    let mut zdotw: f64 = 0.;
    let mut zdvabs: f64 = 0.;
    let mut zdwabs: f64 = 0.;
    let mut i: i64 = 0;
    let mut icon: i64 = 0;
    let mut icount: i64 = 0;
    let mut isave: i64 = 0;
    let mut j: i64 = 0;
    let mut k: i64 = 0;
    let mut kk: i64 = 0;
    let mut kl: i64 = 0;
    let mut kp: i64 = 0;
    let mut kw: i64 = 0;
    let mut mcon: i64 = 0;
    let mut nact: i64 = 0;
    let mut nactx: i64 = 0;
    *ifull = 1 as i32 as i64;
    icon = 0 as i32 as i64;
    mcon = m;
    nact = 0 as i32 as i64;
    resmax = zero;
    resold = zero;
    i = 1 as i32 as i64;
    while i <= n {
        j = 1 as i32 as i64;
        while j <= n {
            *z.offset(
                (i - 1 as i32 as i64 + n * (j - 1 as i32 as i64))
                    as isize,
            ) = zero;
            j += 1;
        }
        *z.offset(
            (i - 1 as i32 as i64 + n * (i - 1 as i32 as i64))
                as isize,
        ) = one;
        *dx.offset((i - 1 as i32 as i64) as isize) = zero;
        i += 1;
    }
    if m >= 1 as i32 as i64 {
        k = 1 as i32 as i64;
        while k <= m {
            if *b.offset((k - 1 as i32 as i64) as isize) > resmax {
                resmax = *b.offset((k - 1 as i32 as i64) as isize);
                icon = k;
            }
            k += 1;
        }
        k = 1 as i32 as i64;
        while k <= m {
            *iact.offset((k - 1 as i32 as i64) as isize) = k;
            *vmultc.offset((k - 1 as i32 as i64) as isize) =
                resmax - *b.offset((k - 1 as i32 as i64) as isize);
            k += 1;
        }
    }
    if resmax == zero {
        current_block = 16746262731645592041;
    } else {
        i = 1 as i32 as i64;
        while i <= n {
            *sdirn.offset((i - 1 as i32 as i64) as isize) = zero;
            i += 1;
        }
        current_block = 10785318058620693244;
    }
    'c_5472: loop {
        match current_block {
            16746262731645592041 => {
                mcon = m + 1 as i32 as i64;
                icon = mcon;
                *iact.offset((mcon - 1 as i32 as i64) as isize) = mcon;
                *vmultc.offset((mcon - 1 as i32 as i64) as isize) = zero;
                current_block = 10785318058620693244;
            }
            _ => {
                optold = zero;
                icount = 0 as i32 as i64;
                loop {
                    if mcon == m {
                        optnew = resmax;
                    } else {
                        optnew = zero;
                        i = 1 as i32 as i64;
                        while i <= n {
                            optnew -= *dx.offset((i - 1 as i32 as i64) as isize)
                                * *a.offset(
                                    (i - 1 as i32 as i64
                                        + n * (mcon - 1 as i32 as i64))
                                        as isize,
                                );
                            i += 1;
                        }
                    }
                    if icount == 0 as i32 as i64 || optnew < optold {
                        optold = optnew;
                        nactx = nact;
                        icount = 3 as i32 as i64;
                    } else if nact > nactx {
                        nactx = nact;
                        icount = 3 as i32 as i64;
                    } else {
                        icount -= 1;
                        if icount == 0 as i32 as i64 {
                            break;
                        }
                    }
                    if icon <= nact {
                        if icon < nact {
                            isave =
                                *iact.offset((icon - 1 as i32 as i64) as isize);
                            vsave =
                                *vmultc.offset((icon - 1 as i32 as i64) as isize);
                            k = icon;
                            loop {
                                kp = k + 1 as i32 as i64;
                                kk = *iact.offset((kp - 1 as i32 as i64) as isize);
                                sp = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    sp += *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) * *a.offset(
                                        (i - 1 as i32 as i64
                                            + n * (kk - 1 as i32 as i64))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                temp = (sp * sp
                                    + *zdota
                                        .offset((kp - 1 as i32 as i64) as isize)
                                        * *zdota.offset(
                                            (kp - 1 as i32 as i64) as isize,
                                        ))
                                .sqrt();
                                alpha = *zdota
                                    .offset((kp - 1 as i32 as i64) as isize)
                                    / temp;
                                beta = sp / temp;
                                *zdota.offset((kp - 1 as i32 as i64) as isize) =
                                    alpha
                                        * *zdota.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                *zdota.offset((k - 1 as i32 as i64) as isize) =
                                    temp;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp = alpha
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (kp - 1 as i32 as i64))
                                                as isize,
                                        )
                                        + beta
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (k - 1 as i32 as i64))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (kp - 1 as i32 as i64))
                                            as isize,
                                    ) = alpha
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        )
                                        - beta
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (kp - 1 as i32 as i64))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                                *iact.offset((k - 1 as i32 as i64) as isize) = kk;
                                *vmultc.offset((k - 1 as i32 as i64) as isize) =
                                    *vmultc
                                        .offset((kp - 1 as i32 as i64) as isize);
                                k = kp;
                                if !(k < nact) {
                                    break;
                                }
                            }
                            *iact.offset((k - 1 as i32 as i64) as isize) = isave;
                            *vmultc.offset((k - 1 as i32 as i64) as isize) = vsave;
                        }
                        nact -= 1;
                        if mcon > m {
                            current_block = 14831642685178214089;
                        } else {
                            temp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp += *sdirn
                                    .offset((i - 1 as i32 as i64) as isize)
                                    * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nact + 1 as i32 as i64
                                                - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                            i = 1 as i32 as i64;
                            while i <= n {
                                *sdirn.offset((i - 1 as i32 as i64) as isize) -=
                                    temp * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nact + 1 as i32 as i64
                                                - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                            current_block = 18040478258510813106;
                        }
                    } else {
                        kk = *iact.offset((icon - 1 as i32 as i64) as isize);
                        i = 1 as i32 as i64;
                        while i <= n {
                            *dxnew.offset((i - 1 as i32 as i64) as isize) = *a
                                .offset(
                                    (i - 1 as i32 as i64
                                        + n * (kk - 1 as i32 as i64))
                                        as isize,
                                );
                            i += 1;
                        }
                        tot = zero;
                        k = n;
                        while k > nact {
                            sp = zero;
                            spabs = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp = *z.offset(
                                    (i - 1 as i32 as i64
                                        + n * (k - 1 as i32 as i64))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as i32 as i64) as isize);
                                sp += temp;
                                spabs += temp.abs();
                                i += 1;
                            }
                            acca = spabs + Op1 * sp.abs();
                            accb = spabs + Op2 * sp.abs();
                            if spabs >= acca || acca >= accb {
                                sp = zero;
                            }
                            if tot == zero {
                                tot = sp;
                            } else {
                                kp = k + 1 as i32 as i64;
                                temp = (sp * sp + tot * tot).sqrt();
                                alpha = sp / temp;
                                beta = tot / temp;
                                tot = temp;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp = alpha
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        )
                                        + beta
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (kp - 1 as i32 as i64))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (kp - 1 as i32 as i64))
                                            as isize,
                                    ) = alpha
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (kp - 1 as i32 as i64))
                                                as isize,
                                        )
                                        - beta
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (k - 1 as i32 as i64))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                            }
                            k -= 1;
                        }
                        if tot != zero {
                            nact += 1;
                            *zdota.offset((nact - 1 as i32 as i64) as isize) = tot;
                            *vmultc.offset((icon - 1 as i32 as i64) as isize) =
                                *vmultc.offset((nact - 1 as i32 as i64) as isize);
                            *vmultc.offset((nact - 1 as i32 as i64) as isize) =
                                zero;
                        } else {
                            ratio = -one;
                            k = nact;
                            loop {
                                zdotv = zero;
                                zdvabs = zero;
                                i = 1 as i32 as i64;
                                while i <= n {
                                    temp = *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    ) * *dxnew
                                        .offset((i - 1 as i32 as i64) as isize);
                                    zdotv += temp;
                                    zdvabs += temp.abs();
                                    i += 1;
                                }
                                acca = zdvabs + Op1 * zdotv.abs();
                                accb = zdvabs + Op2 * zdotv.abs();
                                if zdvabs < acca && acca < accb {
                                    temp = zdotv
                                        / *zdota.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                    if temp > zero
                                        && *iact
                                            .offset((k - 1 as i32 as i64) as isize)
                                            <= m
                                    {
                                        tempa = *vmultc.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        ) / temp;
                                        if ratio < zero || tempa < ratio {
                                            ratio = tempa;
                                        }
                                    }
                                    if k >= 2 as i32 as i64 {
                                        kw = *iact.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                        i = 1 as i32 as i64;
                                        while i <= n {
                                            *dxnew.offset(
                                                (i - 1 as i32 as i64) as isize,
                                            ) -= temp
                                                * *a.offset(
                                                    (i - 1 as i32 as i64
                                                        + n * (kw
                                                            - 1 as i32 as i64))
                                                        as isize,
                                                );
                                            i += 1;
                                        }
                                    }
                                    *vmultd
                                        .offset((k - 1 as i32 as i64) as isize) =
                                        temp;
                                } else {
                                    *vmultd
                                        .offset((k - 1 as i32 as i64) as isize) =
                                        zero;
                                }
                                k -= 1;
                                if !(k > 0 as i32 as i64) {
                                    break;
                                }
                            }
                            if ratio < zero {
                                break;
                            }
                            k = 1 as i32 as i64;
                            while k <= nact {
                                tempb = *vmultc
                                    .offset((k - 1 as i32 as i64) as isize)
                                    - ratio
                                        * *vmultd.offset(
                                            (k - 1 as i32 as i64) as isize,
                                        );
                                *vmultc.offset((k - 1 as i32 as i64) as isize) =
                                    if zero >= tempb { zero } else { tempb };
                                k += 1;
                            }
                            if icon < nact {
                                isave = *iact
                                    .offset((icon - 1 as i32 as i64) as isize);
                                vsave = *vmultc
                                    .offset((icon - 1 as i32 as i64) as isize);
                                k = icon;
                                loop {
                                    kp = k + 1 as i32 as i64;
                                    kw = *iact
                                        .offset((kp - 1 as i32 as i64) as isize);
                                    sp = zero;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        sp += *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) * *a.offset(
                                            (i - 1 as i32 as i64
                                                + n * (kw - 1 as i32 as i64))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    temp = (sp * sp
                                        + *zdota.offset(
                                            (kp - 1 as i32 as i64) as isize,
                                        ) * *zdota.offset(
                                            (kp - 1 as i32 as i64) as isize,
                                        ))
                                    .sqrt();
                                    alpha = *zdota
                                        .offset((kp - 1 as i32 as i64) as isize)
                                        / temp;
                                    beta = sp / temp;
                                    *zdota
                                        .offset((kp - 1 as i32 as i64) as isize) =
                                        alpha
                                            * *zdota.offset(
                                                (k - 1 as i32 as i64) as isize,
                                            );
                                    *zdota
                                        .offset((k - 1 as i32 as i64) as isize) =
                                        temp;
                                    i = 1 as i32 as i64;
                                    while i <= n {
                                        temp = alpha
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (kp - 1 as i32 as i64))
                                                    as isize,
                                            )
                                            + beta
                                                * *z.offset(
                                                    (i - 1 as i32 as i64
                                                        + n * (k - 1 as i32
                                                            as i64))
                                                        as isize,
                                                );
                                        *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (kp - 1 as i32 as i64))
                                                as isize,
                                        ) = alpha
                                            * *z.offset(
                                                (i - 1 as i32 as i64
                                                    + n * (k - 1 as i32 as i64))
                                                    as isize,
                                            )
                                            - beta
                                                * *z.offset(
                                                    (i - 1 as i32 as i64
                                                        + n * (kp
                                                            - 1 as i32 as i64))
                                                        as isize,
                                                );
                                        *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        ) = temp;
                                        i += 1;
                                    }
                                    *iact.offset((k - 1 as i32 as i64) as isize) =
                                        kw;
                                    *vmultc
                                        .offset((k - 1 as i32 as i64) as isize) =
                                        *vmultc.offset(
                                            (kp - 1 as i32 as i64) as isize,
                                        );
                                    k = kp;
                                    if !(k < nact) {
                                        break;
                                    }
                                }
                                *iact.offset((k - 1 as i32 as i64) as isize) =
                                    isave;
                                *vmultc.offset((k - 1 as i32 as i64) as isize) =
                                    vsave;
                            }
                            temp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp += *z.offset(
                                    (i - 1 as i32 as i64
                                        + n * (nact - 1 as i32 as i64))
                                        as isize,
                                ) * *a.offset(
                                    (i - 1 as i32 as i64
                                        + n * (kk - 1 as i32 as i64))
                                        as isize,
                                );
                                i += 1;
                            }
                            if temp == zero {
                                break;
                            }
                            *zdota.offset((nact - 1 as i32 as i64) as isize) =
                                temp;
                            *vmultc.offset((icon - 1 as i32 as i64) as isize) =
                                zero;
                            *vmultc.offset((nact - 1 as i32 as i64) as isize) =
                                ratio;
                        }
                        *iact.offset((icon - 1 as i32 as i64) as isize) =
                            *iact.offset((nact - 1 as i32 as i64) as isize);
                        *iact.offset((nact - 1 as i32 as i64) as isize) = kk;
                        if mcon > m && kk != mcon {
                            k = nact - 1 as i32 as i64;
                            sp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                sp += *z.offset(
                                    (i - 1 as i32 as i64
                                        + n * (k - 1 as i32 as i64))
                                        as isize,
                                ) * *a.offset(
                                    (i - 1 as i32 as i64
                                        + n * (kk - 1 as i32 as i64))
                                        as isize,
                                );
                                i += 1;
                            }
                            temp = (sp * sp
                                + *zdota
                                    .offset((nact - 1 as i32 as i64) as isize)
                                    * *zdota.offset(
                                        (nact - 1 as i32 as i64) as isize,
                                    ))
                            .sqrt();
                            alpha = *zdota
                                .offset((nact - 1 as i32 as i64) as isize)
                                / temp;
                            beta = sp / temp;
                            *zdota.offset((nact - 1 as i32 as i64) as isize) =
                                alpha
                                    * *zdota
                                        .offset((k - 1 as i32 as i64) as isize);
                            *zdota.offset((k - 1 as i32 as i64) as isize) = temp;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp = alpha
                                    * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nact - 1 as i32 as i64))
                                            as isize,
                                    )
                                    + beta
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (k - 1 as i32 as i64))
                                                as isize,
                                        );
                                *z.offset(
                                    (i - 1 as i32 as i64
                                        + n * (nact - 1 as i32 as i64))
                                        as isize,
                                ) = alpha
                                    * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (k - 1 as i32 as i64))
                                            as isize,
                                    )
                                    - beta
                                        * *z.offset(
                                            (i - 1 as i32 as i64
                                                + n * (nact - 1 as i32 as i64))
                                                as isize,
                                        );
                                *z.offset(
                                    (i - 1 as i32 as i64
                                        + n * (k - 1 as i32 as i64))
                                        as isize,
                                ) = temp;
                                i += 1;
                            }
                            *iact.offset((nact - 1 as i32 as i64) as isize) =
                                *iact.offset((k - 1 as i32 as i64) as isize);
                            *iact.offset((k - 1 as i32 as i64) as isize) = kk;
                            temp = *vmultc.offset((k - 1 as i32 as i64) as isize);
                            *vmultc.offset((k - 1 as i32 as i64) as isize) =
                                *vmultc.offset((nact - 1 as i32 as i64) as isize);
                            *vmultc.offset((nact - 1 as i32 as i64) as isize) =
                                temp;
                        }
                        if mcon > m {
                            current_block = 14831642685178214089;
                        } else {
                            kk = *iact.offset((nact - 1 as i32 as i64) as isize);
                            temp = zero;
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp += *sdirn
                                    .offset((i - 1 as i32 as i64) as isize)
                                    * *a.offset(
                                        (i - 1 as i32 as i64
                                            + n * (kk - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                            temp = (temp - one)
                                / *zdota.offset((nact - 1 as i32 as i64) as isize);
                            i = 1 as i32 as i64;
                            while i <= n {
                                *sdirn.offset((i - 1 as i32 as i64) as isize) -=
                                    temp * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nact - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                            current_block = 18040478258510813106;
                        }
                    }
                    match current_block {
                        14831642685178214089 => {
                            temp = one
                                / *zdota.offset((nact - 1 as i32 as i64) as isize);
                            i = 1 as i32 as i64;
                            while i <= n {
                                *sdirn.offset((i - 1 as i32 as i64) as isize) =
                                    temp * *z.offset(
                                        (i - 1 as i32 as i64
                                            + n * (nact - 1 as i32 as i64))
                                            as isize,
                                    );
                                i += 1;
                            }
                        }
                        _ => {}
                    }
                    dd = rho * rho;
                    sd = zero;
                    ss = zero;
                    i = 1 as i32 as i64;
                    while i <= n {
                        if (*dx.offset((i - 1 as i32 as i64) as isize)).abs()
                            >= tiny * rho
                        {
                            dd -= *dx.offset((i - 1 as i32 as i64) as isize)
                                * *dx.offset((i - 1 as i32 as i64) as isize);
                        }
                        sd += *sdirn.offset((i - 1 as i32 as i64) as isize)
                            * *dx.offset((i - 1 as i32 as i64) as isize);
                        ss += *sdirn.offset((i - 1 as i32 as i64) as isize)
                            * *sdirn.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    if dd <= zero {
                        break;
                    }
                    temp = (ss * dd).sqrt();
                    if sd.abs() >= tiny * temp {
                        temp = (ss * dd + sd * sd).sqrt();
                    }
                    stpful = dd / (temp + sd);
                    step = stpful;
                    if mcon == m {
                        acca = step + Op1 * resmax;
                        accb = step + Op2 * resmax;
                        if step >= acca || acca >= accb {
                            current_block = 16746262731645592041;
                            continue 'c_5472;
                        }
                        step = if step <= resmax { step } else { resmax };
                    }
                    i = 1 as i32 as i64;
                    while i <= n {
                        *dxnew.offset((i - 1 as i32 as i64) as isize) = *dx
                            .offset((i - 1 as i32 as i64) as isize)
                            + step * *sdirn.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    if mcon == m {
                        resold = resmax;
                        resmax = zero;
                        k = 1 as i32 as i64;
                        while k <= nact {
                            kk = *iact.offset((k - 1 as i32 as i64) as isize);
                            temp = *b.offset((kk - 1 as i32 as i64) as isize);
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp -= *a.offset(
                                    (i - 1 as i32 as i64
                                        + n * (kk - 1 as i32 as i64))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as i32 as i64) as isize);
                                i += 1;
                            }
                            resmax = if resmax >= temp { resmax } else { temp };
                            k += 1;
                        }
                    }
                    k = nact;
                    loop {
                        zdotw = zero;
                        zdwabs = zero;
                        i = 1 as i32 as i64;
                        while i <= n {
                            temp = *z.offset(
                                (i - 1 as i32 as i64
                                    + n * (k - 1 as i32 as i64))
                                    as isize,
                            ) * *dxnew
                                .offset((i - 1 as i32 as i64) as isize);
                            zdotw += temp;
                            zdwabs += temp.abs();
                            i += 1;
                        }
                        acca = zdwabs + Op1 * zdotw.abs();
                        accb = zdwabs + Op2 * zdotw.abs();
                        if zdwabs >= acca || acca >= accb {
                            zdotw = zero;
                        }
                        *vmultd.offset((k - 1 as i32 as i64) as isize) =
                            zdotw / *zdota.offset((k - 1 as i32 as i64) as isize);
                        if k < 2 as i32 as i64 {
                            break;
                        }
                        kk = *iact.offset((k - 1 as i32 as i64) as isize);
                        i = 1 as i32 as i64;
                        while i <= n {
                            *dxnew.offset((i - 1 as i32 as i64) as isize) -=
                                *vmultd.offset((k - 1 as i32 as i64) as isize)
                                    * *a.offset(
                                        (i - 1 as i32 as i64
                                            + n * (kk - 1 as i32 as i64))
                                            as isize,
                                    );
                            i += 1;
                        }
                        k -= 1;
                    }
                    if mcon > m
                        && *vmultd.offset((nact - 1 as i32 as i64) as isize) < zero
                    {
                        *vmultd.offset((nact - 1 as i32 as i64) as isize) = zero;
                    }
                    i = 1 as i32 as i64;
                    while i <= n {
                        *dxnew.offset((i - 1 as i32 as i64) as isize) = *dx
                            .offset((i - 1 as i32 as i64) as isize)
                            + step * *sdirn.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    if mcon > nact {
                        kl = nact + 1 as i32 as i64;
                        k = kl;
                        while k <= mcon {
                            kk = *iact.offset((k - 1 as i32 as i64) as isize);
                            sum = resmax
                                - *b.offset((kk - 1 as i32 as i64) as isize);
                            sumabs = resmax
                                + (*b.offset((kk - 1 as i32 as i64) as isize))
                                    .abs();
                            i = 1 as i32 as i64;
                            while i <= n {
                                temp = *a.offset(
                                    (i - 1 as i32 as i64
                                        + n * (kk - 1 as i32 as i64))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as i32 as i64) as isize);
                                sum += temp;
                                sumabs += temp.abs();
                                i += 1;
                            }
                            acca = sumabs + Op1 * sum.abs();
                            accb = sumabs + Op2 * sum.abs();
                            if sumabs >= acca || acca >= accb {
                                sum = zero;
                            }
                            *vmultd.offset((k - 1 as i32 as i64) as isize) = sum;
                            k += 1;
                        }
                    }
                    ratio = one;
                    icon = 0 as i32 as i64;
                    k = 1 as i32 as i64;
                    while k <= mcon {
                        if *vmultd.offset((k - 1 as i32 as i64) as isize) < zero {
                            temp = *vmultc.offset((k - 1 as i32 as i64) as isize)
                                / (*vmultc.offset((k - 1 as i32 as i64) as isize)
                                    - *vmultd
                                        .offset((k - 1 as i32 as i64) as isize));
                            if temp < ratio {
                                ratio = temp;
                                icon = k;
                            }
                        }
                        k += 1;
                    }
                    temp = one - ratio;
                    i = 1 as i32 as i64;
                    while i <= n {
                        *dx.offset((i - 1 as i32 as i64) as isize) = temp
                            * *dx.offset((i - 1 as i32 as i64) as isize)
                            + ratio
                                * *dxnew.offset((i - 1 as i32 as i64) as isize);
                        i += 1;
                    }
                    k = 1 as i32 as i64;
                    while k <= mcon {
                        tempb = temp
                            * *vmultc.offset((k - 1 as i32 as i64) as isize)
                            + ratio
                                * *vmultd.offset((k - 1 as i32 as i64) as isize);
                        *vmultc.offset((k - 1 as i32 as i64) as isize) =
                            if zero >= tempb { zero } else { tempb };
                        k += 1;
                    }
                    if mcon == m {
                        resmax = resold + ratio * (resmax - resold);
                    }
                    if icon > 0 as i32 as i64 {
                        continue;
                    }
                    if step == stpful {
                        break 'c_5472;
                    } else {
                        current_block = 16746262731645592041;
                        continue 'c_5472;
                    }
                }
                if mcon == m {
                    current_block = 16746262731645592041;
                    continue;
                }
                *ifull = 0 as i32 as i64;
                break;
            }
        }
    }
}

pub(crate) fn cobyla_reason(mut status: i32) -> &'static str {
    match status {
        1 => "user requested to compute F(X) and C(X)",
        0 => "algorithm was successful",
        -1 => "rounding errors are becoming damaging",
        -2 => "MAXFUN limit has been reached",
        -3 => "illegal NULL address",
        -4 => "unexpected parameter or corrupted workspace",
        _ => "unknown status",
    }
}
unsafe fn print_calcfc(
    mut n: i64,
    mut nfvals: i64,
    mut f: f64,
    mut maxcv: f64,
    mut x: *const f64,
) {
    let mut i: i64 = 0;
    println!(
        "\n   NFVALS ={}   F ={}    MAXCV ={}\n   X ={}",
        nfvals as i32,
        f,
        maxcv,
        *x.offset(0 as i32 as isize),
    );
    i = 1 as i32 as i64;
    while i < n {
        // let fmt = if i % 5 as i32 as i64 == 0 as i32 as i64 {
        //     "\n{}"
        // } else {
        //     "{}"
        // };
        println!("{}", *x.offset(i as isize));
        i += 1;
    }
    println!("\n");
}
