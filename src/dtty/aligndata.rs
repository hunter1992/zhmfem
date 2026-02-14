use std::alloc::{alloc, dealloc, handle_alloc_error, Layout};
use std::arch::x86_64::*;
use std::fmt;
use std::ops::{Deref, DerefMut};
use std::ptr::{self, NonNull};

// ============================================================================
// 1. FastCast 协议：极致性能版
// ============================================================================

pub trait FastCast<T> {
    unsafe fn fast_cast(src: *const Self, dst: *mut T, len: usize);
}

// --- f32 -> f64 (Expansion) ---
impl FastCast<f64> for f32 {
    #[inline(always)]
    unsafe fn fast_cast(src: *const f32, dst: *mut f64, len: usize) {
        let mut i = 0;

        #[cfg(target_feature = "avx2")]
        {
            unsafe {
                while i + 8 <= len {
                    let s_vec = _mm256_load_ps(src.add(i));
                    let d_low = _mm256_cvtps_pd(_mm256_castps256_ps128(s_vec));
                    let d_high = _mm256_cvtps_pd(_mm256_extractf128_ps(s_vec, 1));

                    _mm256_store_pd(dst.add(i), d_low);
                    _mm256_store_pd(dst.add(i + 4), d_high);
                    i += 8;
                }
            }
        }

        if i < len {
            for j in i..len {
                unsafe {
                    *dst.add(j) = *src.add(j) as f64;
                }
            }
        }
    }
}

// --- f64 -> f32 (Compression) ---
impl FastCast<f32> for f64 {
    #[inline(always)]
    unsafe fn fast_cast(src: *const f64, dst: *mut f32, len: usize) {
        let mut i = 0;

        #[cfg(target_feature = "avx2")]
        unsafe {
            while i + 16 <= len {
                let d0 = _mm256_load_pd(src.add(i));
                let d1 = _mm256_load_pd(src.add(i + 4));
                let d2 = _mm256_load_pd(src.add(i + 8));
                let d3 = _mm256_load_pd(src.add(i + 12));

                let f0 = _mm256_cvtpd_ps(d0);
                let f1 = _mm256_cvtpd_ps(d1);
                let res0 = _mm256_insertf128_ps(_mm256_castps128_ps256(f0), f1, 1);
                _mm256_store_ps(dst.add(i), res0);

                let f2 = _mm256_cvtpd_ps(d2);
                let f3 = _mm256_cvtpd_ps(d3);
                let res1 = _mm256_insertf128_ps(_mm256_castps128_ps256(f2), f3, 1);
                _mm256_store_ps(dst.add(i + 8), res1);

                i += 16;
            }
        }

        for j in i..len {
            unsafe {
                *dst.add(j) = *src.add(j) as f32;
            }
        }
    }
}

// --- i32 -> usize (u64) ---
impl FastCast<usize> for i32 {
    #[inline(always)]
    unsafe fn fast_cast(src: *const i32, dst: *mut usize, len: usize) {
        let mut i = 0;
        let d_ptr = dst as *mut u64;

        #[cfg(target_feature = "avx2")]
        unsafe {
            while i + 8 <= len {
                let v_i32 = _mm256_load_si256(src.add(i) as *const __m256i);
                let v_u64_lo = _mm256_cvtepu32_epi64(_mm256_castsi256_si128(v_i32));
                let v_u64_hi = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(v_i32, 1));

                _mm256_store_si256(d_ptr.add(i) as *mut __m256i, v_u64_lo);
                _mm256_store_si256(d_ptr.add(i + 4) as *mut __m256i, v_u64_hi);
                i += 8;
            }
        }

        for j in i..len {
            unsafe {
                *d_ptr.add(j) = *src.add(j) as u64;
            }
        }
    }
}

macro_rules! impl_copy {
    ($($t:ty),*) => {
        $(
            impl FastCast<$t> for $t {
                #[inline(always)]
                unsafe fn fast_cast(src: *const $t, dst: *mut $t, len: usize) {
                    unsafe { ptr::copy_nonoverlapping(src, dst, len); }
                }
            }
        )*
    };
}
impl_copy!(f32, f64, i32, usize, u64, i64);

// ============================================================================
// 2. AlignedVec：容器部分
// ============================================================================

pub struct AlignedVec<T> {
    ptr: NonNull<T>,
    len: usize,
}

impl<T> AlignedVec<T> {
    const ALIGN: usize = 64;

    pub fn new(len: usize) -> Self {
        if len == 0 {
            return Self {
                ptr: NonNull::dangling(),
                len: 0,
            };
        }
        unsafe {
            let layout =
                Layout::from_size_align_unchecked(len * std::mem::size_of::<T>(), Self::ALIGN);
            let ptr = alloc(layout) as *mut T;
            if ptr.is_null() {
                handle_alloc_error(layout);
            }
            Self {
                ptr: NonNull::new_unchecked(ptr),
                len,
            }
        }
    }

    // 在将数据交给 PARDISO 之前，先让 CPU “感知”这块内存
    #[inline(always)]
    pub fn warm_up(&self) {
        let p = self.as_ptr() as *const i8;
        let len = self.len() * std::mem::size_of::<T>();
        let mut i = 0;
        unsafe {
            while i < len {
                // 读取每个缓存行的一个字节，强制触发硬件预取和页面加载
                std::ptr::read_volatile(p.add(i));
                i += 64;
            }
        }
    }

    #[inline(always)]
    pub fn from_slice<S>(src: &[S]) -> Self
    where
        S: FastCast<T>,
    {
        let dst = Self::new(src.len());
        if !src.is_empty() {
            unsafe {
                S::fast_cast(src.as_ptr(), dst.ptr.as_ptr(), src.len());
            }
        }
        dst
    }

    #[inline(always)]
    pub fn into_vec<Target>(self) -> Vec<Target>
    where
        T: FastCast<Target>,
    {
        let mut out = Vec::with_capacity(self.len);
        unsafe {
            T::fast_cast(self.ptr.as_ptr(), out.as_mut_ptr(), self.len);
            out.set_len(self.len);
        }
        out
    }

    #[inline(always)]
    pub fn as_ptr(&self) -> *const T {
        self.ptr.as_ptr()
    }
    #[inline(always)]
    pub fn as_mut_ptr(&mut self) -> *mut T {
        self.ptr.as_ptr()
    }
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }
}

impl<T> Drop for AlignedVec<T> {
    #[inline]
    fn drop(&mut self) {
        if self.len > 0 {
            unsafe {
                let layout = Layout::from_size_align_unchecked(
                    self.len * std::mem::size_of::<T>(),
                    Self::ALIGN,
                );
                dealloc(self.ptr.as_ptr() as *mut u8, layout);
            }
        }
    }
}

impl<T: fmt::Debug> fmt::Debug for AlignedVec<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "AlignedVec(len: {}, ptr: {:p})", self.len, self.ptr)
    }
}

unsafe impl<T: Send> Send for AlignedVec<T> {}
unsafe impl<T: Sync> Sync for AlignedVec<T> {}

impl<T> Deref for AlignedVec<T> {
    type Target = [T];
    #[inline(always)]
    fn deref(&self) -> &[T] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl<T> DerefMut for AlignedVec<T> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [T] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}
