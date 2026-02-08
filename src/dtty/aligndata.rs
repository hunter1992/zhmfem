use crate::dtty::basic::Dtype;
use std::arch::x86_64::*;
use std::fmt;
use std::ops::{Deref, DerefMut};

pub struct AlignedF64 {
    pub ptr: *mut f64,
    pub len: usize,
    pub layout: std::alloc::Layout,
}

/// A 64-byte memory-aligned buffer
impl AlignedF64 {
    const ALIGN: usize = 64; //byte size to be aligned

    /// All elements initialized to 0
    pub fn new(len: usize) -> Self {
        if len == 0 {
            return Self::empty();
        }

        // 计算总字节数并设定对齐，f64 占用 8 字节
        let layout = std::alloc::Layout::from_size_align(len * 8, Self::ALIGN)
            .expect("!!! 无效的内存布局设定！error from: src/dtty/aligndata");

        // 分配内存
        // 使用 alloc_zeroed 确保内存是干净的（对有限元组装很关键）
        // 如果你确定马上就会用 fast_convert 覆盖它，也可以改用 alloc(layout)
        // let ptr = unsafe { alloc(layout) as *mut f64 };
        let ptr = unsafe { std::alloc::alloc_zeroed(layout) as *mut f64 };

        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }

        Self { ptr, len, layout }
    }

    /// All elements uninitialized
    pub unsafe fn new_uninit(len: usize) -> Self {
        if len == 0 {
            return Self::empty();
        }

        // 计算总字节数并设定对齐，f64 占用 8 字节
        let layout = std::alloc::Layout::from_size_align(len * 8, Self::ALIGN)
            .expect("!!! 无效的内存布局设定！error from: src/dtty/aligndata");

        let ptr = unsafe { std::alloc::alloc(layout) as *mut f64 };

        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }

        Self { ptr, len, layout }
    }

    /// Return an empty AlignedF64 when arg len is 0
    fn empty() -> Self {
        Self {
            ptr: std::ptr::null_mut(),
            len: 0,
            layout: std::alloc::Layout::from_size_align(0, Self::ALIGN).unwrap(),
        }
    }

    /// Construct from Vec<Dtype>
    #[target_feature(enable = "avx")]
    pub fn from_slice(src: &[Dtype]) -> Self {
        let len = src.len();
        if len == 0 {
            return Self::empty();
        }

        let layout = std::alloc::Layout::from_size_align(len * 8, Self::ALIGN)
            .expect("!!! Invalid layout! error from src/dtty/aligndata: from_vec func");

        unsafe {
            let dst_ptr = std::alloc::alloc(layout) as *mut f64;
            if dst_ptr.is_null() {
                std::alloc::handle_alloc_error(layout);
            }

            if std::mem::size_of::<Dtype>() == 8 {
                // case 1: f64 -> f64
                std::ptr::copy_nonoverlapping(src.as_ptr() as *const f64, dst_ptr, len);
            } else {
                // case 2: f32 -> f64
                let src_const_ptr = src.as_ptr() as *const f32;
                let mut i = 0;
                while i + 15 < len {
                    let s_ptr = src_const_ptr.add(i);
                    let d_ptr = dst_ptr.add(i);

                    // 加载与转换 (利用流水线并行)
                    let v0 = _mm256_cvtps_pd(_mm_loadu_ps(s_ptr));
                    let v1 = _mm256_cvtps_pd(_mm_loadu_ps(s_ptr.add(4)));
                    let v2 = _mm256_cvtps_pd(_mm_loadu_ps(s_ptr.add(8)));
                    let v3 = _mm256_cvtps_pd(_mm_loadu_ps(s_ptr.add(12)));

                    // 对齐存储 (Aligned Store)
                    _mm256_store_pd(d_ptr, v0);
                    _mm256_store_pd(d_ptr.add(4), v1);
                    _mm256_store_pd(d_ptr.add(8), v2);
                    _mm256_store_pd(d_ptr.add(12), v3);

                    i += 16;
                }

                // 剩余部分的 SIMD 处理 (每次 4 个)
                while i + 3 < len {
                    _mm256_store_pd(
                        dst_ptr.add(i),
                        _mm256_cvtps_pd(_mm_loadu_ps(src_const_ptr.add(i))),
                    );
                    i += 4;
                }

                // 标量尾部处理
                while i < len {
                    *dst_ptr.add(i) = *src_const_ptr.add(i) as f64;
                    i += 1;
                }
            }

            AlignedF64 {
                ptr: dst_ptr,
                len,
                layout,
            }
        }
    }

    /// 将对齐的 f64 数据转换回普通的 Vec<Dtype>
    #[target_feature(enable = "avx")]
    pub fn into_vec(self) -> Vec<Dtype> {
        let len = self.len;
        let mut result = Vec::with_capacity(len);

        unsafe {
            let src_ptr = self.ptr;
            let dst_ptr = result.as_mut_ptr() as *mut Dtype;

            if std::mem::size_of::<Dtype>() == 8 {
                // 情况 1: f64 -> f64 (极致拷贝)
                std::ptr::copy_nonoverlapping(src_ptr as *const f64, dst_ptr as *mut f64, len);
            } else {
                // 情况 2: f64 -> f32 (AVX 循环展开窄化)
                let dst = dst_ptr as *mut f32;
                let mut i = 0;
                // 4x 循环展开：每个循环处理 16 个元素 (使用 4 个 AVX 寄存器)
                while i + 15 < len {
                    // 加载 4 个 256-bit f64 块
                    // 转换并将 4 个 f32 压缩进 128-bit 块
                    let v0 = _mm256_cvtpd_ps(_mm256_load_pd(src_ptr.add(i)));
                    let v1 = _mm256_cvtpd_ps(_mm256_load_pd(src_ptr.add(i + 4)));
                    let v2 = _mm256_cvtpd_ps(_mm256_load_pd(src_ptr.add(i + 8)));
                    let v3 = _mm256_cvtpd_ps(_mm256_load_pd(src_ptr.add(i + 12)));

                    // 存储到 f32 目标内存 (unaligned store，因为 Vec<f32> 只有 4 字节对齐)
                    _mm_storeu_ps(dst.add(i), v0);
                    _mm_storeu_ps(dst.add(i + 4), v1);
                    _mm_storeu_ps(dst.add(i + 8), v2);
                    _mm_storeu_ps(dst.add(i + 12), v3);

                    i += 16;
                }

                // 剩余部分的单步 SIMD
                while i + 3 < len {
                    _mm_storeu_ps(dst.add(i), _mm256_cvtpd_ps(_mm256_load_pd(src_ptr.add(i))));
                    i += 4;
                }

                // 标量尾部
                while i < len {
                    *dst.add(i) = *src_ptr.add(i) as f32;
                    i += 1;
                }
            }

            // 手动设置 Vec 的长度，因为是直接写内存的
            result.set_len(len);
        }
        result
        // 此处 self 离开作用域，会触发 Drop 释放对齐内存
    }

    pub fn as_ptr(&self) -> *mut f64 {
        self.ptr
    }

    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut f64 {
        self.ptr
    }

    pub fn clear(&self) {
        unsafe {
            std::ptr::write_bytes(self.as_ptr(), 0, self.len());
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }
}

impl Deref for AlignedF64 {
    type Target = [f64];
    fn deref(&self) -> &[f64] {
        if self.ptr.is_null() {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
        }
    }
}

impl DerefMut for AlignedF64 {
    fn deref_mut(&mut self) -> &mut [f64] {
        if self.ptr.is_null() {
            &mut []
        } else {
            unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
        }
    }
}

impl Drop for AlignedF64 {
    fn drop(&mut self) {
        if !self.ptr.is_null() && self.layout.size() > 0 {
            unsafe {
                std::alloc::dealloc(self.ptr as *mut u8, self.layout);
            }
        }
    }
}

impl fmt::Display for AlignedF64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let len = self.len();
        if len == 0 {
            return write!(f, "[] (empty)");
        }

        let slice = &self[..]; // 利用 Deref 特性获取切片
        let display_limit = 16;

        write!(f, "[")?;

        // 打印逻辑
        if len <= display_limit {
            // 如果长度小于等于 16，全量打印
            for (i, val) in slice.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{:.6e}", val)?; // 使用科学计数法，保留 6 位小数
            }
        } else {
            // 如果长度超过 16，只打印前 16 个
            for i in 0..display_limit {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{:.6e}", slice[i])?;
            }
            write!(f, ", ...")?;
        }

        // 打印后缀信息，标注长度和 64 字节对齐
        write!(f, "] (len: {}, aligned: 64B)", len)
    }
}

impl fmt::Debug for AlignedF64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let display_len = self.len.min(16); //只查看前16个
        f.debug_struct("AlignedF64")
            .field("address", &self.ptr)
            .field("len", &self.len)
            .field("data_preview", &&self[..display_len])
            .finish()
    }
}

pub struct AlignedI32 {
    pub ptr: *mut i32,
    pub len: usize,
    pub layout: std::alloc::Layout,
}

impl AlignedI32 {
    const ALIGN: usize = 64;

    /// 从切片构造 AlignedI32
    pub fn from_slice(src: &[i32]) -> Self {
        let len = src.len();
        if len == 0 {
            return Self::empty();
        }

        // i32 占 4 字节
        let size = len * 4;
        let layout =
            std::alloc::Layout::from_size_align(size, Self::ALIGN).expect("Invalid layout");

        unsafe {
            let dst_ptr = std::alloc::alloc(layout) as *mut i32;
            if dst_ptr.is_null() {
                std::alloc::handle_alloc_error(layout);
            }

            // 直接使用系统级极致拷贝
            // 对于 i32 -> i32，编译器会根据 i5-8265U 的特性自动选择最优的 AVX2 拷贝指令
            std::ptr::copy_nonoverlapping(src.as_ptr(), dst_ptr, len);

            AlignedI32 {
                ptr: dst_ptr,
                len,
                layout,
            }
        }
    }

    pub fn into_vec_i32(self) -> Vec<i32> {
        let len = self.len;
        if len == 0 {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(len);

        unsafe {
            let dst_ptr = result.as_mut_ptr();
            let src_ptr = self.ptr;

            std::ptr::copy_nonoverlapping(src_ptr, dst_ptr, len);

            result.set_len(len);
        }

        result
    }

    pub fn into_vec_usize(self) -> Vec<usize> {
        let len = self.len;
        if len == 0 {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(len);

        unsafe {
            let src_ptr = self.ptr;
            let dst_ptr = result.as_mut_ptr() as *mut i64; // usize 在 64 位下等同于 u64/i64
            let mut i = 0;

            #[cfg(target_feature = "avx2")]
            {
                // 4x 循环展开：每次处理 16 个元素 (4 个 AVX 寄存器)
                while i + 15 < len {
                    // _mm_loadu_si128 加载 4 个 i32 (128位)
                    // _mm256_cvtepi32_epi64 拓宽为 4 个 i64 (256位)
                    let v0 =
                        _mm256_cvtepi32_epi64(_mm_loadu_si128(src_ptr.add(i) as *const __m128i));
                    let v1 = _mm256_cvtepi32_epi64(_mm_loadu_si128(
                        src_ptr.add(i + 4) as *const __m128i
                    ));
                    let v2 = _mm256_cvtepi32_epi64(_mm_loadu_si128(
                        src_ptr.add(i + 8) as *const __m128i
                    ));
                    let v3 = _mm256_cvtepi32_epi64(_mm_loadu_si128(
                        src_ptr.add(i + 12) as *const __m128i
                    ));

                    // 存储到 Vec 的连续空间
                    _mm256_storeu_si256(dst_ptr.add(i) as *mut __m256i, v0);
                    _mm256_storeu_si256(dst_ptr.add(i + 4) as *mut __m256i, v1);
                    _mm256_storeu_si256(dst_ptr.add(i + 8) as *mut __m256i, v2);
                    _mm256_storeu_si256(dst_ptr.add(i + 12) as *mut __m256i, v3);

                    i += 16;
                }

                // 处理剩余的 4 的倍数
                while i + 3 < len {
                    let v =
                        _mm256_cvtepi32_epi64(_mm_loadu_si128(src_ptr.add(i) as *const __m128i));
                    _mm256_storeu_si256(dst_ptr.add(i) as *mut __m256i, v);
                    i += 4;
                }
            }

            // 2. 标量尾部处理
            while i < len {
                // 直接转换，由于是正索引，as usize 是安全的
                *dst_ptr.add(i) = *src_ptr.add(i) as i64;
                i += 1;
            }

            // 3. 设置 Vec 正确长度
            result.set_len(len);
        }

        result
    }

    pub fn empty() -> Self {
        Self {
            ptr: std::ptr::null_mut(),
            len: 0,
            layout: std::alloc::Layout::from_size_align(0, Self::ALIGN).unwrap(),
        }
    }

    /// 专门为 PARDISO 准备的指针获取方法
    #[inline]
    pub fn as_ptr(&self) -> *const i32 {
        self.ptr
    }
    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut i32 {
        self.ptr
    }
}

impl Drop for AlignedI32 {
    fn drop(&mut self) {
        if !self.ptr.is_null() && self.layout.size() > 0 {
            unsafe {
                std::alloc::dealloc(self.ptr as *mut u8, self.layout);
            }
        }
    }
}

impl Deref for AlignedI32 {
    type Target = [i32];
    fn deref(&self) -> &[i32] {
        if self.ptr.is_null() {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
        }
    }
}

impl DerefMut for AlignedI32 {
    fn deref_mut(&mut self) -> &mut [i32] {
        if self.ptr.is_null() {
            &mut []
        } else {
            unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
        }
    }
}

impl fmt::Display for AlignedI32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let len = self.len();
        let slice = &self[..]; // 利用之前实现的 Deref

        if len == 0 {
            return write!(f, "[] (empty)");
        }

        // 定义打印阈值，超过 16 个元素就进行缩略
        let threshold = 16;

        write!(f, "[")?;
        if len <= threshold {
            for (i, val) in slice.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", val)?;
            }
        } else {
            // 打印前 5 个
            for i in 0..5 {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", slice[i])?;
            }
            write!(f, ", ... ")?;
            // 打印后 5 个
            for i in (len - 5)..len {
                write!(f, "{}", slice[i])?;
                if i < len - 1 {
                    write!(f, ", ")?;
                }
            }
        }
        write!(f, "] (len: {}, aligned: 64B)", len)
    }
}

impl fmt::Debug for AlignedI32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // 调试模式下可以复用 Display 的逻辑，或者增加更多原始指针信息
        write!(
            f,
            "AlignedI32 {{ ptr: {:?}, len: {} }}",
            self.as_ptr(),
            self.len()
        )?;
        fmt::Display::fmt(self, f)
    }
}
