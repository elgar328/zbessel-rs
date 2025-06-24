use num_complex::Complex64;
use std::os::raw::{c_double, c_int};

// 생성된 바인딩을 포함합니다
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// 복소수 베셀 함수 결과를 나타내는 구조체
#[derive(Debug, Clone)]
pub struct BesselResult {
    /// 계산된 함수 값들
    pub values: Vec<Complex64>,
    /// 언더플로우가 발생한 함수 값의 개수
    pub underflow_count: i32,
}

/// 에러 유형
#[derive(Debug, Clone)]
pub enum BesselError {
    /// 잘못된 입력 매개변수
    InvalidParameter(String),
    /// 계산 오류
    ComputationError(String),
}

impl std::fmt::Display for BesselError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BesselError::InvalidParameter(msg) => write!(f, "잘못된 매개변수: {}", msg),
            BesselError::ComputationError(msg) => write!(f, "계산 오류: {}", msg),
        }
    }
}

impl std::error::Error for BesselError {}

/// 복소수 베셀 함수 J_ν(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(-abs(Im(z))) 스케일링)
/// * `n` - 계산할 함수 값의 개수
pub fn bessel_j(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n은 0보다 커야 합니다".to_string(),
        ));
    }

    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0i32;

    let result = unsafe {
        zbesj(
            z.re as c_double,
            z.im as c_double,
            nu as c_double,
            kode as c_int,
            n as c_int,
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zbesj 오류 코드: {}",
            result
        )));
    }

    let values = cyr
        .into_iter()
        .zip(cyi.into_iter())
        .map(|(r, i)| Complex64::new(r, i))
        .collect();

    Ok(BesselResult {
        values,
        underflow_count: nz,
    })
}

/// 복소수 베셀 함수 Y_ν(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(-abs(Im(z))) 스케일링)
/// * `n` - 계산할 함수 값의 개수
pub fn bessel_y(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n은 0보다 커야 합니다".to_string(),
        ));
    }

    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut cwrkr = vec![0.0; n];
    let mut cwrki = vec![0.0; n];
    let mut nz = 0i32;

    let result = unsafe {
        zbesy(
            z.re as c_double,
            z.im as c_double,
            nu as c_double,
            kode as c_int,
            n as c_int,
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
            cwrkr.as_mut_ptr(),
            cwrki.as_mut_ptr(),
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zbesy 오류 코드: {}",
            result
        )));
    }

    let values = cyr
        .into_iter()
        .zip(cyi.into_iter())
        .map(|(r, i)| Complex64::new(r, i))
        .collect();

    Ok(BesselResult {
        values,
        underflow_count: nz,
    })
}

/// 복소수 수정 베셀 함수 I_ν(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(-abs(Re(z))) 스케일링)
/// * `n` - 계산할 함수 값의 개수
pub fn bessel_i(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n은 0보다 커야 합니다".to_string(),
        ));
    }

    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0i32;

    let result = unsafe {
        zbesi(
            z.re as c_double,
            z.im as c_double,
            nu as c_double,
            kode as c_int,
            n as c_int,
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zbesi 오류 코드: {}",
            result
        )));
    }

    let values = cyr
        .into_iter()
        .zip(cyi.into_iter())
        .map(|(r, i)| Complex64::new(r, i))
        .collect();

    Ok(BesselResult {
        values,
        underflow_count: nz,
    })
}

/// 복소수 수정 베셀 함수 K_ν(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(z) 스케일링)
/// * `n` - 계산할 함수 값의 개수
pub fn bessel_k(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n은 0보다 커야 합니다".to_string(),
        ));
    }

    let mut cyr = vec![0.0; n];
    let mut cyi = vec![0.0; n];
    let mut nz = 0i32;

    let result = unsafe {
        zbesk(
            z.re as c_double,
            z.im as c_double,
            nu as c_double,
            kode as c_int,
            n as c_int,
            cyr.as_mut_ptr(),
            cyi.as_mut_ptr(),
            &mut nz,
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zbesk 오류 코드: {}",
            result
        )));
    }

    let values = cyr
        .into_iter()
        .zip(cyi.into_iter())
        .map(|(r, i)| Complex64::new(r, i))
        .collect();

    Ok(BesselResult {
        values,
        underflow_count: nz,
    })
}

/// 복소수 에어리 함수 Ai(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `id` - 미분 옵션 (0: Ai(z), 1: Ai'(z))
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(|Re(z)|*2/3) 스케일링)
pub fn airy_ai(z: Complex64, id: i32, kode: i32) -> Result<Complex64, BesselError> {
    let mut air = 0.0;
    let mut aii = 0.0;
    let mut nz = 0i32;

    let result = unsafe {
        zairy(
            z.re as c_double,
            z.im as c_double,
            id as c_int,
            kode as c_int,
            &mut air,
            &mut aii,
            &mut nz,
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zairy 오류 코드: {}",
            result
        )));
    }

    Ok(Complex64::new(air, aii))
}

/// 복소수 에어리 함수 Bi(z)를 계산합니다
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `id` - 미분 옵션 (0: Bi(z), 1: Bi'(z))
/// * `kode` - 스케일링 옵션 (1: 스케일링 안함, 2: exp(-|Re(z)|*2/3) 스케일링)
pub fn airy_bi(z: Complex64, id: i32, kode: i32) -> Result<Complex64, BesselError> {
    let mut bir = 0.0;
    let mut bii = 0.0;

    let result = unsafe {
        zbiry(
            z.re as c_double,
            z.im as c_double,
            id as c_int,
            kode as c_int,
            &mut bir,
            &mut bii,
        )
    };

    if result != 0 {
        return Err(BesselError::ComputationError(format!(
            "zbiry 오류 코드: {}",
            result
        )));
    }

    Ok(Complex64::new(bir, bii))
}

// ========================================
// 간단한 단일 값 계산 함수들
// ========================================

/// 베셀 함수 J_ν(z)를 계산합니다 (단일 값, 스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
///
/// # 반환값
/// J_ν(z)의 복소수 값
#[allow(non_snake_case)]
pub fn J(z: Complex64, nu: f64) -> Result<Complex64, BesselError> {
    let result = bessel_j(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// 베셀 함수 Y_ν(z)를 계산합니다 (단일 값, 스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
///
/// # 반환값
/// Y_ν(z)의 복소수 값
#[allow(non_snake_case)]
pub fn Y(z: Complex64, nu: f64) -> Result<Complex64, BesselError> {
    let result = bessel_y(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// 수정 베셀 함수 I_ν(z)를 계산합니다 (단일 값, 스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
///
/// # 반환값
/// I_ν(z)의 복소수 값
#[allow(non_snake_case)]
pub fn I(z: Complex64, nu: f64) -> Result<Complex64, BesselError> {
    let result = bessel_i(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// 수정 베셀 함수 K_ν(z)를 계산합니다 (단일 값, 스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
/// * `nu` - 차수 (실수)
///
/// # 반환값
/// K_ν(z)의 복소수 값
#[allow(non_snake_case)]
pub fn K(z: Complex64, nu: f64) -> Result<Complex64, BesselError> {
    let result = bessel_k(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// 에어리 함수 Ai(z)를 계산합니다 (스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
///
/// # 반환값
/// Ai(z)의 복소수 값
#[allow(non_snake_case)]
pub fn Ai(z: Complex64) -> Result<Complex64, BesselError> {
    airy_ai(z, 0, 1)
}

/// 에어리 함수 Bi(z)를 계산합니다 (스케일링 없음)
///
/// # 매개변수
/// * `z` - 복소수 인수
///
/// # 반환값
/// Bi(z)의 복소수 값
#[allow(non_snake_case)]
pub fn Bi(z: Complex64) -> Result<Complex64, BesselError> {
    airy_bi(z, 0, 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_j() {
        let z = Complex64::new(1.0, 0.5);
        let result = bessel_j(z, 0.0, 1, 1).unwrap();
        assert_eq!(result.values.len(), 1);
        // 기본적인 유효성 검사
        assert!(result.values[0].norm() > 0.0);
    }

    #[test]
    fn test_airy_ai() {
        let z = Complex64::new(1.0, 0.0);
        let result = airy_ai(z, 0, 1).unwrap();
        // Ai(1) ≈ 0.135... (실제 값과 비교)
        assert!((result.re - 0.135).abs() < 0.01);
    }

    // 간단한 API 테스트
    #[test]
    fn test_simple_j() {
        let z = Complex64::new(1.0, 0.5);
        let result = J(z, 0.0).unwrap();
        assert!(result.norm() > 0.0);
    }

    #[test]
    fn test_simple_i() {
        let z = Complex64::new(1.0, 0.5);
        let result = I(z, 0.0).unwrap();
        assert!(result.norm() > 0.0);
    }

    #[test]
    fn test_simple_ai() {
        let z = Complex64::new(1.0, 0.0);
        let result = Ai(z).unwrap();
        // Ai(1) ≈ 0.135...
        assert!((result.re - 0.135).abs() < 0.01);
        assert!(result.im.abs() < 1e-10); // 실수 입력에서는 허수부가 0에 가까워야 함
    }
}
