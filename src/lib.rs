//! # zbessel-rs
//!
//! A library that provides complex Bessel functions and Airy functions for Rust.
//!
//! ## Simple Usage
//!
//! ```rust
//! use num_complex::Complex64;
//! use zbessel_rs::{J, Y, I, K, Ai, Bi, J_scaled, Y_scaled, I_scaled, K_scaled, Ai_scaled, Bi_scaled};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let z = Complex64::new(1.0, 0.5);
//!     
//!     // Bessel functions (order, variable order)
//!     let j0 = J(0.0, z)?;  // J_0(z)
//!     let j1 = J(1.0, z)?;  // J_1(z)
//!     let y0 = Y(0.0, z)?;  // Y_0(z)
//!     let i0 = I(0.0, z)?;  // I_0(z)
//!     let k0 = K(0.0, z)?;  // K_0(z)
//!     
//!     // Airy functions
//!     let ai_val = Ai(z)?;  // Ai(z)
//!     let bi_val = Bi(z)?;  // Bi(z)
//!     
//!     // Scaled versions (useful for large arguments)
//!     let j0_scaled = J_scaled(0.0, z)?;  // J_0(z) with scaling
//!     let i0_scaled = I_scaled(0.0, z)?;  // I_0(z) with scaling
//!     let ai_scaled = Ai_scaled(z)?;      // Ai(z) with scaling
//!     
//!     println!("J_0({}) = {}", z, j0);
//!     println!("Y_0({}) = {}", z, y0);
//!     println!("I_0({}) = {}", z, i0);
//!     println!("K_0({}) = {}", z, k0);
//!     println!("Ai({}) = {}", z, ai_val);
//!     println!("Bi({}) = {}", z, bi_val);
//!     
//!     println!("J_0({}) (scaled) = {}", z, j0_scaled);
//!     println!("I_0({}) (scaled) = {}", z, i0_scaled);
//!     println!("Ai({}) (scaled) = {}", z, ai_scaled);
//!     
//!     Ok(())
//! }
//! ```
//!
//! ## Advanced Usage
//!
//! For multiple value calculations or scaling:
//!
//! ```rust
//! use num_complex::Complex64;
//! use zbessel_rs::{bessel_j, bessel_i, airy_ai};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let z = Complex64::new(2.0, 1.0);
//!     
//!     // Calculate multiple orders at once: J_0(z), J_1(z), J_2(z)
//!     let result = bessel_j(z, 0.0, 1, 3)?;
//!     for (n, value) in result.values.iter().enumerate() {
//!         println!("J_{}({}) = {}", n, z, value);
//!     }
//!     
//!     // Scaled result (kode=2)
//!     let scaled = bessel_i(z, 0.0, 2, 1)?;
//!     println!("I_0({}) (scaled) = {}", z, scaled.values[0]);
//!     
//!     // Airy function derivative (id=1)
//!     let ai_prime = airy_ai(z, 1, 1)?;
//!     println!("Ai'({}) = {}", z, ai_prime);
//!     
//!     Ok(())
//! }

use num_complex::Complex64;
use std::os::raw::{c_double, c_int};

// Include the generated bindings
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// Structure representing the result of complex Bessel function calculations
#[derive(Debug, Clone)]
pub struct BesselResult {
    /// Calculated function values
    pub values: Vec<Complex64>,
    /// Number of function values that experienced underflow
    pub underflow_count: i32,
}

/// Error types
#[derive(Debug, Clone)]
pub enum BesselError {
    /// Invalid input parameters
    InvalidParameter(String),
    /// Computation error
    ComputationError(String),
}

impl std::fmt::Display for BesselError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BesselError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
            BesselError::ComputationError(msg) => write!(f, "Computation error: {}", msg),
        }
    }
}

impl std::error::Error for BesselError {}

/// Calculate complex Bessel function J_ν(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `nu` - Order (real number)
/// * `kode` - Scaling option (1: no scaling, 2: exp(-abs(Im(z))) scaling)
/// * `n` - Number of function values to calculate
pub fn bessel_j(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n must be greater than 0".to_string(),
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
            "zbesj error code: {}",
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

/// Calculate complex Bessel function Y_ν(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `nu` - Order (real number)
/// * `kode` - Scaling option (1: no scaling, 2: exp(-abs(Im(z))) scaling)
/// * `n` - Number of function values to calculate
pub fn bessel_y(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n must be greater than 0".to_string(),
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
            "zbesy error code: {}",
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

/// Calculate complex modified Bessel function I_ν(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `nu` - Order (real number)
/// * `kode` - Scaling option (1: no scaling, 2: exp(-abs(Re(z))) scaling)
/// * `n` - Number of function values to calculate
pub fn bessel_i(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n must be greater than 0".to_string(),
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
            "zbesi error code: {}",
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

/// Calculate complex modified Bessel function K_ν(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `nu` - Order (real number)
/// * `kode` - Scaling option (1: no scaling, 2: exp(z) scaling)
/// * `n` - Number of function values to calculate
pub fn bessel_k(z: Complex64, nu: f64, kode: i32, n: usize) -> Result<BesselResult, BesselError> {
    if n == 0 {
        return Err(BesselError::InvalidParameter(
            "n must be greater than 0".to_string(),
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
            "zbesk error code: {}",
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

/// Calculate complex Airy function Ai(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `id` - Differentiation option (0: Ai(z), 1: Ai'(z))
/// * `kode` - Scaling option (1: no scaling, 2: exp(zeta) scaling where zeta=(2/3)*z^(3/2))
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
            "zairy error code: {}",
            result
        )));
    }

    Ok(Complex64::new(air, aii))
}

/// Calculate complex Airy function Bi(z)
///
/// # Parameters
/// * `z` - Complex argument
/// * `id` - Differentiation option (0: Bi(z), 1: Bi'(z))
/// * `kode` - Scaling option (1: no scaling, 2: exp(-abs(Re(zeta))) scaling where zeta=(2/3)*z^(3/2))
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
            "zbiry error code: {}",
            result
        )));
    }

    Ok(Complex64::new(bir, bii))
}

// ========================================
// Simple single-value calculation functions
// ========================================

/// Calculate Bessel function J_ν(z) (single value, no scaling)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of J_ν(z)
#[allow(non_snake_case)]
pub fn J(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_j(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// Calculate Bessel function Y_ν(z) (single value, no scaling)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Y_ν(z)
#[allow(non_snake_case)]
pub fn Y(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_y(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// Calculate modified Bessel function I_ν(z) (single value, no scaling)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of I_ν(z)
#[allow(non_snake_case)]
pub fn I(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_i(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// Calculate modified Bessel function K_ν(z) (single value, no scaling)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of K_ν(z)
#[allow(non_snake_case)]
pub fn K(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_k(z, nu, 1, 1)?;
    Ok(result.values[0])
}

/// Calculate Airy function Ai(z) (no scaling)
///
/// # Parameters
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Ai(z)
#[allow(non_snake_case)]
pub fn Ai(z: Complex64) -> Result<Complex64, BesselError> {
    airy_ai(z, 0, 1)
}

/// Calculate Airy function Bi(z) (no scaling)
///
/// # Parameters
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Bi(z)
#[allow(non_snake_case)]
pub fn Bi(z: Complex64) -> Result<Complex64, BesselError> {
    airy_bi(z, 0, 1)
}

// ========================================
// Scaled single-value calculation functions
// ========================================

/// Calculate Bessel function J_ν(z) with scaling (single value)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of J_ν(z) with exp(-abs(Im(z))) scaling
#[allow(non_snake_case)]
pub fn J_scaled(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_j(z, nu, 2, 1)?;
    Ok(result.values[0])
}

/// Calculate Bessel function Y_ν(z) with scaling (single value)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Y_ν(z) with exp(-abs(Im(z))) scaling
#[allow(non_snake_case)]
pub fn Y_scaled(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_y(z, nu, 2, 1)?;
    Ok(result.values[0])
}

/// Calculate modified Bessel function I_ν(z) with scaling (single value)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of I_ν(z) with exp(-abs(Re(z))) scaling
#[allow(non_snake_case)]
pub fn I_scaled(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_i(z, nu, 2, 1)?;
    Ok(result.values[0])
}

/// Calculate modified Bessel function K_ν(z) with scaling (single value)
///
/// # Parameters
/// * `nu` - Order (real number)
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of K_ν(z) with exp(z) scaling
#[allow(non_snake_case)]
pub fn K_scaled(nu: f64, z: Complex64) -> Result<Complex64, BesselError> {
    let result = bessel_k(z, nu, 2, 1)?;
    Ok(result.values[0])
}

/// Calculate Airy function Ai(z) with scaling
///
/// # Parameters
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Ai(z) with exp(zeta) scaling where zeta=(2/3)*z^(3/2)
#[allow(non_snake_case)]
pub fn Ai_scaled(z: Complex64) -> Result<Complex64, BesselError> {
    airy_ai(z, 0, 2)
}

/// Calculate Airy function Bi(z) with scaling
///
/// # Parameters
/// * `z` - Complex argument
///
/// # Returns
/// Complex value of Bi(z) with exp(-abs(Re(zeta))) scaling where zeta=(2/3)*z^(3/2)
#[allow(non_snake_case)]
pub fn Bi_scaled(z: Complex64) -> Result<Complex64, BesselError> {
    airy_bi(z, 0, 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Simple function tests
    #[test]
    fn test_simple_j() {
        let z = Complex64::new(10.0, 20.0);
        let nu = 1.0;
        let result = J(nu, z).unwrap();

        // Expected value
        let expected = Complex64::new(-1.3869150348751683e7, -3.785203660760742e7);

        let diff = (result - expected).norm();
        assert!(diff < 1e-8, "J test failed: diff = {}", diff);
    }

    #[test]
    fn test_simple_y() {
        let z = Complex64::new(10.0, 20.0);
        let nu = 1.0;
        let result = Y(nu, z).unwrap();

        // Expected value
        let expected = Complex64::new(3.785203660760742e7, -1.3869150348751685e7);

        let diff = (result - expected).norm();
        assert!(diff < 1e-8, "Y test failed: diff = {}", diff);
    }

    #[test]
    fn test_simple_i() {
        let z = Complex64::new(10.0, 20.0);
        let nu = 1.0;
        let result = I(nu, z).unwrap();

        // Expected value
        let expected = Complex64::new(1509.8283640825061, 1060.1232442308794);

        let diff = (result - expected).norm();
        assert!(diff < 1e-10, "I test failed: diff = {}", diff);
    }

    #[test]
    fn test_simple_k() {
        let z = Complex64::new(10.0, 20.0);
        let nu = 1.0;
        let result = K(nu, z).unwrap();

        // Expected value
        let expected = Complex64::new(-1.7871627759052974e-6, -1.1993686627062101e-5);

        let diff = (result - expected).norm();
        assert!(diff < 1e-10, "K test failed: diff = {}", diff);
    }

    #[test]
    fn test_simple_ai() {
        let z = Complex64::new(10.0, 20.0);
        let result = Ai(z).unwrap();

        // Expected value
        let expected = Complex64::new(14.71701664241453, -71.3378944410467);

        let diff = (result - expected).norm();
        assert!(diff < 1e-10, "Ai test failed: diff = {}", diff);
    }

    #[test]
    fn test_simple_bi() {
        let z = Complex64::new(10.0, 20.0);
        let result = Bi(z).unwrap();

        // Expected value
        let expected = Complex64::new(71.33821176573869, 14.71735250989695);

        let diff = (result - expected).norm();
        assert!(diff < 1e-10, "Bi test failed: diff = {}", diff);
    }

    #[test]
    fn test_j_scaling_consistency() {
        let z = Complex64::new(-100.0, 200.0);
        let nu = 1.0;

        // Calculate using regular J function
        let j_regular = J(nu, z).unwrap();

        // Calculate using J_scaled function
        let j_scaled = J_scaled(nu, z).unwrap();

        // For J functions, the scaling factor is exp(-abs(Im(z)))
        let scale_factor = (-z.im.abs()).exp();
        let j_regular_scaled = j_regular * scale_factor;

        // Check if J_scaled result matches J result multiplied by scale factor
        let diff = (j_scaled - j_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "J scaling consistency failed: diff = {}",
            diff
        );
    }

    #[test]
    fn test_y_scaling_consistency() {
        let z = Complex64::new(200.0, -150.0);
        let nu = 1.0;

        // Calculate using regular Y function
        let y_regular = Y(nu, z).unwrap();

        // Calculate using Y_scaled function
        let y_scaled = Y_scaled(nu, z).unwrap();

        // For Y functions, the scaling factor is exp(-abs(Im(z)))
        let scale_factor = (-z.im.abs()).exp();
        let y_regular_scaled = y_regular * scale_factor;

        // Check if Y_scaled result matches Y result multiplied by scale factor
        let diff = (y_scaled - y_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "Y scaling consistency failed: diff = {}",
            diff
        );
    }

    #[test]
    fn test_i_scaling_consistency() {
        let z = Complex64::new(200.0, -100.0);
        let nu = 1.0;

        // Calculate using regular I function
        let i_regular = I(nu, z).unwrap();

        // Calculate using I_scaled function
        let i_scaled = I_scaled(nu, z).unwrap();

        // For I functions, the scaling factor is exp(-abs(Re(z)))
        let scale_factor = (-z.re.abs()).exp();
        let i_regular_scaled = i_regular * scale_factor;

        // Check if I_scaled result matches I result multiplied by scale factor
        let diff = (i_scaled - i_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "I scaling consistency failed: diff = {}",
            diff
        );
    }

    #[test]
    fn test_k_scaling_consistency() {
        let z = Complex64::new(100.0, -50.0);
        let nu = 1.0;

        // Calculate using regular K function
        let k_regular = K(nu, z).unwrap();

        // Calculate using K_scaled function
        let k_scaled = K_scaled(nu, z).unwrap();

        // For K functions, the scaling factor is exp(z)
        let scale_factor = z.exp();
        let k_regular_scaled = k_regular * scale_factor;

        // Check if K_scaled result matches K result multiplied by scale factor
        let diff = (k_scaled - k_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "K scaling consistency failed: diff = {}",
            diff
        );
    }

    #[test]
    fn test_ai_scaling_consistency() {
        let z = Complex64::new(100.0, -50.0);

        // Calculate using regular Ai function
        let ai_regular = Ai(z).unwrap();

        // Calculate using Ai_scaled function
        let ai_scaled = Ai_scaled(z).unwrap();

        // For Ai functions, the scaling factor is exp(zeta) where zeta = (2/3)*z^(3/2)
        let z_sqrt = z.sqrt();
        let zeta = (2.0 / 3.0) * z * z_sqrt;
        let scale_factor = zeta.exp();
        let ai_regular_scaled = ai_regular * scale_factor;

        // Check if Ai_scaled result matches Ai result multiplied by scale factor
        let diff = (ai_scaled - ai_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "Ai scaling consistency failed: diff = {}",
            diff
        );
    }

    #[test]
    fn test_bi_scaling_consistency() {
        let z = Complex64::new(100.0, -50.0);

        // Calculate using regular Bi function
        let bi_regular = Bi(z).unwrap();

        // Calculate using Bi_scaled function
        let bi_scaled = Bi_scaled(z).unwrap();

        // For Bi functions, the scaling factor is exp(-abs(Re(zeta))) where zeta = (2/3)*z^(3/2)
        let z_sqrt = z.sqrt();
        let zeta = (2.0 / 3.0) * z * z_sqrt;
        let scale_factor = (-zeta.re.abs()).exp();
        let bi_regular_scaled = bi_regular * scale_factor;

        // Check if Bi_scaled result matches Bi result multiplied by scale factor
        let diff = (bi_scaled - bi_regular_scaled).norm();
        assert!(
            diff < 1e-10,
            "Bi scaling consistency failed: diff = {}",
            diff
        );
    }
}
