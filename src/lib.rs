//! # zbessel-rs
//!
//! A library that provides complex Bessel functions and Airy functions for Rust.
//!
//! ## Simple Usage
//!
//! ```rust
//! use num_complex::Complex64;
//! use zbessel_rs::{J, Y, I, K, Ai, Bi};
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
//!     println!("J_0({}) = {}", z, j0);
//!     println!("Y_0({}) = {}", z, y0);
//!     println!("I_0({}) = {}", z, i0);
//!     println!("K_0({}) = {}", z, k0);
//!     println!("Ai({}) = {}", z, ai_val);
//!     println!("Bi({}) = {}", z, bi_val);
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
/// * `kode` - Scaling option (1: no scaling, 2: exp(|Re(z)|*2/3) scaling)
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
/// * `kode` - Scaling option (1: no scaling, 2: exp(-|Re(z)|*2/3) scaling)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_j() {
        let z = Complex64::new(1.0, 0.5);
        let result = bessel_j(z, 0.0, 1, 1).unwrap();
        assert_eq!(result.values.len(), 1);
        // Basic validity check
        assert!(result.values[0].norm() > 0.0);
    }

    #[test]
    fn test_airy_ai() {
        let z = Complex64::new(1.0, 0.0);
        let result = airy_ai(z, 0, 1).unwrap();
        // Ai(1) ≈ 0.135... (compare with actual value)
        assert!((result.re - 0.135).abs() < 0.01);
    }

    // Simple API tests
    #[test]
    fn test_simple_j() {
        let z = Complex64::new(1.0, 0.5);
        let result = J(0.0, z).unwrap();
        assert!(result.norm() > 0.0);
    }

    #[test]
    fn test_simple_i() {
        let z = Complex64::new(1.0, 0.5);
        let result = I(0.0, z).unwrap();
        assert!(result.norm() > 0.0);
    }

    #[test]
    fn test_simple_ai() {
        let z = Complex64::new(1.0, 0.0);
        let result = Ai(z).unwrap();
        // Ai(1) ≈ 0.135...
        assert!((result.re - 0.135).abs() < 0.01);
        assert!(result.im.abs() < 1e-10); // For real input, imaginary part should be close to zero
    }
}
