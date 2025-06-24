# zbessel-rs

A library that provides complex Bessel functions and Airy functions for Rust. This library is a Rust binding for the C++ zbessel library, which is based on the original Fortran 77 implementation designed and implemented by D.E. Amos.

## Key Features

- **Complex Bessel Functions**: J_ν(z), Y_ν(z), I_ν(z), K_ν(z)
- **Complex Airy Functions**: Ai(z), Bi(z)
- **Safe Rust API**: Error handling using Result types
- **Auto Build**: Automatic C binding generation using bindgen

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
zbessel-rs = { git = "https://github.com/elgar328/zbessel-rs.git" }
num-complex = "0.4"
```

## Usage

### Simple API (Recommended)

Simple functions for the most common case of single value calculations:

```rust
use num_complex::Complex64;
use zbessel_rs::{J, Y, I, K, Ai, Bi};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let z = Complex64::new(1.0, 0.5);
    
    // Bessel functions
    let j0 = J(0.0, z)?;  // J_0(z)
    let j1 = J(1.0, z)?;  // J_1(z)
    let y0 = Y(0.0, z)?;  // Y_0(z)
    let i0 = I(0.0, z)?;  // I_0(z)
    let k0 = K(0.0, z)?;  // K_0(z)
    
    // Airy functions
    let ai_val = Ai(z)?;  // Ai(z)
    let bi_val = Bi(z)?;  // Bi(z)
    
    println!("J_0({}) = {}", z, j0);
    println!("Y_0({}) = {}", z, y0);
    println!("I_0({}) = {}", z, i0);
    println!("K_0({}) = {}", z, k0);
    println!("Ai({}) = {}", z, ai_val);
    println!("Bi({}) = {}", z, bi_val);
    
    Ok(())
}
```

### Advanced API

For multiple value calculations or scaling:

```rust
use num_complex::Complex64;
use zbessel_rs::{bessel_j, bessel_i, airy_ai};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let z = Complex64::new(2.0, 1.0);
    
    // Calculate multiple orders at once: J_0(z), J_1(z), J_2(z)
    let result = bessel_j(z, 0.0, 1, 3)?;
    for (n, value) in result.values.iter().enumerate() {
        println!("J_{}({}) = {}", n, z, value);
    }
    
    // Scaled result (kode=2)
    let scaled = bessel_i(z, 0.0, 2, 1)?;
    println!("I_0({}) (scaled) = {}", z, scaled.values[0]);
    
    // Airy function derivative (id=1)
    let ai_prime = airy_ai(z, 1, 1)?;
    println!("Ai'({}) = {}", z, ai_prime);
    
    Ok(())
}
```

## API Reference

### Simple API

#### `J(nu, z) -> Result<Complex64, BesselError>`
Calculate Bessel function J_ν(z) (single value, no scaling).

- `nu`: Order (real number)
- `z`: Complex argument

#### `Y(nu, z) -> Result<Complex64, BesselError>`
Calculate Bessel function Y_ν(z) (single value, no scaling).

#### `I(nu, z) -> Result<Complex64, BesselError>`
Calculate modified Bessel function I_ν(z) (single value, no scaling).

#### `K(nu, z) -> Result<Complex64, BesselError>`
Calculate modified Bessel function K_ν(z) (single value, no scaling).

#### `Ai(z) -> Result<Complex64, BesselError>`
Calculate Airy function Ai(z) (no scaling).

#### `Bi(z) -> Result<Complex64, BesselError>`
Calculate Airy function Bi(z) (no scaling).

### Advanced API

#### `bessel_j(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex Bessel function J_ν(z).

- `z`: Complex argument
- `nu`: Order (real number)
- `kode`: Scaling option (1: no scaling, 2: exp(-abs(Im(z))) scaling)
- `n`: Number of function values to calculate

#### `bessel_y(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex Bessel function Y_ν(z).

#### `bessel_i(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex modified Bessel function I_ν(z).

#### `bessel_k(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex modified Bessel function K_ν(z).

### Airy Functions

#### `airy_ai(z, id, kode) -> Result<Complex64, BesselError>`
Calculate complex Airy function Ai(z).

- `z`: Complex argument
- `id`: Differentiation option (0: Ai(z), 1: Ai'(z))
- `kode`: Scaling option (1: no scaling, 2: exp(|Re(z)|*2/3) scaling)

#### `airy_bi(z, id, kode) -> Result<Complex64, BesselError>`
Calculate complex Airy function Bi(z).

## Build and Run

### Build Library
```bash
cargo build
```

### Run Tests
```bash
cargo test
```

## Data Types

### `BesselResult`
```rust
pub struct BesselResult {
    /// Calculated function values
    pub values: Vec<Complex64>,
    /// Number of function values that experienced underflow
    pub underflow_count: i32,
}
```

### `BesselError`
```rust
pub enum BesselError {
    /// Invalid input parameters
    InvalidParameter(String),
    /// Computation error
    ComputationError(String),
}
```

## Original Library Information

This Rust binding is based on:
- **Original Implementation**: D.E. Amos's Fortran 77 implementation (part of the SLATEC mathematical library)
- **C++ Port**: Modern C++ implementation based on f2c conversion
- **Features**: Thread-safe, no f2c/gfortran runtime dependencies

## License

Follows the license of the original zbessel library. 