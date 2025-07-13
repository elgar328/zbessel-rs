# zbessel-rs

A Rust binding for the [zbessel](https://github.com/jpcima/zbessel) library that provides complex Bessel functions and Airy functions. This library is based on the original Fortran 77 implementation designed and implemented by D.E. Amos, which is part of the SLATEC mathematical library.

## Features

- **Complex Bessel Functions**: J_ν(z), Y_ν(z), I_ν(z), K_ν(z)
- **Complex Airy Functions**: Ai(z), Bi(z)
- **Scaled Functions**: All functions available with appropriate scaling factors
- **Safe Rust API**: Error handling using Result types
- **Auto Build**: Automatic C binding generation using bindgen
- **Thread-safe**: Based on the original library's stateless design
- **No Runtime Dependencies**: No f2c or gfortran runtime dependencies

## Supported Platforms

- **macOS** (x86_64, aarch64)
- **Windows** (MSVC)
- **Linux** (x86_64, aarch64)

## Installation

Add the following to your `Cargo.toml`:

```toml
[dependencies]
zbessel-rs = "0.1"
num-complex = "^0.4"
```

## Usage

### Simple API (Recommended)

Simple functions for the most common case of single value calculations:

```rust
use num_complex::Complex64;
use zbessel_rs::{J, Y, I, K, Ai, Bi, J_scaled, Y_scaled, I_scaled, K_scaled, Ai_scaled, Bi_scaled};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let z = Complex64::new(1.0, 0.5);
    
    // Bessel functions (no scaling)
    let j0 = J(0.0, z)?;  // J_0(z)
    let j1 = J(1.0, z)?;  // J_1(z)
    let y0 = Y(0.0, z)?;  // Y_0(z)
    let i0 = I(0.0, z)?;  // I_0(z)
    let k0 = K(0.0, z)?;  // K_0(z)
    
    // Airy functions (no scaling)
    let ai_val = Ai(z)?;  // Ai(z)
    let bi_val = Bi(z)?;  // Bi(z)
    
    // Scaled versions (useful for large arguments)
    let j0_scaled = J_scaled(0.0, z)?;  // J_0(z) with exp(-|Im(z)|) scaling
    let y0_scaled = Y_scaled(0.0, z)?;  // Y_0(z) with exp(-|Im(z)|) scaling
    let i0_scaled = I_scaled(0.0, z)?;  // I_0(z) with exp(-|Re(z)|) scaling
    let k0_scaled = K_scaled(0.0, z)?;  // K_0(z) with exp(z) scaling
    let ai_scaled = Ai_scaled(z)?;      // Ai(z) with exp(zeta) scaling, where zeta = (2/3)*z^(3/2)
    let bi_scaled = Bi_scaled(z)?;      // Bi(z) with exp(-|Re(zeta)|) scaling, where zeta = (2/3)*z^(3/2)
    
    println!("J_0({}) = {}", z, j0);
    println!("Y_0({}) = {}", z, y0);
    println!("I_0({}) = {}", z, i0);
    println!("K_0({}) = {}", z, k0);
    println!("Ai({}) = {}", z, ai_val);
    println!("Bi({}) = {}", z, bi_val);
    
    println!("J_0({}) (scaled) = {}", z, j0_scaled);
    println!("Y_0({}) (scaled) = {}", z, y0_scaled);
    println!("I_0({}) (scaled) = {}", z, i0_scaled);
    println!("K_0({}) (scaled) = {}", z, k0_scaled);
    println!("Ai({}) (scaled) = {}", z, ai_scaled);
    println!("Bi({}) (scaled) = {}", z, bi_scaled);
    
    Ok(())
}
```

### Low-level API

> **Low-level API**
> These functions directly expose the original AMOS library interface for advanced use cases.
> You can specify the scaling option and compute multiple orders at once, just like the underlying Fortran/C routines.
> This API is recommended for users who need full control or want to match AMOS documentation and behavior exactly.

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

### Scaled Simple API

#### `J_scaled(nu, z) -> Result<Complex64, BesselError>`
Calculate Bessel function J_ν(z) with scaling.

- `nu`: Order (real number)
- `z`: Complex argument
- **Scaling factor**: `exp(-|Im(z)|)`

#### `Y_scaled(nu, z) -> Result<Complex64, BesselError>`
Calculate Bessel function Y_ν(z) with scaling.

- `nu`: Order (real number)
- `z`: Complex argument
- **Scaling factor**: `exp(-|Im(z)|)`

#### `I_scaled(nu, z) -> Result<Complex64, BesselError>`
Calculate modified Bessel function I_ν(z) with scaling.

- `nu`: Order (real number)
- `z`: Complex argument
- **Scaling factor**: `exp(-|Re(z)|)`

#### `K_scaled(nu, z) -> Result<Complex64, BesselError>`
Calculate modified Bessel function K_ν(z) with scaling.

- `nu`: Order (real number)
- `z`: Complex argument
- **Scaling factor**: `exp(z)`

#### `Ai_scaled(z) -> Result<Complex64, BesselError>`
Calculate Airy function Ai(z) with scaling.

- `z`: Complex argument
- **Scaling factor**: `exp(zeta)` where `zeta = (2/3)*z^(3/2)`

#### `Bi_scaled(z) -> Result<Complex64, BesselError>`
Calculate Airy function Bi(z) with scaling.

- `z`: Complex argument
- **Scaling factor**: `exp(-|Re(zeta)|)` where `zeta = (2/3)*z^(3/2)`

### Low-level API

#### `bessel_j(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex Bessel function J_ν(z).

- `z`: Complex argument
- `nu`: Order (real number)
- `kode`: Scaling option (1: no scaling, 2: exp(-|Im(z)|) scaling)
- `n`: Number of function values to calculate

#### `bessel_y(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex Bessel function Y_ν(z).

- `z`: Complex argument
- `nu`: Order (real number)
- `kode`: Scaling option (1: no scaling, 2: exp(-|Im(z)|) scaling)
- `n`: Number of function values to calculate

#### `bessel_i(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex modified Bessel function I_ν(z).

- `z`: Complex argument
- `nu`: Order (real number)
- `kode`: Scaling option (1: no scaling, 2: exp(-|Re(z)|) scaling)
- `n`: Number of function values to calculate

#### `bessel_k(z, nu, kode, n) -> Result<BesselResult, BesselError>`
Calculate complex modified Bessel function K_ν(z).

- `z`: Complex argument
- `nu`: Order (real number)
- `kode`: Scaling option (1: no scaling, 2: exp(z) scaling)
- `n`: Number of function values to calculate

### Airy Functions

#### `airy_ai(z, id, kode) -> Result<Complex64, BesselError>`
Calculate complex Airy function Ai(z).

- `z`: Complex argument
- `id`: Differentiation option (0: Ai(z), 1: Ai'(z))
- `kode`: Scaling option (1: no scaling, 2: exp(zeta) scaling where zeta=(2/3)*z^(3/2))

#### `airy_bi(z, id, kode) -> Result<Complex64, BesselError>`
Calculate complex Airy function Bi(z).

- `z`: Complex argument
- `id`: Differentiation option (0: Bi(z), 1: Bi'(z))
- `kode`: Scaling option (1: no scaling, 2: exp(-|Re(zeta)|) scaling where zeta=(2/3)*z^(3/2))

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 