# zbessel-rs

Rust에서 복소수 베셀 함수와 에어리 함수를 사용할 수 있게 해주는 라이브러리입니다. 이 라이브러리는 D.E. Amos가 설계하고 구현한 원본 Fortran 77 구현을 기반으로 한 C++ zbessel 라이브러리의 Rust 바인딩입니다.

## 주요 기능

- **복소수 베셀 함수**: J_ν(z), Y_ν(z), I_ν(z), K_ν(z)
- **복소수 에어리 함수**: Ai(z), Bi(z)
- **안전한 Rust API**: Result 타입을 사용한 오류 처리
- **자동 빌드**: bindgen을 사용한 자동 C 바인딩 생성

## 설치

`Cargo.toml`에 다음을 추가하세요:

```toml
[dependencies]
zbessel-rs = { path = "." }
num-complex = "0.4"
```

## 사용법

### 간단한 API (추천)

가장 일반적인 경우인 단일 값 계산을 위한 간단한 함수들:

```rust
use num_complex::Complex64;
use zbessel_rs::{J, Y, I, K, Ai, Bi};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let z = Complex64::new(1.0, 0.5);
    
    // 베셀 함수들
    let j0 = J(0.0, z)?;  // J_0(z)
    let j1 = J(1.0, z)?;  // J_1(z)
    let y0 = Y(0.0, z)?;  // Y_0(z)
    let i0 = I(0.0, z)?;  // I_0(z)
    let k0 = K(0.0, z)?;  // K_0(z)
    
    // 에어리 함수들
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

### 고급 API

다중 값 계산이나 스케일링이 필요한 경우:

```rust
use num_complex::Complex64;
use zbessel_rs::{bessel_j, bessel_i, airy_ai};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let z = Complex64::new(2.0, 1.0);
    
    // 여러 차수를 한 번에 계산: J_0(z), J_1(z), J_2(z)
    let result = bessel_j(z, 0.0, 1, 3)?;
    for (n, value) in result.values.iter().enumerate() {
        println!("J_{}({}) = {}", n, z, value);
    }
    
    // 스케일링된 결과 (kode=2)
    let scaled = bessel_i(z, 0.0, 2, 1)?;
    println!("I_0({}) (scaled) = {}", z, scaled.values[0]);
    
    // 에어리 함수의 도함수 (id=1)
    let ai_prime = airy_ai(z, 1, 1)?;
    println!("Ai'({}) = {}", z, ai_prime);
    
    Ok(())
}
```

## API 참조

### 간단한 API

#### `J(nu, z) -> Result<Complex64, BesselError>`
베셀 함수 J_ν(z)를 계산합니다 (단일 값, 스케일링 없음).

- `nu`: 차수 (실수)
- `z`: 복소수 인수

#### `Y(nu, z) -> Result<Complex64, BesselError>`
베셀 함수 Y_ν(z)를 계산합니다 (단일 값, 스케일링 없음).

#### `I(nu, z) -> Result<Complex64, BesselError>`
수정 베셀 함수 I_ν(z)를 계산합니다 (단일 값, 스케일링 없음).

#### `K(nu, z) -> Result<Complex64, BesselError>`
수정 베셀 함수 K_ν(z)를 계산합니다 (단일 값, 스케일링 없음).

#### `Ai(z) -> Result<Complex64, BesselError>`
에어리 함수 Ai(z)를 계산합니다 (스케일링 없음).

#### `Bi(z) -> Result<Complex64, BesselError>`
에어리 함수 Bi(z)를 계산합니다 (스케일링 없음).

### 고급 API

#### `bessel_j(z, nu, kode, n) -> Result<BesselResult, BesselError>`
복소수 베셀 함수 J_ν(z)를 계산합니다.

- `z`: 복소수 인수
- `nu`: 차수 (실수)
- `kode`: 스케일링 옵션 (1: 스케일링 안함, 2: exp(-abs(Im(z))) 스케일링)
- `n`: 계산할 함수 값의 개수

#### `bessel_y(z, nu, kode, n) -> Result<BesselResult, BesselError>`
복소수 베셀 함수 Y_ν(z)를 계산합니다.

#### `bessel_i(z, nu, kode, n) -> Result<BesselResult, BesselError>`
복소수 수정 베셀 함수 I_ν(z)를 계산합니다.

#### `bessel_k(z, nu, kode, n) -> Result<BesselResult, BesselError>`
복소수 수정 베셀 함수 K_ν(z)를 계산합니다.

### 에어리 함수

#### `airy_ai(z, id, kode) -> Result<Complex64, BesselError>`
복소수 에어리 함수 Ai(z)를 계산합니다.

- `z`: 복소수 인수
- `id`: 미분 옵션 (0: Ai(z), 1: Ai'(z))
- `kode`: 스케일링 옵션 (1: 스케일링 안함, 2: exp(|Re(z)|*2/3) 스케일링)

#### `airy_bi(z, id, kode) -> Result<Complex64, BesselError>`
복소수 에어리 함수 Bi(z)를 계산합니다.

## 빌드 및 실행

### 라이브러리 빌드
```bash
cargo build
```

### 테스트 실행
```bash
cargo test
```

## 데이터 타입

### `BesselResult`
```rust
pub struct BesselResult {
    /// 계산된 함수 값들
    pub values: Vec<Complex64>,
    /// 언더플로우가 발생한 함수 값의 개수
    pub underflow_count: i32,
}
```

### `BesselError`
```rust
pub enum BesselError {
    /// 잘못된 입력 매개변수
    InvalidParameter(String),
    /// 계산 오류
    ComputationError(String),
}
```

## 원본 라이브러리 정보

이 Rust 바인딩은 다음을 기반으로 합니다:
- **원본 구현**: D.E. Amos의 Fortran 77 구현 (SLATEC 수학 라이브러리의 일부)
- **C++ 포트**: f2c 변환을 기반으로 한 현대적인 C++ 구현
- **특징**: 스레드 안전, f2c/gfortran 런타임 의존성 없음

## 라이센스

원본 zbessel 라이브러리의 라이센스를 따릅니다. 