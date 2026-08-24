# Numerical Methods C++ Library

작은 C++17 수치해석 라이브러리입니다. 학부 수치해석에서 자주 만나는 개념을 코드로 바로 실험할 수 있도록 구성했습니다.

## 구성

- `include/numerical/numerical.hpp`: `src`의 기능별 헤더를 모아주는 얇은 umbrella header
- `src/types.hpp`: 공통 타입, 결과 구조체, enum
- `src/errors.hpp`: 오차와 norm
- `src/roots.hpp`: 1변수 방정식의 근 찾기
- `src/differentiation_integration.hpp`: 수치미분/수치적분
- `src/special_functions.hpp`: 베셀 함수
- `src/interpolation.hpp`: 보간/스플라인
- `src/linear_algebra.hpp`: 선형시스템
- `src/eigen_approx.hpp`: 고유값/QR/SVD
- `src/approximation.hpp`: 근사이론/FFT
- `src/nonlinear.hpp`: 비선형 시스템
- `src/ode.hpp`: ODE 초기값/경계값 문제
- `src/pde.hpp`: PDE 유한차분
- `src/quantum.hpp`: 양자역학 학습용 선형대수/연산자/모형계 연산
- `src/errors.cpp`: 오차와 norm 구현
- `src/roots.cpp`: bisection, Newton, secant 구현
- `src/nonlinear.cpp`: 고정점, Steffensen, Muller, 비선형 시스템
- `src/differentiation_integration.cpp`: 수치미분/수치적분 구현
- `src/special_functions.cpp`: 베셀 함수 구현
- `src/interpolation.cpp`: 보간/스플라인 구현
- `src/linear_algebra.cpp`: 선형시스템 구현
- `src/eigen_approx.cpp`: power method, QR, SVD
- `src/approximation.cpp`: least squares, Chebyshev, FFT
- `src/ode.cpp`: ODE 초기값/경계값 문제 구현
- `src/pde.cpp`: BVP finite difference, heat/wave/Poisson solvers
- `src/quantum.cpp`: quantum states, operators, spin, wells, oscillators, finite-difference Hamiltonians
- `src/internal.hpp`: source 내부 공통 헬퍼
- `examples/demo.cpp`: 사용 예시
- `examples/cli.cpp`: 터미널 계산기형 실행 파일
- `tests/test_numerical.cpp`: 기본 검증 테스트
- `docs/book_implementation_map.md`: 교재 목차 기준 구현 범위 추적
- `CMakeLists.txt`: 라이브러리, 예제, 테스트 빌드 설정

## 기능별 직접 컴파일 파일

CMake를 쓰면 아래 파일 관계를 직접 신경 쓰지 않아도 됩니다. `cmake --build build`가 `CMakeLists.txt`에 등록된 모든 구현 파일을 정적 라이브러리 `libnumerical_methods.a`로 묶어줍니다.

하지만 `clang++`나 `g++`로 직접 컴파일할 때는, 사용한 함수가 정의된 `.cpp` 파일을 함께 넘겨야 합니다. 헤더만 include하면 선언만 보이는 것이고, 실제 함수 본문은 대응 `.cpp`에 있습니다.

기본 옵션:

```bash
clang++ -std=c++17 -Iinclude -Isrc <필요한 .cpp 파일들> your_program.cpp -o your_program
```

전체 라이브러리를 한 번에 포함하려면:

```bash
clang++ -std=c++17 -Iinclude -Isrc src/*.cpp your_program.cpp -o your_program
```

기능별 필요한 파일은 아래와 같습니다.

| 쓰는 기능 | include할 헤더 | 같이 컴파일할 `.cpp` |
| --- | --- | --- |
| 공통 타입만 사용 | `src/types.hpp` | 없음 |
| 오차, `l2_norm` | `src/errors.hpp` | `src/errors.cpp` |
| 이분법, Newton, secant | `src/roots.hpp` | `src/roots.cpp` |
| fixed-point, Steffensen, Muller, nonlinear Newton | `src/nonlinear.hpp` | `src/nonlinear.cpp`, `src/roots.cpp`, `src/linear_algebra.cpp`, `src/errors.cpp` |
| 유한차분 미분, 사다리꼴, Simpson, Richardson, Romberg, adaptive Simpson, Gaussian quadrature | `src/differentiation_integration.hpp` | `src/differentiation_integration.cpp` |
| 베셀 함수 | `src/special_functions.hpp` | `src/special_functions.cpp` |
| Lagrange, Neville, divided differences, Hermite, cubic spline | `src/interpolation.hpp` | `src/interpolation.cpp`, `src/errors.cpp` |
| Gaussian elimination, LU, determinant, inverse, Jacobi, Gauss-Seidel, SOR, conjugate gradient, iterative refinement | `src/linear_algebra.hpp` | `src/linear_algebra.cpp`, `src/errors.cpp` |
| power method, QR, QR eigenvalues, SVD | `src/eigen_approx.hpp` | `src/eigen_approx.cpp`, `src/linear_algebra.cpp`, `src/errors.cpp` |
| least squares, polynomial evaluation, Chebyshev, FFT | `src/approximation.hpp` | `src/approximation.cpp`, `src/linear_algebra.cpp`, `src/errors.cpp` |
| Euler, RK4, RKF45, Adams-Bashforth-Moulton, shooting BVP | `src/ode.hpp` | `src/ode.cpp` |
| finite-difference BVP, heat/wave equation, Poisson solver | `src/pde.hpp` | `src/pde.cpp`, `src/linear_algebra.cpp`, `src/errors.cpp` |
| 양자역학 상태, 연산자, spin, well, oscillator, Hamiltonian | `src/quantum.hpp` | `src/quantum.cpp`, `src/differentiation_integration.cpp`, `src/linear_algebra.cpp`, `src/errors.cpp` |
| umbrella header로 전체 API 사용 | `include/numerical/numerical.hpp` | 사용하는 기능에 맞는 위 `.cpp`들 또는 간단히 `src/*.cpp` |

예를 들어 Newton 방법만 쓰는 `my_root.cpp`를 컴파일하려면:

```bash
clang++ -std=c++17 -Iinclude -Isrc src/roots.cpp my_root.cpp -o my_root
```

선형시스템 풀이를 쓰는 `my_linear.cpp`라면:

```bash
clang++ -std=c++17 -Iinclude -Isrc src/errors.cpp src/linear_algebra.cpp my_linear.cpp -o my_linear
```

양자역학 spin 연산자를 쓰는 `my_quantum.cpp`라면:

```bash
clang++ -std=c++17 -Iinclude -Isrc src/errors.cpp src/linear_algebra.cpp src/differentiation_integration.cpp src/quantum.cpp my_quantum.cpp -o my_quantum
```

## 포함된 개념

### 오차와 수렴

수치해석은 정확한 해 대신 근사값을 다룹니다. 그래서 절대오차, 상대오차, 잔차, 반복 횟수, 수렴 여부를 함께 보는 것이 중요합니다.

- 절대오차: `|exact - approximate|`
- 상대오차: `|exact - approximate| / |exact|`
- 잔차: 방정식이나 선형시스템에 근사해를 대입했을 때 남는 값

지원 API:

- `absolute_error`
- `relative_error`
- `l2_norm`
- `Tolerance`
- `IterationResult`
- `LinearSolveResult`

### 비선형 방정식의 근

`f(x) = 0`을 만족하는 값을 찾는 방법입니다.

- 이분법: 구간 양 끝의 부호가 달라야 하며, 안정적입니다.
- 고정점 반복: `x = g(x)` 형태로 바꾸어 반복합니다.
- Newton 방법: 도함수를 사용해 빠르게 수렴할 수 있지만 초기값과 도함수에 민감합니다.
- 할선법: 도함수 없이 Newton 방법과 비슷한 형태로 근을 찾습니다.
- Steffensen 방법: 고정점 반복을 Aitken 가속으로 빠르게 만듭니다.
- Muller 방법: 2차 보간을 이용해 복소근까지 찾을 수 있습니다.

지원 API:

- `bisection`
- `fixed_point_iteration`
- `newton`
- `secant`
- `steffensen`
- `muller`

### 수치미분

도함수를 직접 구하기 어렵거나 함수값만 알 때 차분으로 미분값을 근사합니다.

- 전진차분
- 후진차분
- 중앙차분
- Richardson 외삽

지원 API:

- `finite_difference`
- `richardson_derivative`
- `DifferenceMethod`

### 수치적분

정적분을 닫힌형으로 계산하기 어렵거나 데이터 기반 함수일 때 적분값을 근사합니다.

- 사다리꼴 공식
- Simpson 공식
- Romberg 적분
- 적응 Simpson 적분
- Gauss-Legendre 구적법

지원 API:

- `trapezoidal_rule`
- `simpson_rule`
- `romberg_integration`
- `adaptive_simpson`
- `gaussian_quadrature`

### 베셀 함수

원통 대칭 문제, 진동, 파동, 열전도 문제에서 자주 등장하는 제1종 베셀 함수 `J_n(x)`를 계산합니다.

- `J_0(x)`, `J_1(x)`: 멱급수 전개로 계산합니다.
- `J_n(x)`: `J_0`, `J_1`에서 시작해 점화식으로 계산합니다.
- 음수 차수는 `J_{-n}(x) = (-1)^n J_n(x)` 관계를 사용합니다.

지원 API:

- `bessel_j0`
- `bessel_j1`
- `bessel_j`

### 보간

주어진 데이터 점들을 지나는 다항식을 만들어 중간 값을 예측합니다.

- Lagrange 보간
- Neville 보간
- Newton 분할차분 보간
- Hermite 보간
- 자연/고정단 3차 스플라인

지원 API:

- `lagrange_interpolate`
- `neville_interpolate`
- `divided_difference_coefficients`
- `evaluate_newton_polynomial`
- `hermite_interpolate`
- `natural_cubic_spline`
- `clamped_cubic_spline`
- `evaluate_spline`

### 선형시스템

`Ax = b` 형태의 연립일차방정식을 풉니다.

- Gaussian elimination: 직접법입니다. 부분 피벗팅을 사용합니다.
- LU factorization: 행렬을 하삼각/상삼각 행렬로 분해합니다.
- determinant/inverse: LU 분해 기반으로 계산합니다.
- Jacobi method: 반복법입니다. 대각우세 행렬에서 특히 잘 작동합니다.
- Gauss-Seidel/SOR: Jacobi보다 빠르게 수렴할 수 있는 반복법입니다.
- Conjugate gradient: 대칭 양의정부호 시스템을 효율적으로 풉니다.
- Iterative refinement: 이미 구한 해를 잔차 기반으로 개선합니다.
- power method: 지배 고유값과 고유벡터를 근사합니다.
- QR 분해, QR 고유값 반복, 교육용 SVD도 제공합니다.

지원 API:

- `gaussian_elimination`
- `lu_decomposition`
- `solve_lu`
- `determinant`
- `inverse`
- `jacobi`
- `gauss_seidel`
- `sor`
- `conjugate_gradient`
- `iterative_refinement`
- `power_method`
- `qr_decomposition`
- `qr_eigenvalues`
- `svd`

### 근사이론과 FFT

데이터를 다항식으로 근사하거나, Chebyshev 다항식과 FFT로 근사/변환 문제를 다룹니다.

지원 API:

- `least_squares_polynomial`
- `evaluate_polynomial`
- `chebyshev_nodes`
- `chebyshev_polynomial`
- `fft`

### 비선형 시스템

여러 변수의 연립 비선형 방정식 `F(x) = 0`을 Newton 방법으로 풉니다. Jacobian은 사용자가 제공합니다.

지원 API:

- `nonlinear_newton`

### 상미분방정식 초기값 문제

`y' = f(t, y)`, `y(t0) = y0` 형태의 문제를 시간 격자 위에서 근사합니다.

- Euler 방법
- 4차 Runge-Kutta 방법
- Runge-Kutta-Fehlberg 4(5) 적응 스텝 방법
- Adams-Bashforth-Moulton 4단계 predictor-corrector 방법

지원 API:

- `euler_method`
- `runge_kutta4`
- `runge_kutta_fehlberg45`
- `adams_bashforth_moulton4`

### 경계값 문제와 PDE

상미분방정식 경계값 문제와 열/파동/Poisson 방정식의 기본 유한차분 스킴을 제공합니다.

지원 API:

- `linear_shooting_bvp`
- `finite_difference_bvp`
- `heat_equation_explicit`
- `wave_equation_explicit`
- `poisson_dirichlet`

## 빌드와 테스트

CMake가 설치되어 있다면 아래 방식이 가장 편합니다.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

예제 실행:

```bash
./build/numerical_demo
```

계산기형 CLI 실행:

```bash
./build/numerical_cli help
./build/numerical_cli bessel j 2 1
./build/numerical_cli diff sin 0
./build/numerical_cli integrate simpson sin 0 3.141592653589793 100
./build/numerical_cli root newton sqrt2 1
./build/numerical_cli solve2 3 2 1 2 5 5
./build/numerical_cli quantum spin-z 1 0 0 0
./build/numerical_cli quantum well-energy 1 1 1
```

CMake가 없다면 `clang++` 또는 `g++`로 직접 컴파일할 수 있습니다.

```bash
clang++ -std=c++17 -Iinclude -Isrc src/*.cpp examples/demo.cpp -o numerical_demo
./numerical_demo
```

테스트 직접 컴파일:

```bash
clang++ -std=c++17 -Iinclude -Isrc src/*.cpp tests/test_numerical.cpp -o numerical_tests
./numerical_tests
```

`g++`를 사용한다면 위 명령에서 `clang++`만 `g++`로 바꾸면 됩니다.

파일이 여러 개로 나뉘었고 헤더도 `src`에 있으므로 직접 컴파일할 때는 모든 `src/*.cpp` 파일과 `-Isrc`를 함께 넘겨야 합니다.

```bash
clang++ -std=c++17 -Iinclude -Isrc src/*.cpp examples/cli.cpp -o numerical_cli
```

## 사용 예시

```cpp
#include "numerical/numerical.hpp"

#include <cmath>
#include <iostream>

int main()
{
    auto root = numerical::newton(
        [](double x) { return x * x - 2.0; },
        [](double x) { return 2.0 * x; },
        1.0);

    std::cout << root.value << "\n";

    double integral = numerical::simpson_rule(
        [](double x) { return std::sin(x); },
        0.0,
        3.14159265358979323846,
        100);

    std::cout << integral << "\n";
}
```
