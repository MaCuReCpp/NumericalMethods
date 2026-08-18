cmake_minimum_required(VERSION 3.16)

project(NumericalMethodsCpp
    VERSION 0.1.0
    DESCRIPTION "A compact C++ numerical analysis library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(numerical_methods
    src/approximation.cpp
    src/differentiation_integration.cpp
    src/eigen_approx.cpp
    src/interpolation.cpp
    src/linear_algebra.cpp
    src/nonlinear.cpp
    src/numerical.cpp
    src/ode.cpp
    src/pde.cpp
    src/quantum.cpp)

target_include_directories(numerical_methods
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>)

target_compile_features(numerical_methods PUBLIC cxx_std_17)

add_library(NumericalMethods::numerical_methods ALIAS numerical_methods)

add_executable(numerical_demo examples/demo.cpp)
target_link_libraries(numerical_demo PRIVATE numerical_methods)

add_executable(numerical_cli examples/cli.cpp)
target_link_libraries(numerical_cli PRIVATE numerical_methods)

enable_testing()
add_executable(numerical_tests tests/test_numerical.cpp)
target_link_libraries(numerical_tests PRIVATE numerical_methods)
add_test(NAME numerical_tests COMMAND numerical_tests)
