cmake_minimum_required(VERSION 3.20)
project(FSMELD VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# The inherited numerical kernels use variable-length stack arrays (GNU extension).
set(CMAKE_CXX_EXTENSIONS ON)
if(MSVC)
    message(FATAL_ERROR "Use GCC/MinGW or Clang: the numerical kernels require GNU array extensions.")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(FFTW3_ROOT "" CACHE PATH "FFTW3 installation prefix")
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h
    HINTS "${FFTW3_ROOT}" ENV FFTW3_ROOT PATH_SUFFIXES include)
find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3-3 fftw3-3
    HINTS "${FFTW3_ROOT}" ENV FFTW3_ROOT PATH_SUFFIXES lib lib64)
if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
    message(FATAL_ERROR "Install FFTW3 (double precision), or set FFTW3_INCLUDE_DIR and FFTW3_LIBRARY.")
endif()

add_executable(fsmeld
    src/main.cpp
    src/fsmeld.cpp
    src/find_t_lp.cpp
    src/lb.cpp
    src/lower_upper_lemire.cpp
    src/q_to_h.cpp
    src/sip.cpp
    src/lb_normalize.cpp
    src/mvmean.cpp
    src/ElasticMeasures.cpp
    src/glb.cpp
    src/normalization_mode.cpp
    src/search_result.cpp)
target_include_directories(fsmeld PRIVATE include "${FFTW3_INCLUDE_DIR}")
target_link_libraries(fsmeld PRIVATE "${FFTW3_LIBRARY}")
# Keep the MinGW compiler runtime portable; FFTW remains an external dependency.
if(MINGW)
    target_link_options(fsmeld PRIVATE -static-libgcc -static-libstdc++)
endif()
install(TARGETS fsmeld RUNTIME DESTINATION bin)
