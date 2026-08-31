cmake_minimum_required(VERSION 3.16)
project(FSMELD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(fsmeld
    main.cpp
    FSMELD.cpp
    find_t_lp.cpp
    lb.cpp
    lower_upper_lemire.cpp
    q_to_h.cpp
    sip.cpp
    lb_normalize.cpp
    mvmean.cpp
    ElasticMeasures.cpp
    glb.cpp
)

target_include_directories(fsmeld PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(fsmeld PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/libfftw3-3.lib")

if(MSVC)
    target_compile_options(fsmeld PRIVATE /O2)
else()
    target_compile_options(fsmeld PRIVATE -O2)
endif()

if(MINGW)
    target_link_options(fsmeld PRIVATE -static)
endif()

if(WIN32 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/libfftw3-3.dll")
    add_custom_command(TARGET fsmeld POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/lib/libfftw3-3.dll"
            "$<TARGET_FILE_DIR:fsmeld>/libfftw3-3.dll"
    )
endif()
