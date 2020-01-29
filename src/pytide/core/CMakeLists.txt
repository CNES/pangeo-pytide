# Copyright (c) 2020 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
file(GLOB_RECURSE WRAPPED_SOURCES "*.cpp")
pybind11_add_module(core ${WRAPPED_SOURCES})

# BLAS
find_package(BLAS)

# MKL
if (DEFINED ENV{MKLROOT})
    find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS $ENV{MKLROOT}/include)
    if (MKL_INCLUDE_DIR)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        include_directories(${MKL_INCLUDE_DIR})
        target_link_libraries(core PRIVATE "${BLAS_LIBRARIES}")
    endif()
endif()
