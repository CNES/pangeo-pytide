# Copyright (c) 2022 CNES
#
# All rights reserved. Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.
file(GLOB_RECURSE WRAPPED_SOURCES "*.cpp")
pybind11_add_module(core ${WRAPPED_SOURCES})

# BLAS
if(DEFINED ENV{MKLROOT})
  set(BLA_VENDOR Intel10_64_dyn)
endif()
find_package(BLAS)
if(NOT BLAS_FOUND AND DEFINED ENV{MKLROOT})
  set(BLA_VENDOR Intel10_64ilp)
  find_package(BLAS)
endif()

# MKL
if(DEFINED ENV{MKLROOT})
  find_path(
    MKL_INCLUDE_DIR
    NAMES mkl.h
    HINTS $ENV{MKLROOT}/include)
  if(MKL_INCLUDE_DIR)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    add_definitions(-DMKL_ILP64)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(core PRIVATE "${BLAS_LIBRARIES}")
  endif()
endif()
