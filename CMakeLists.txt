cmake_minimum_required(VERSION 3.16.3)

if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()

project(ANTsWasm)

# Prefer an installed ANTs (>= the release containing ANTsX/ANTs#1973,
# which adds find_package(ANTS) support). Fall back to FetchContent-cloning
# the ANTs source tree for users on older ANTs that does not yet provide
# ANTSConfig.cmake.
find_package(ANTS CONFIG QUIET)
if(ANTS_FOUND)
  message(STATUS "ANTsWasm: using installed ANTs ${ANTS_VERSION} from ${ANTS_DIR}")
  set(_antswasm_ants_include_dirs ${ANTS_INCLUDE_DIRS})
else()
  message(STATUS "ANTsWasm: ANTSConfig.cmake not found; "
                 "falling back to FetchContent of ANTs source tree.")
  include(FetchContent)
  FetchContent_Declare(
    ants_inner
    GIT_REPOSITORY https://github.com/ANTsX/ANTs
    GIT_TAG        0e0e488fe17aa96dc25a8e88d2db42e5bdfefce7  # 2024-09-27
  )
  FetchContent_GetProperties(ants_inner)
  if(NOT ants_inner_POPULATED)
    FetchContent_Populate(ants_inner)
  endif()
  set(_antswasm_ants_include_dirs
    ${ants_inner_SOURCE_DIR}/Utilities
    ${ants_inner_SOURCE_DIR}/ImageRegistration
    ${ants_inner_SOURCE_DIR}/Tensor
    ${ants_inner_SOURCE_DIR}/Temporary
    ${ants_inner_SOURCE_DIR}/Examples
    ${ants_inner_SOURCE_DIR}/Examples/include
  )
endif()

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
endif()

set(ANTsWasm_INCLUDE_DIRS
  ${ANTsWasm_SOURCE_DIR}/include
  ${_antswasm_ants_include_dirs}
)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK 5.4 REQUIRED)  # we need at least 71c5c83 from 2024-02-02 for debug mode
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()
