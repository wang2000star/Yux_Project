cmake_minimum_required(VERSION 3.16)

# 设置编译器
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(Yux LANGUAGES CXX)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置构建类型为 Release，启用优化
set(CMAKE_BUILD_TYPE "Release")

# 启用编译优化选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mavx2")

# # 启用 AVX-512 支持（如果支持）
# if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512vl")
# endif()

# # 启用 OpenMP 支持
# find_package(OpenMP REQUIRED)
# if (OpenMP_CXX_FOUND)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

# 查找 SEAL 库
find_package(SEAL 4.1.2 REQUIRED)
if (SEAL_FOUND)
    message("SEAL ${SEAL_VERSION} found.")
else()
    message(FATAL_ERROR "SEAL not found")
endif()

# Add include directories for SEAL
include_directories(${SEAL_INCLUDE_DIRS})

# 添加包含目录
include_directories(${CMAKE_SOURCE_DIR}/keccak)
include_directories(${CMAKE_SOURCE_DIR}/keccak/avx2)
include_directories(${CMAKE_SOURCE_DIR}/keccak/opt64)

# 添加 keccak 子目录
add_subdirectory(${CMAKE_SOURCE_DIR}/keccak)

# add UTILS library
aux_source_directory("${CMAKE_SOURCE_DIR}/utils" UTILS_SRC_LIST)
add_library(Utils ${UTILS_SRC_LIST})
# Link Keccak library to Utils
target_link_libraries(Utils PRIVATE keccak)
# Link SEAL library to Utils
target_link_libraries(Utils PRIVATE SEAL::seal)
# Output UTILS_SRC_LIST for debugging
message(STATUS "UTILS_SRC_LIST: ${UTILS_SRC_LIST}")

# add symmetric library
aux_source_directory("${CMAKE_SOURCE_DIR}/Symmetric" SYMMETRIC_SRC_LIST)
add_library(Symmetric ${SYMMETRIC_SRC_LIST})
# Link Utils library to Symmetric
target_link_libraries(Symmetric PRIVATE Utils)
# Link Keccak library to Symmetric
target_link_libraries(Symmetric PRIVATE keccak)
# Link SEAL library to Symmetric
target_link_libraries(Symmetric PRIVATE SEAL::seal)
# Output SYMMETRIC_SRC_LIST for debugging
message(STATUS "SYMMETRIC_SRC_LIST: ${SYMMETRIC_SRC_LIST}")

# # 针对特定文件启用 AVX2 优化选项
# set_source_files_properties(
#     ${CMAKE_SOURCE_DIR}/tests/test_Yus_p_C36_ClientAndServer5_2.cpp
#     PROPERTIES COMPILE_FLAGS "-mavx2"
# )

add_subdirectory(tests)

# 输出编译器标志以进行调试
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")