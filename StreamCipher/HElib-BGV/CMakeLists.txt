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

# 添加编译优化选项
add_compile_options(
    -O3                # 启用最大优化
    -march=native      # 针对本机 CPU 架构优化
)

# # 启用 AVX2 优化选项
# if(CMAKE_COMPILER_IS_GNUCXX)
#     add_compile_options(
#         -mavx2           # 启用 AVX2 支持
#     )
# endif()

# # 启用 AVX-512 支持（如果支持）
# if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
#     add_compile_options(
#         -mavx512vl        # 启用 AVX-512 整数向量运算支持
#     )
# endif()


# 启用 OpenMP 支持
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# 查找 HELIB 库
find_package(helib 2.3.0)#or 2.2.0
if (helib_NOTFOUND)
    message("helib not found")
else()
    message("helib ${HELIB_VERSION} found.")
endif()

# add symmetric library
aux_source_directory("${CMAKE_SOURCE_DIR}/Symmetric" SYMMETRIC_SRC_LIST)
add_library(Symmetric ${SYMMETRIC_SRC_LIST})
# Output SYMMETRIC_SRC_LIST for debugging
message(STATUS "SYMMETRIC_SRC_LIST: ${SYMMETRIC_SRC_LIST}")

# add UTILS library
aux_source_directory("${CMAKE_SOURCE_DIR}/utils" UTILS_SRC_LIST)
add_library(Utils ${UTILS_SRC_LIST})
# Output UTILS_SRC_LIST for debugging
message(STATUS "UTILS_SRC_LIST: ${UTILS_SRC_LIST}")

# # add TRANSCIPHER library
# aux_source_directory("${CMAKE_SOURCE_DIR}/transciphering" TRANSCIPHER_SRC_LIST)
# add_library(Transcipher ${TRANSCIPHER_SRC_LIST})
# # add NTL and GMP include directories
# target_include_directories(Transcipher PUBLIC "${CMAKE_SOURCE_DIR}/Symmetric" "${CMAKE_SOURCE_DIR}/utils" "/usr/local/include")
# # Output TRANSCIPHER_SRC_LIST for debugging
# message(STATUS "TRANSCIPHER_SRC_LIST: ${TRANSCIPHER_SRC_LIST}")

add_subdirectory(tests)