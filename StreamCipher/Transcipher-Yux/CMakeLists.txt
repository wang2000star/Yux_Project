cmake_minimum_required(VERSION 3.16)

project(Yux LANGUAGES CXX)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## if Debug 
set(CMAKE_BUILD_TYPE "Debug")

#set(CMAKE_CXX_)
#add_compile_options("-ferror-limit=40")

find_package(helib 2.3.0)
if (helib_NOTFOUND)
    message("helib have not found")
else()
    message("helib ${HELIB_VERSION} have found.)")
endif()

# add symmetric library
aux_source_directory("${CMAKE_SOURCE_DIR}/Yux" SYMMETRIC_SRC_LIST)
add_library(Symmetric ${SYMMETRIC_SRC_LIST})
# add NTL and GMP include directories

#target_include_directories(Symmetric PUBLIC "${CMAKE_SOURCE_DIR}/../thirdparty/ntl-11.5.1/include" "${CMAKE_SOURCE_DIR}/../thirdparty/gmp-6.2.1/include")

# Output SYMMETRIC_SRC_LIST for debugging
message(STATUS "SYMMETRIC_SRC_LIST: ${SYMMETRIC_SRC_LIST}")

aux_source_directory("${CMAKE_SOURCE_DIR}/transciphering" TRANSCIPHER_SRC_LIST)
add_library(Transcipher ${TRANSCIPHER_SRC_LIST})
# add NTL and GMP include directories
target_include_directories(Transcipher PUBLIC "${CMAKE_SOURCE_DIR}/Yux" "/usr/local/include")
#target_include_directories(Transcipher PUBLIC "${CMAKE_SOURCE_DIR}/../thirdparty/ntl-11.5.1/include" "${CMAKE_SOURCE_DIR}/../thirdparty/gmp-6.2.1/include" "${CMAKE_SOURCE_DIR}/Yux")

# Output TRANSCIPHER_SRC_LIST for debugging
message(STATUS "TRANSCIPHER_SRC_LIST: ${TRANSCIPHER_SRC_LIST}")

add_subdirectory(tests)