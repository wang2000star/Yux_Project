# 查找 OpenMP 包
find_package(OpenMP REQUIRED)

# 设置编译器选项以禁止所有警告
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# 添加可执行文件 homoAES，并链接 OpenMP 库
add_executable(homoAES AES.cpp AES_homo.cpp)
target_link_libraries(homoAES tfhe++ OpenMP::OpenMP_CXX)

# 其他注释掉的代码
# file(GLOB test_sources RELATIVE "${CMAKE_CURRENT_LIST_DIR}" "*.cpp")
# foreach(test_source ${test_sources})
#     string( REPLACE ".cpp" "" test_name ${test_source} )
#     add_executable(${test_name} ${test_source})
#     target_link_libraries(${test_name} tfhe++ )
# endforeach(test_source ${test_sources})

# # For GMP used in INTorus test
# find_library(gmp NAMES gmp gmpxx)
# if(gmp)
#     add_compile_definitions(ENABLE_LSHTEST)
#     target_link_libraries(INTorus gmp)
# endif()

# add_subdirectory(axell)