# 包含路径
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/HElib_v2.2.2/include)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/ntl-11.5.1/include)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/gmp/include)

# 库路径
link_directories(${CMAKE_SOURCE_DIR}/thirdparty/HElib_v2.2.2/lib)
link_directories(${CMAKE_SOURCE_DIR}/thirdparty/ntl-11.5.1/lib)
link_directories(${CMAKE_SOURCE_DIR}/thirdparty/gmp/lib)

# 源文件
set(SOURCES
    main.cpp
    tool.cpp
    Yux2_8.cpp
    Client.cpp
    Server.cpp
    FHEtool.cpp
)

# 头文件
set(HEADERS
    tool.hpp
    Yux2_8.hpp
    Client.hpp
    Server.hpp
    FHEtool.hpp
    params.hpp
    random_bit.hpp
)

# 生成可执行文件
add_executable(Yux2_8 ${SOURCES} ${HEADERS})

# 链接库
target_link_libraries(Yux2_8 helib ntl gmp)