# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Yux")
include_directories("${CMAKE_SOURCE_DIR}/transciphering")

add_executable(test_Yux2_8 test_Yux2_8.cpp)

target_link_libraries(test_Yux2_8 PUBLIC Transcipher Symmetric helib)