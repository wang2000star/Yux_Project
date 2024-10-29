# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Yux")
include_directories("${CMAKE_SOURCE_DIR}/transciphering")

add_executable(test_Yux2_8_C1 test_Yux2_8_C1.cpp)

target_link_libraries(test_Yux2_8_C1 PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16 test_Yux2_8_C16.cpp)

target_link_libraries(test_Yux2_8_C16 PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux_p_C16 test_Yux_p_C16.cpp)

target_link_libraries(test_Yux_p_C16 PUBLIC Transcipher Symmetric helib)