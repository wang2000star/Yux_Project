# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Yux")
include_directories("${CMAKE_SOURCE_DIR}/transciphering")

add_executable(test_Yux2_8_C1 test_Yux2_8_C1.cpp)

target_link_libraries(test_Yux2_8_C1 PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux_p_C16 test_Yux_p_C16.cpp)

target_link_libraries(test_Yux_p_C16 PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16_Clientoff test_Yux2_8_C16_Clientoff.cpp)

target_link_libraries(test_Yux2_8_C16_Clientoff PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16_Serveroff test_Yux2_8_C16_Serveroff.cpp)

target_link_libraries(test_Yux2_8_C16_Serveroff PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16_Clienton test_Yux2_8_C16_Clienton.cpp)

target_link_libraries(test_Yux2_8_C16_Clienton PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16_Serveron test_Yux2_8_C16_Serveron.cpp)

target_link_libraries(test_Yux2_8_C16_Serveron PUBLIC Transcipher Symmetric helib)

add_custom_command(TARGET test_Yux2_8_C1 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/run_Yux2_8_C16.sh $<TARGET_FILE_DIR:test_Yux2_8_C1>)