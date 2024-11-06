# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Yux")
include_directories("${CMAKE_SOURCE_DIR}/utils")
include_directories("${CMAKE_SOURCE_DIR}/transciphering")

#add_executable(test_Yux2_8_C1 test_Yux2_8_C1.cpp)

#target_link_libraries(test_Yux2_8_C1 PUBLIC Transcipher Symmetric helib)

#add_executable(test_Yux_p_C16 test_Yux_p_C16.cpp)

#target_link_libraries(test_Yux_p_C16 PUBLIC Transcipher Symmetric helib)

add_executable(test_Yux2_8_C16_Clientoff test_Yux2_8_C16_Clientoff.cpp)

target_link_libraries(test_Yux2_8_C16_Clientoff PUBLIC Transcipher Symmetric Utils helib)

add_executable(test_Yux2_8_C16_Serveroff test_Yux2_8_C16_Serveroff.cpp)

target_link_libraries(test_Yux2_8_C16_Serveroff PUBLIC Transcipher Symmetric Utils helib)

add_executable(test_Yux2_8_C16_Clienton test_Yux2_8_C16_Clienton.cpp)

target_link_libraries(test_Yux2_8_C16_Clienton PUBLIC Transcipher Symmetric Utils helib)

add_executable(test_Yux2_8_C16_Serveron test_Yux2_8_C16_Serveron.cpp)

target_link_libraries(test_Yux2_8_C16_Serveron PUBLIC Transcipher Symmetric Utils helib)

add_custom_command(TARGET test_Yux2_8_C16_Clientoff POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/run_Yux2_8_C16.sh $<TARGET_FILE_DIR:test_Yux2_8_C16_Clientoff>)


add_executable(test_Yus_p_C32_Clientoff test_Yus_p_C32_Clientoff.cpp)

target_link_libraries(test_Yus_p_C32_Clientoff PUBLIC Transcipher Symmetric Utils helib)

add_executable(test_Yus_p_C32_Serveroff test_Yus_p_C32_Serveroff.cpp)

target_link_libraries(test_Yus_p_C32_Serveroff PUBLIC Transcipher Symmetric Utils helib)

add_custom_command(TARGET test_Yus_p_C32_Clientoff POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/tests/run_Yus_p_C32.sh $<TARGET_FILE_DIR:test_Yus_p_C32_Clientoff>)

    
add_executable(test_Yus_p_C32_ClientAndServer test_Yus_p_C32_ClientAndServer.cpp)

target_link_libraries(test_Yus_p_C32_ClientAndServer PUBLIC Transcipher Symmetric Utils helib)