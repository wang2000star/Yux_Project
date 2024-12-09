# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Symmetric")
include_directories("${CMAKE_SOURCE_DIR}/utils")
# include_directories("${CMAKE_SOURCE_DIR}/transciphering")

#add_executable(test_Yux2_8_C1 test_Yux2_8_C1.cpp)

#target_link_libraries(test_Yux2_8_C1 PUBLIC Transcipher Symmetric helib)

#add_executable(test_Yux_p_C16 test_Yux_p_C16.cpp)

#target_link_libraries(test_Yux_p_C16 PUBLIC Transcipher Symmetric helib)

#add_executable(test_Yux2_8_C16_Clientoff test_Yux2_8_C16_Clientoff.cpp)

#target_link_libraries(test_Yux2_8_C16_Clientoff PUBLIC Transcipher Symmetric Utils helib)

#add_executable(test_Yux2_8_C16_Serveroff test_Yux2_8_C16_Serveroff.cpp)

#target_link_libraries(test_Yux2_8_C16_Serveroff PUBLIC Transcipher Symmetric Utils helib)

#add_executable(test_Yux2_8_C16_Clienton test_Yux2_8_C16_Clienton.cpp)

#target_link_libraries(test_Yux2_8_C16_Clienton PUBLIC Transcipher Symmetric Utils helib)

#add_executable(test_Yux2_8_C16_Serveron test_Yux2_8_C16_Serveron.cpp)

#target_link_libraries(test_Yux2_8_C16_Serveron PUBLIC Transcipher Symmetric Utils helib)

#add_custom_command(TARGET test_Yux2_8_C16_Clientoff POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#    ${CMAKE_SOURCE_DIR}/tests/run_Yux2_8_C16.sh $<TARGET_FILE_DIR:test_Yux2_8_C16_Clientoff>)


#add_executable(test_Yus_p_C32_Clientoff test_Yus_p_C32_Clientoff.cpp)

#target_link_libraries(test_Yus_p_C32_Clientoff PUBLIC Transcipher Symmetric Utils helib)

#add_executable(test_Yus_p_C32_Serveroff test_Yus_p_C32_Serveroff.cpp)

#target_link_libraries(test_Yus_p_C32_Serveroff PUBLIC Transcipher Symmetric Utils helib)

#add_custom_command(TARGET test_Yus_p_C32_Clientoff POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#    ${CMAKE_SOURCE_DIR}/tests/run_Yus_p_C32.sh $<TARGET_FILE_DIR:test_Yus_p_C32_Clientoff>)

#add_executable(test_Yus_p_C64_ClientAndServer test_Yus_p_C64_ClientAndServer.cpp)
#target_link_libraries(test_Yus_p_C64_ClientAndServer PUBLIC Transcipher Symmetric Utils helib)

# add_executable(test_Yus_p_C32_ClientAndServer test_Yus_p_C32_ClientAndServer.cpp)
# target_link_libraries(test_Yus_p_C32_ClientAndServer PUBLIC  Symmetric Utils helib)

# add_executable(test_Yus_p_C32_ClientAndServer2 test_Yus_p_C32_ClientAndServer2.cpp)
# target_link_libraries(test_Yus_p_C32_ClientAndServer2 PUBLIC  Symmetric Utils helib)

#add_executable(test_Yux_p_C16_ClientAndServer test_Yux_p_C16_ClientAndServer.cpp)
#target_link_libraries(test_Yux_p_C16_ClientAndServer PUBLIC Transcipher Symmetric Utils helib)

#add_executable(test_Yux_p_C16_ClientAndServer2 test_Yux_p_C16_ClientAndServer2.cpp)
#target_link_libraries(test_Yux_p_C16_ClientAndServer2 PUBLIC Transcipher Symmetric Utils helib)

add_executable(test_Hera test_Hera.cpp)
target_link_libraries(test_Hera PUBLIC  Symmetric Utils helib)

add_executable(test_context test_context.cpp)
target_link_libraries(test_context PUBLIC  helib)

add_executable(test_Yus_p_C48_ClientAndServer3 test_Yus_p_C48_ClientAndServer3.cpp)
target_link_libraries(test_Yus_p_C48_ClientAndServer3 PUBLIC Symmetric Utils helib)

add_executable(test_Yus_p_C64_ClientAndServer4 test_Yus_p_C64_ClientAndServer4.cpp)
target_link_libraries(test_Yus_p_C64_ClientAndServer4 PUBLIC Symmetric Utils helib)

# add_custom_command(TARGET test_Yus_p_C64_ClientAndServer4 POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#    ${CMAKE_SOURCE_DIR}/tests/run_test_Yus_p_C64_ClientAndServer4.sh $<TARGET_FILE_DIR:test_Yus_p_C64_ClientAndServer4>)