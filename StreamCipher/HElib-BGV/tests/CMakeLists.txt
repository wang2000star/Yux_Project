# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Symmetric")
include_directories("${CMAKE_SOURCE_DIR}/utils")
include_directories("${CMAKE_SOURCE_DIR}/keccak")


# add_executable(test_context test_context.cpp)
# target_link_libraries(test_context PUBLIC  helib)

add_executable(test_HERA_4 test_HERA_4.cpp)
target_link_libraries(test_HERA_4 PUBLIC Symmetric Utils helib)

add_executable(test_HERA_5 test_HERA_5.cpp)
target_link_libraries(test_HERA_5 PUBLIC Symmetric Utils helib)

add_executable(test_pasta2_3 test_pasta2_3.cpp)
target_link_libraries(test_pasta2_3 PUBLIC Symmetric Utils keccak helib)

add_executable(test_pasta2_4 test_pasta2_4.cpp)
target_link_libraries(test_pasta2_4 PUBLIC Symmetric Utils keccak helib)

add_executable(test_pasta_3 test_pasta_3.cpp)
target_link_libraries(test_pasta_3 PUBLIC Symmetric Utils keccak helib)

add_executable(test_pasta_4 test_pasta_4.cpp)
target_link_libraries(test_pasta_4 PUBLIC Symmetric Utils keccak helib)

add_executable(test_Yus_4 test_Yus_4.cpp)
target_link_libraries(test_Yus_4 PUBLIC Symmetric Utils keccak helib)

add_executable(test_Yus_5 test_Yus_5.cpp)
target_link_libraries(test_Yus_5 PUBLIC Symmetric Utils keccak helib)

add_executable(test_Yus_6 test_Yus_6.cpp)
target_link_libraries(test_Yus_6 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_HERA_4_33 test_HERA_4_33.cpp)
# target_link_libraries(test_HERA_4_33 PUBLIC Symmetric Utils helib)

# add_executable(test_HERA_5_33 test_HERA_5_33.cpp)
# target_link_libraries(test_HERA_5_33 PUBLIC Symmetric Utils helib)

# add_executable(test_pasta2_3_33 test_pasta2_3_33.cpp)
# target_link_libraries(test_pasta2_3_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_pasta2_4_33 test_pasta2_4_33.cpp)
# target_link_libraries(test_pasta2_4_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_pasta_3_33 test_pasta_3_33.cpp)
# target_link_libraries(test_pasta_3_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_pasta_4_33 test_pasta_4_33.cpp)
# target_link_libraries(test_pasta_4_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_Yus_4_33 test_Yus_4_33.cpp)
# target_link_libraries(test_Yus_4_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_Yus_5_33 test_Yus_5_33.cpp)
# target_link_libraries(test_Yus_5_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_Yus_6_33 test_Yus_6_33.cpp)
# target_link_libraries(test_Yus_6_33 PUBLIC Symmetric Utils keccak helib)

# add_executable(test_Yux_14 test_Yux_14.cpp)
# target_link_libraries(test_Yux_14 PUBLIC Symmetric Utils helib)



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

# add_executable(test_Hera test_Hera.cpp)
# target_link_libraries(test_Hera PUBLIC  Symmetric Utils helib)



# add_executable(test_Yus_p_C48_ClientAndServer3 test_Yus_p_C48_ClientAndServer3.cpp)
# target_link_libraries(test_Yus_p_C48_ClientAndServer3 PUBLIC Symmetric Utils helib)

# add_executable(test_Yus_p_C64_ClientAndServer4 test_Yus_p_C64_ClientAndServer4.cpp)
# target_link_libraries(test_Yus_p_C64_ClientAndServer4 PUBLIC Symmetric Utils helib)

# add_executable(test_Yus_p_C36_ClientAndServer5_2 test_Yus_p_C36_ClientAndServer5_2.cpp)
# target_link_libraries(test_Yus_p_C36_ClientAndServer5_2 PUBLIC Symmetric Utils helib ${OpenMP_CXX_LIBRARIES})
# add_executable(test_Yus_p_C36_ClientAndServer7 test_Yus_p_C36_ClientAndServer7.cpp)
# target_link_libraries(test_Yus_p_C36_ClientAndServer7 PUBLIC Symmetric Utils helib)


# add_executable(test_Yus_p_C18_ClientAndServer8 test_Yus_p_C18_ClientAndServer8.cpp)
# target_link_libraries(test_Yus_p_C18_ClientAndServer8 PUBLIC Symmetric Utils helib)

# add_executable(test_Yus_p_C18_ClientAndServer8_2 test_Yus_p_C18_ClientAndServer8_2.cpp)
# target_link_libraries(test_Yus_p_C18_ClientAndServer8_2 PUBLIC Symmetric Utils helib ${OpenMP_CXX_LIBRARIES})
# # 定义参数数对
# set(param_pairs
#   "4 0"
#   "4 1"
#   "4 6"
#   "4 7"
#   "5 0"
#   "5 1"
#   "5 2"
#   "5 3"
#   "6 0"
#   "6 1"
#   "6 2"
#   "7 0"
#   "8 0"
# )

# # 创建子文件夹
# set(subdir "${CMAKE_SOURCE_DIR}/tests/generated")
# file(MAKE_DIRECTORY ${subdir})

# # 遍历所有参数数对
# foreach(pair ${param_pairs})
#   string(REPLACE " " ";" pair_list ${pair})
#   list(GET pair_list 0 Nr)
#   list(GET pair_list 1 idx)

#   # 设置源文件和可执行文件名称
#   set(temp_file "${subdir}/test_Yus_p_C64_ClientAndServer6_${Nr}_${idx}.cpp")
#   set(executable_name "test_Yus_p_C64_ClientAndServer6_${Nr}_${idx}")

#   # 读取原始源文件内容
#   file(READ "${CMAKE_SOURCE_DIR}/tests/test_Yus_p_C64_ClientAndServer6.cpp" file_content)

#   # 替换 Nr 和 idx 的值
#   string(REPLACE "constexpr unsigned Nr = 4" "constexpr unsigned Nr = ${Nr}" file_content "${file_content}")
#   string(REPLACE "constexpr long idx = 0" "constexpr long idx = ${idx}" file_content "${file_content}")

#   # 写入修改后的内容到新的文件
#   file(WRITE "${temp_file}" "${file_content}")

#   # 添加可执行文件和链接库
#   add_executable(${executable_name} ${temp_file})
#   target_link_libraries(${executable_name} PUBLIC Symmetric Utils helib)
# endforeach()

# add_custom_command(TARGET test_Yus_p_C64_ClientAndServer4 POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#    ${CMAKE_SOURCE_DIR}/tests/run_test_Yus_p_C64_ClientAndServer4.sh $<TARGET_FILE_DIR:test_Yus_p_C64_ClientAndServer4>)

# add_custom_command(TARGET test_Yus_p_C64_ClientAndServer4 POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy
#    ${CMAKE_SOURCE_DIR}/tests/run_test_Yus_p_C64_ClientAndServer4.sh $<TARGET_FILE_DIR:test_Yus_p_C64_ClientAndServer4>)