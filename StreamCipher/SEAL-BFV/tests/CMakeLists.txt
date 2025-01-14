# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Symmetric")
include_directories("${CMAKE_SOURCE_DIR}/utils")
include_directories("${CMAKE_SOURCE_DIR}/keccak")

add_executable(test_Yus_4 test_Yus_4.cpp)
target_link_libraries(test_Yus_4 PUBLIC Symmetric Utils SEAL::seal)

add_executable(test_Yus_5 test_Yus_5.cpp)
target_link_libraries(test_Yus_5 PUBLIC Symmetric Utils SEAL::seal)

add_executable(test_Yus_6 test_Yus_6.cpp)
target_link_libraries(test_Yus_6 PUBLIC Symmetric Utils SEAL::seal)

add_executable(test_HERA_4 test_HERA_4.cpp)
target_link_libraries(test_HERA_4 PUBLIC Symmetric Utils SEAL::seal)

add_executable(test_HERA_5 test_HERA_5.cpp)
target_link_libraries(test_HERA_5 PUBLIC Symmetric Utils SEAL::seal)


