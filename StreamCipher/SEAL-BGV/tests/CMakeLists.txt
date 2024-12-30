# Add include directories for tests
include_directories(${CMAKE_SOURCE_DIR})
include_directories("${CMAKE_SOURCE_DIR}/Symmetric")
include_directories("${CMAKE_SOURCE_DIR}/utils")


add_executable(test_Yus test_Yus.cpp)
target_link_libraries(test_Yus PUBLIC Symmetric Utils SEAL::seal ${OpenMP_CXX_LIBRARIES})


