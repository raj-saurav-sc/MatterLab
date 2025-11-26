# Add test executable
add_executable(test_stability test_stability.cpp)
target_link_libraries(test_stability)
target_include_directories(test_stability PRIVATE ${CMAKE_SOURCE_DIR})

# Enable testing
enable_testing()
add_test(NAME StabilityTests COMMAND test_stability)
