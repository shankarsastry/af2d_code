add_executable(af2d_tests
    test_pslg.cpp
    test_distance.cpp
    test_lower_envelope.cpp
    test_ode_solver.cpp
    test_segment_splitter.cpp
    test_quadtree.cpp
    test_triangle_wrapper.cpp
    test_end_to_end.cpp
    test_io.cpp
)

target_link_libraries(af2d_tests PRIVATE af2d_lib GTest::gtest_main)
target_compile_definitions(af2d_tests PRIVATE CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

if(OpenMP_CXX_FOUND)
    target_link_libraries(af2d_tests PRIVATE OpenMP::OpenMP_CXX)
endif()

include(GoogleTest)
gtest_discover_tests(af2d_tests)
