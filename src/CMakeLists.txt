add_library(af2d_lib STATIC
    pslg.cpp
    distance.cpp
    lower_envelope.cpp
    ode_solver.cpp
    segment_splitter.cpp
    triangle_wrapper.cpp
    io.cpp
)

target_include_directories(af2d_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(af2d_lib PUBLIC triangle)

if(OpenMP_CXX_FOUND)
    target_link_libraries(af2d_lib PRIVATE OpenMP::OpenMP_CXX)
endif()
