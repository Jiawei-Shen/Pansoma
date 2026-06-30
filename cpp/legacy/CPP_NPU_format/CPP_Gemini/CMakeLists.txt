cmake_minimum_required(VERSION 3.15)
project(GamProcessor CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required libraries
find_package(Protobuf REQUIRED)
find_package(ZLIB REQUIRED)

# Add our source directory and include directory
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}
)

# ---- Protobuf Compilation ----
# Specify the location of the .proto file
set(PROTO_FILES ${CMAKE_CURRENT_SOURCE_DIR}/proto/vg.proto)

# Generate C++ sources from the .proto file
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${PROTO_FILES})

# Add the executable target
add_executable(gam_processor
    src/main.cpp
    src/pipeline.cpp
    ${PROTO_SRCS}
    ${PROTO_HDRS}
)

# --- CORRECTED SECTION ---
# Tell the compiler where to find the Protobuf headers
target_include_directories(gam_processor
    PRIVATE
    ${Protobuf_INCLUDE_DIRS}
)
# --- END CORRECTED SECTION ---


# Link libraries to the executable
target_link_libraries(gam_processor
    PRIVATE
    ${PROTOBUF_LIBRARIES}
    ${ZLIB_LIBRARIES}
)

# Add a hint for the compiler to optimize
target_compile_options(gam_processor PRIVATE -O3)