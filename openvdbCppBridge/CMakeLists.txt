cmake_minimum_required(VERSION 3.18)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CXX_WARNINGS "-Wall -Wextra -Wpointer-arith -Wcast-align -fstrict-aliasing -Wno-unused-local-typedefs -Wno-misleading-indentation -Wno-maybe-uninitialized -Wno-int-in-bool-context -Wno-implicit-fallthrough -Wno-class-memaccess")

project(OpenVdbBridge)

add_executable(${PROJECT_NAME} src/main.cpp src/particleList.hpp src/rasterize.hpp src/utils.hpp src/pointGrid.hpp src/exportToObj.h)

# add openvdb library
list(APPEND CMAKE_MODULE_PATH $ENV{HOME}/openvdb/cmake/)
find_package(OpenVDB REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenVDB::openvdb)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_DL_LIBS})

# copy the executable to the source directory
add_custom_command(TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_CURRENT_SOURCE_DIR})
