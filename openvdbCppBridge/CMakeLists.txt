cmake_minimum_required(VERSION 3.18)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(OpenVdbBridge)

add_executable(${PROJECT_NAME}  src/main.cpp)

# add openvdb library
list(APPEND CMAKE_MODULE_PATH $ENV{HOME}/openvdb/cmake/)
find_package(OpenVDB REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenVDB::openvdb)

find_package(OpenMP)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_DL_LIBS})

# copy the executable to the source directory
add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_CURRENT_SOURCE_DIR})
