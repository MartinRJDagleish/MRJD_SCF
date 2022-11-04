cmake_minimum_required(VERSION 3.0.0)
project(MRJD_SCF VERSION 0.1.0)

include(CTest)
enable_testing()

# # stackoverlow 
# set(SOURCE_FILES main.cpp)
# link_directories(include)

# add_executable(molecule.o molecule.cc)
# add_executable(water.o water.cc)
# add_executable(water water.o molecule.o)
add_executable(MRJD_SCF main.cpp)
add_library(molcule molecule.cc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(MRJD_SCF molecule) 

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
