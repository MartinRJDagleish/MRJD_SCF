cmake_minimum_required(VERSION 3.0.0)
project(MRJD_SCF)

include(CTest)
enable_testing()

# # stackoverlow 
# set(SOURCE_FILES main.cpp)
# link_directories(include)

# add_executable(molecule.o molecule.cc)
# add_executable(water.o water.cc)
# add_executable(water water.o molecule.o)
# add_library(molecule molecule.cc)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_link_libraries(MRJD_SCF molecule) 

add_executable(${PROJECT_NAME} main.cpp)
#* wihtout EIGEN
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_link_libraries(${PROJECT_NAME}  molecule)
add_subdirectory(include) 
#* with Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
# target_include_directories(${PROJECT_NAME} PUBLIC include include/eigen3/Eigen) 
# target_link_libraries(${PROJECT_NAME} PRIVATE molecule)
# add_subdirectory(include include/eigen3/Eigen) 

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
