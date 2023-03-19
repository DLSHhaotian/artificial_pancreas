cmake_minimum_required(VERSION 3.10)
project(CasADi_C_test)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(CASADI REQUIRED)
include_directories(${CASADI_INCLUDE_DIR})
include_directories(/usr/include/eigen3)
add_executable(${PROJECT_NAME} main.cpp MVPmodel.cpp HovorkaModel.cpp CasADiSolver.cpp explicitRK4.cpp functionUsed.cpp)
target_link_libraries(${PROJECT_NAME} ${CASADI_LIBRARIES})

set(BOOST_ROOT "/usr/local/include/boost")

include_directories(/usr/local/include)

link_directories(/usr/local/lib)

aux_source_directory(. DIR_SRCS)

target_link_libraries(CasADi_C_test boost_thread boost_system) 
#add_executable(CasADi_C_test main.cpp)