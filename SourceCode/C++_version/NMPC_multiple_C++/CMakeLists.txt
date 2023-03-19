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

#添加头文件搜索路径
include_directories(/usr/local/include)

#添加库文件搜索路径
link_directories(/usr/local/lib)

#用于将当前目录下的所有源文件的名字保存在变量 DIR_SRCS 中
aux_source_directory(. DIR_SRCS)

#在这里根据名字boost_thread去寻找libboost_thread.a文件
target_link_libraries(CasADi_C_test boost_thread boost_system) #改为自己的项目生成的可执行文件

#add_executable(CasADi_C_test main.cpp)