cmake_minimum_required (VERSION 3.16)
project(CNNmodel)

MATH(EXPR stack_size "256 * 1024 * 1024")
SET(CMAKE_C_COMPILER g++)
add_compile_options(-std=c++11)
set(SRC_LIST CNNmodel.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#需要的openCV库
FIND_PACKAGE(OpenCV REQUIRED)
#显示OpenCV的版本号
MESSAGE("OpenCV version : ${OpenCV_VERSION}")
MESSAGE("OpenCV include : ${OpenCV_INCLUDE_DIRS}")
MESSAGE("OpenCV libs11 : ${OpenCV_LIBS}")
MESSAGE("OpenCV DIR : ${OpenCV_DIR}")
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBS})
MESSAGE("OpenCV libs : ${OpenCV_LIBS}")

# openMP 配置
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

add_compile_options(-mavx2)
add_compile_options(-fopenmp)
set(SOURCE matclass.hpp CNNmodel.cpp matclass.cpp)
# set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,${stack_size}")
add_executable (${PROJECT_NAME} ${SOURCE})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${OpenCV_LIBS})