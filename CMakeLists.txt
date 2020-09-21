cmake_minimum_required(VERSION 3.10)
project(tftest)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++14")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
#-DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")
set(LIB_PREFIX g2o_)
SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )
# 의존성 패키지.
find_package(Eigen3)
FIND_PACKAGE (PCL 1.2 REQUIRED)  # 프로그램 실행시 필요한 패키지
# 없을 경우 에러 발생

# 헤더 디렉토리 지정 (-I)
INCLUDE_DIRECTORIES(
        ${PCL_INCLUDE_DIRS}
)

# 라이브러리 디렉토리 지정 (-L)
LINK_DIRECTORIES(${PCL_LIBRARY_DIRS})
set(mylibSRCS main.cpp geofunc.cpp)
add_library(nav SHARED ${mylibSRCS} ${PCL_INCLUDE_DIRS})
#add_library(pcllib_apps SHARED IMPORTED)
#add_library(pcllib_search SHARED IMPORTED)

# 전처리기 매크로 추가 (-D)
ADD_DEFINITIONS(${COMPILE_FLAGS})
ADD_DEFINITIONS(${PCL_DEFINITIONS})
add_executable(tftest main.cpp geofunc.cpp geofunc.h)
TARGET_LINK_LIBRARIES(tftest ${PCL_LIBRARIES})    #( <Target_이름> <라이브러리> <라이브러리> ... )
# 실행 파일생성하기에 앞서 링크 해야 하는 라이브러리와 실행 파일 링크