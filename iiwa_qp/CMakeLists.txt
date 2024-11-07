cmake_minimum_required(VERSION 3.0.2)
project(iiwa_qp)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages

add_compile_options(-std=c++11)
find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_BUILD_TYPE Release) #From terminal: catkin config -DCMAKE_BUILD_TYPE=Release (or RelWithDebInfo)
include_directories(${catkin_INCLUDE_DIRS} )

set(OsqpEigen_DIR /catkin_ws/src/osqp-eigen/bin/lib/cmake/OsqpEigen)


find_package(catkin REQUIRED COMPONENTS
  message_generation
  nav_msgs
  geometry_msgs
  actionlib_msgs
  roscpp
  std_msgs
  tf
  tf_conversions
  eigen_conversions
  kdl_parser
  OsqpEigen
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")



catkin_package(
    INCLUDE_DIRS include ${catkin_INCLUDE_DIRS} 
    LIBRARIES separator)
    
include_directories( include ${catkin_INCLUDE_DIRS} )

#-----------------------------------------------------------------------------------------------------------------
###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
 
)


add_executable(iiwa_control src/osqpcontrol.cpp src/lbr_iiwa_7_robot.cpp )
target_link_libraries(iiwa_control PUBLIC OsqpEigen::OsqpEigen ${catkin_LIBRARIES} )
#add_dependencies(iiwa_control iiwa_qp)

