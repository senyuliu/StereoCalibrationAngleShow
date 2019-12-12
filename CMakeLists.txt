cmake_minimum_required(VERSION 2.6)

PROJECT(StereoCalibrationAngle)

add_compile_options(-g)
add_compile_options(-std=c++11)

############Boost####################
find_package( Boost 1.58.0 COMPONENTS filesystem )

FIND_PACKAGE(OpenCV 3 REQUIRED )
INCLUDE_DIRECTORIES( ${OPENCV_INCLUDE_DIR} )

####################################
include_directories(./)
ADD_EXECUTABLE(StereoCalibrationAngle main.cpp CeresBa.cpp shellCommand.cpp)

####################ceres###############
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

###########pcl#######################
find_package(PCL REQUIRED)

################Eigen################
find_package(Eigen3)

include_directories(${PCL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${PCL_LIBRARIES_DIRS})
add_definitions(${PCL_DEFINITIONS})
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")


TARGET_LINK_LIBRARIES(StereoCalibrationAngle opencv_core opencv_highgui opencv_flann opencv_imgproc opencv_highgui opencv_ml opencv_video opencv_objdetect opencv_features2d opencv_calib3d opencv_legacy opencv_contrib -lboost_system -lboost_regex -lboost_filesystem ${CERES_LIBRARIES} ${PCL_LIBRARIES})
