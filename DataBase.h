#ifndef __DATABASE_H__
#define __DATABASE_H__
#include <opencv2/opencv.hpp>
#include <fstream>
#include <iostream>
#include <string>


#define PI 3.1415926

using namespace std; 
using namespace cv; 

//#define IMU 10000

const double pose[6] = {544383.0779,3374202.4282,35.2317, 0, 0, 0};//{544383.0779,3374202.4282,35.2317, 0, 0, 0};//{499501.337341,3374115.9899,44.1299, 0, 0, 0};//{544388.9955,3374216.1275,32.033,0.33456,-0.06835,-146.6993};//{544353.5895,3374207.5119,32.130,0.0680,0.51399,137.60565};//{499564.147548,3374109.520344,35.459608, 0, 0,0};//{544389.0781,3374216.0524,32.025, 0.0158,0.30185,-0.07244};//{545679.738,3366650.734,13.329,0.875,-0.931,61.296}; //{545800.631,3366719.911,12.856,0.197,-0.854,60.819};{546100.7107,3370454.8788,19.503, -0.16620,1.13705,-93.02298}; 

namespace IBD
{
    struct DataUV
    {
        double ptx; 
        double pty; 
        double Pitch;    //pay attention to this points 
        double Roll; 
        double Yaw; 
        double XShift; 
        double YShift; 
        double ZShift; 
    }; 
    struct DataUVXYZ
    {
        double ptx; 
        double pty;
        double ptxr; 
        double ptyr;   
        double XGlobal; 
        double YGlobal; 
        double ZGlobal; 
        std::string name; 
        double XShift; 
        double YShift; 
        double ZShift;
        double Roll;
        double Pitch; 
        double Yaw; 
    };

    bool loadUVleft(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints);
    bool loadUVleft(std::string path, std::vector<IBD::DataUV>& imagePoints); 
    bool loadUVleft(std::string path, std::vector<cv::Point2d>& imagePoints); 
    bool loadXYZleft(std::string path, std::vector<cv::Point3d>& objectPoints); 
    bool loadXYZleft(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints, std::vector<cv::Point3d>& objectPoints);

    bool loadUVright(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints);
    bool loadUVright(std::string path, std::vector<IBD::DataUV>& imagePoints); 
    bool loadUVright(std::string path, std::vector<cv::Point2d>& imagePoints); 
    bool loadXYZright(std::string path, std::vector<cv::Point3d>& objectPoints);
    bool loadXYZright(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints, std::vector<cv::Point3d>& objectPoints);

    bool calcLeftRT(Mat cameraMatrix, Mat distCoeffs, string pathuv, string pathxyz,std::vector<cv::Point3d>& objPointsLeft, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageLeft);
    bool calcRightRT(Mat cameraMatrix, Mat distCoeffs, string pathuv, string pathxyz,std::vector<cv::Point3d>& objPointsRight, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageRight);
    bool calcLeftRT(Mat cameraMatrix, Mat distCoeffs, string path, std::vector<cv::Point3d>& objPointsLeft, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageLeft);
    bool calcRightRT(Mat cameraMatrix, Mat distCoeffs, string path, std::vector<cv::Point3d>& objPointsRight, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageRight);

};

bool IBD::calcLeftRT(Mat cameraMatrix, Mat distCoeffs, string pathuv, string pathxyz,std::vector<cv::Point3d>& objPointsLeft, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageLeft)
{
  // Read points
  std::vector<IBD::DataUV> imagePoints ;
  std::vector<cv::Point3d> objectPoints;

  std::string path_uv = pathuv; 
  std::string path_xyz= pathxyz; 

  //return 0;

  IBD::loadUVleft(path_uv,  imagePoints);
  IBD::loadXYZleft(path_xyz, objectPoints);  
 
  std::cout << "There are " << imagePoints.size() << " imagePoints and " << objectPoints.size() << " objectPoints." << std::endl;

  //objPoints New 
  std::vector<cv::Point3d> objPointsLocal ;
  std::vector<cv::Point2d> imgPointsLocal ;

  assert(imagePoints.size() == objectPoints.size());

  for(int i = 0; i < imagePoints.size() ; i++)
  {
      double uimage = imagePoints[i].ptx; 
      double vimage = imagePoints[i].pty; 

      imageLeft.push_back(cv::Point2d(uimage,vimage)); 

      double Xshift = imagePoints[i].XShift; 
      double Yshift = imagePoints[i].YShift;
      double Zshift = imagePoints[i].ZShift;
      double rolls   = -imagePoints[i].Roll;
      double pitchs  = -imagePoints[i].Pitch;
      double yaws    = imagePoints[i].Yaw;

      //RT
      Mat RX = (Mat_<double>(3,3)<<1, 0,                  0,
                                   0, cos(pitchs*PI/180), -sin(pitchs*PI/180),
                                   0, sin(pitchs*PI/180), cos(pitchs*PI/180) );
        
      Mat RY = (Mat_<double>(3,3)<<cos(rolls*PI/180),0, sin(rolls*PI/180),
                                   0,                1,   0,
                                   -sin(rolls*PI/180),0, cos(rolls*PI/180) );
      Mat RZ = (Mat_<double>(3,3)<<cos(yaws*PI/180), -sin(yaws*PI/180),  0,
                                   sin(yaws*PI/180), cos(yaws*PI/180), 0,
                                   0, 0,                 1 ); 
        
      Mat TT = (Mat_<double>(3,1)<<Xshift, Yshift, Zshift); 
      Mat cloudPt = (Mat_<double>(3,1) << objectPoints[i].x , objectPoints[i].y, objectPoints[i].z);
      cloudPt = cloudPt - TT; 

      transpose(RX, RX);
      transpose(RY, RY);
      transpose(RZ, RZ); 
      Mat cloudPtNew = RY*RX*RZ*cloudPt; 
   
      cv::Point3d pointslocal; 
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = cloudPtNew.at<double>(1, 0);
      pointslocal.z = cloudPtNew.at<double>(2, 0);

      std::cout<<fixed<<setprecision(10)<<pointslocal.x << "," <<pointslocal.y << "," <<pointslocal.z << std::endl; 

      cv::Point2d pointsuv(uimage, vimage);
 
      objPointsLeft.push_back(pointslocal); 
      objPointsLocal.push_back(pointslocal);    
      imgPointsLocal.push_back(pointsuv);   
  }

  bool useExtrinsicGuess = true; 
  int  iterationsCount   = 50000; 
  double reprojectionError= 2.0; 
  double  minInliersCount = 0.96; 
  cv::Mat inliers ; 
  int flags              = cv::SOLVEPNP_ITERATIVE;

  std::cout <<"cv::SOLVEPNP_ITERATIVE: " << cv::SOLVEPNP_ITERATIVE << std::endl; 
  flags = 0; 

  //cv::solvePnPRansac(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec, useExtrinsicGuess,iterationsCount,reprojectionError,minInliersCount,inliers,flags);

  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec, tvec);
 
  std::cout << "rvec: " << rvec << std::endl;
  std::cout << "tvec: " << tvec << std::endl;
 
  std::vector<cv::Point2d> projectedPoints; 
  cv::projectPoints(objPointsLocal, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);
 
  int rows = 2160; 
  int cols = 4096; 
  
  cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
 
  for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  {
    std::cout << "Image point: " << imagePoints[i].ptx<<"," << imagePoints[i].pty << " Projected to " << projectedPoints[i] << std::endl;
     circle(matShow, cv::Point(imagePoints[i].ptx, imagePoints[i].pty), 10,  Scalar(255,0,0), 2, 8); 
     circle(matShow, projectedPoints[i], 8,  Scalar(0,0,255), 2, 8); 
  }
  
  namedWindow("show", 2); 
  imshow("show", matShow);
  waitKey(0); 
  return true; 
}

bool IBD::calcRightRT(Mat cameraMatrix, Mat distCoeffs, string pathuv, string pathxyz,std::vector<cv::Point3d>& objPointsRight, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageRight)
{
  // Read points
  std::vector<IBD::DataUV> imagePoints ;
  std::vector<cv::Point3d> objectPoints;

  std::string path_uv = pathuv; 
  std::string path_xyz= pathxyz; 

  IBD::loadUVright(path_uv,  imagePoints);
  IBD::loadXYZright(path_xyz,objectPoints);  
 
  std::cout << "There are " << imagePoints.size() << " imagePoints and " << objectPoints.size() << " objectPoints." << std::endl;

  //objPoints New 
  std::vector<cv::Point3d> objPointsLocal ;
  std::vector<cv::Point2d> imgPointsLocal ;

  assert(imagePoints.size() == objectPoints.size());

  for(int i = 0; i < imagePoints.size() ; i++)
  {
      double uimage = imagePoints[i].ptx; 
      double vimage = imagePoints[i].pty; 
      imageRight.push_back(cv::Point2d(uimage,vimage)); 

      double Xshift = imagePoints[i].XShift; 
      double Yshift = imagePoints[i].YShift;
      double Zshift = imagePoints[i].ZShift;
      double rolls   = -imagePoints[i].Roll;
      double pitchs  = -imagePoints[i].Pitch;
      double yaws    = imagePoints[i].Yaw;

      //RT
      Mat RX = (Mat_<double>(3,3)<<1, 0,                  0,
                                   0, cos(pitchs*PI/180), sin(pitchs*PI/180),
                                   0, -sin(pitchs*PI/180), cos(pitchs*PI/180) );
        
      Mat RY = (Mat_<double>(3,3)<<cos(rolls*PI/180),0, -sin(rolls*PI/180),
                                   0,                1,   0,
                                   sin(rolls*PI/180),0, cos(rolls*PI/180) );
      Mat RZ = (Mat_<double>(3,3)<<cos(yaws*PI/180), sin(yaws*PI/180),  0,
                                   -sin(yaws*PI/180), cos(yaws*PI/180), 0,
                                   0, 0,                 1 ); 
        
      Mat TT = (Mat_<double>(3,1)<<Xshift, Yshift, Zshift); 
      Mat cloudPt = (Mat_<double>(3,1) << objectPoints[i].x , objectPoints[i].y, objectPoints[i].z);
      cloudPt = cloudPt - TT; 

      transpose(RX, RX);
      transpose(RY, RY);
      transpose(RZ, RZ); 
      Mat cloudPtNew = RY*RX*RZ*cloudPt; 
   
      cv::Point3d pointslocal; 
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = cloudPtNew.at<double>(1, 0);
      pointslocal.z = cloudPtNew.at<double>(2, 0);

      std::cout<<fixed<<setprecision(10)<<pointslocal.x << "," <<pointslocal.y << "," <<pointslocal.z << std::endl; 

      cv::Point2d pointsuv(uimage, vimage);
 
      objPointsRight.push_back(pointslocal); 
      objPointsLocal.push_back(pointslocal);    
      imgPointsLocal.push_back(pointsuv);   
  }

  bool useExtrinsicGuess = true; 
  int  iterationsCount   = 5000; 
  double reprojectionError= 2.0; 
  double  minInliersCount = 0.94; 
  cv::Mat inliers ; 
  int flags              = cv::SOLVEPNP_ITERATIVE;

  std::cout <<"cv::SOLVEPNP_ITERATIVE: " << cv::SOLVEPNP_ITERATIVE << std::endl; 
  flags = 0; 

  cv::solvePnPRansac(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec, useExtrinsicGuess,iterationsCount,reprojectionError,minInliersCount,inliers,flags);

  //cv::solvePnP(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec);
 
  std::cout << "rvec: " << rvec << std::endl;
  std::cout << "tvec: " << tvec << std::endl;
 
  std::vector<cv::Point2d> projectedPoints; 
  cv::projectPoints(objPointsLocal, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);
 
  int rows = 2160; 
  int cols = 4096; 
  
  cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
 
  for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  {
    std::cout << "Image point: " << imagePoints[i].ptx<<"," << imagePoints[i].pty << " Projected to " << projectedPoints[i] << std::endl;
     circle(matShow, cv::Point(imagePoints[i].ptx, imagePoints[i].pty), 10,  Scalar(255,0,0), 2, 8); 
     circle(matShow, projectedPoints[i], 8,  Scalar(0,0,255), 2, 8); 
  }
  
  /*
  namedWindow("show", 2); 
  imshow("show", matShow);
  //imwrite("right.jpg", matShow);
  waitKey(200); 
  */

  return true; 
}

//active
bool IBD::calcLeftRT(Mat cameraMatrix, Mat distCoeffs, string path,std::vector<cv::Point3d>& objPointsLeft, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageLeft)
{
  // Read points
  std::vector<IBD::DataUVXYZ> imagePoints ;
  std::vector<cv::Point3d> objectPoints;

  std::string path_uv = path; 
  std::string path_xyz= path; 

  //return 0;

  IBD::loadUVleft(path_uv,  imagePoints);
  //return 0; 
  IBD::loadXYZleft(path_xyz,imagePoints, objectPoints);  
 
  std::cout << "There are " << imagePoints.size() << " imagePoints and " << objectPoints.size() << " objectPoints." << std::endl;

  //objPoints New 
  std::vector<cv::Point3d> objPointsLocal ;
  std::vector<cv::Point2d> imgPointsLocal ;

  assert(imagePoints.size() == objectPoints.size());

  for(int i = 0; i < imagePoints.size() ; i++)
  {
      double uimage = imagePoints[i].ptx; 
      double vimage = imagePoints[i].pty; 

      imageLeft.push_back(cv::Point2d(uimage,vimage)); 

      double Xshift = imagePoints[i].XShift; 
      double Yshift = imagePoints[i].YShift;
      double Zshift = imagePoints[i].ZShift;
      double rolls   = -imagePoints[i].Roll;
      double pitchs  = -imagePoints[i].Pitch;
      double yaws    = imagePoints[i].Yaw;

      //RT
      Mat RX = (Mat_<double>(3,3)<<1, 0,                  0,
                                   0, cos(pitchs*PI/180), sin(pitchs*PI/180),
                                   0, -sin(pitchs*PI/180), cos(pitchs*PI/180) );
        
      Mat RY = (Mat_<double>(3,3)<<cos(rolls*PI/180),0, -sin(rolls*PI/180),
                                   0,                1,   0,
                                   sin(rolls*PI/180),0, cos(rolls*PI/180) );
      Mat RZ = (Mat_<double>(3,3)<<cos(yaws*PI/180), sin(yaws*PI/180),  0,
                                   -sin(yaws*PI/180), cos(yaws*PI/180), 0,
                                   0, 0,                 1 ); 
        
      Mat TT = (Mat_<double>(3,1)<<Xshift, Yshift, Zshift); 
      Mat cloudPt = (Mat_<double>(3,1) << objectPoints[i].x , objectPoints[i].y, objectPoints[i].z);
      cloudPt = cloudPt - TT; 

      transpose(RX, RX);
      transpose(RY, RY);
      transpose(RZ, RZ); 
      Mat cloudPtNew = RY*RX*RZ*cloudPt; 
   
      cv::Point3d pointslocal; 
#ifndef IMU
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = cloudPtNew.at<double>(1, 0);
      pointslocal.z = cloudPtNew.at<double>(2, 0);
#else
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = -cloudPtNew.at<double>(2, 0);
      pointslocal.z = cloudPtNew.at<double>(1, 0);
#endif 
      std::cout<<fixed<<setprecision(10)<<"local point3d: "<<pointslocal.x << "," <<pointslocal.y << "," <<pointslocal.z << std::endl; 

      cv::Point2d pointsuv(uimage, vimage);
 
      objPointsLeft.push_back(pointslocal); 
      objPointsLocal.push_back(pointslocal);    
      imgPointsLocal.push_back(pointsuv);   
  }

  bool useExtrinsicGuess = true; 
  int  iterationsCount   = 50000; 
  double reprojectionError= 2.0; 
  double  minInliersCount = 0.96; 
  cv::Mat inliers ; 
  int flags              = cv::SOLVEPNP_ITERATIVE;

  std::cout <<"cv::SOLVEPNP_ITERATIVE: " << cv::SOLVEPNP_ITERATIVE << std::endl; 
  flags = 0; 

  cv::solvePnPRansac(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec, useExtrinsicGuess,iterationsCount,reprojectionError,minInliersCount,inliers,flags);

  //cv::solvePnP(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec);
 
  std::cout << "rvec: " << rvec << std::endl;
  std::cout << "tvec: " << tvec << std::endl;
 
  std::vector<cv::Point2d> projectedPoints; 
  cv::projectPoints(objPointsLocal, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);
 
  int rows = 2160; 
  int cols = 4096; 
  
  cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
 
  for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  {
    std::cout << "Image point: " << imagePoints[i].ptx<<"," << imagePoints[i].pty << " Projected to " << projectedPoints[i] << std::endl;
     circle(matShow, cv::Point(imagePoints[i].ptx, imagePoints[i].pty), 10,  Scalar(255,0,0), 2, 8); 
     circle(matShow, projectedPoints[i], 8,  Scalar(0,0,255), 2, 8); 
  }
  
  namedWindow("show", 2); 
  imshow("show", matShow);
  imwrite("firstLeft.jpg", matShow);
  waitKey(0); 
  return true; 
}
//active
bool IBD::calcRightRT(Mat cameraMatrix, Mat distCoeffs, string path, std::vector<cv::Point3d>& objPointsRight, Mat& rvec, Mat& tvec, std::vector<cv::Point2d> &imageRight)
{
  // Read points
  std::vector<IBD::DataUVXYZ> imagePoints ;
  std::vector<cv::Point3d> objectPoints;

  std::string path_uv = path; 
  std::string path_xyz= path; 

  IBD::loadUVright(path_uv,  imagePoints);
  IBD::loadXYZright(path_xyz,imagePoints, objectPoints);  
 
  std::cout << "There are " << imagePoints.size() << " imagePoints and " << objectPoints.size() << " objectPoints." << std::endl;

  //objPoints New 
  std::vector<cv::Point3d> objPointsLocal ;
  std::vector<cv::Point2d> imgPointsLocal ;

  assert(imagePoints.size() == objectPoints.size());

  for(int i = 0; i < imagePoints.size() ; i++)
  {
      double uimage = imagePoints[i].ptxr; 
      double vimage = imagePoints[i].ptyr; 
      imageRight.push_back(cv::Point2d(uimage,vimage)); 

      double Xshift = imagePoints[i].XShift; 
      double Yshift = imagePoints[i].YShift;
      double Zshift = imagePoints[i].ZShift;
      double rolls   = -imagePoints[i].Roll;
      double pitchs  = -imagePoints[i].Pitch;
      double yaws    = imagePoints[i].Yaw;

      //RT
      Mat RX = (Mat_<double>(3,3)<<1, 0,                  0,
                                   0, cos(pitchs*PI/180), -sin(pitchs*PI/180),
                                   0, sin(pitchs*PI/180), cos(pitchs*PI/180) );
        
      Mat RY = (Mat_<double>(3,3)<<cos(rolls*PI/180),0, sin(rolls*PI/180),
                                   0,                1,   0,
                                   -sin(rolls*PI/180),0, cos(rolls*PI/180) );
      Mat RZ = (Mat_<double>(3,3)<<cos(yaws*PI/180), -sin(yaws*PI/180),  0,
                                   sin(yaws*PI/180), cos(yaws*PI/180), 0,
                                   0, 0,                 1 ); 
        
      Mat TT = (Mat_<double>(3,1)<<Xshift, Yshift, Zshift); 
      Mat cloudPt = (Mat_<double>(3,1) << objectPoints[i].x , objectPoints[i].y, objectPoints[i].z);
      cloudPt = cloudPt - TT; 

      transpose(RX, RX);
      transpose(RY, RY);
      transpose(RZ, RZ); 
      Mat cloudPtNew = RY*RX*RZ*cloudPt; 
   
      cv::Point3d pointslocal; 
#ifndef IMU
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = cloudPtNew.at<double>(1, 0);
      pointslocal.z = cloudPtNew.at<double>(2, 0);
#else
      pointslocal.x = cloudPtNew.at<double>(0, 0);  
      pointslocal.y = -cloudPtNew.at<double>(2, 0);
      pointslocal.z = cloudPtNew.at<double>(1, 0);
#endif 
      std::cout<<fixed<<setprecision(10)<<"local point3d: "<<pointslocal.x << "," <<pointslocal.y << "," <<pointslocal.z << std::endl; 

      cv::Point2d pointsuv(uimage, vimage);
 
      objPointsRight.push_back(pointslocal); 
      objPointsLocal.push_back(pointslocal);    
      imgPointsLocal.push_back(pointsuv);   
  }

  bool useExtrinsicGuess = true; 
  int  iterationsCount   = 5000; 
  double reprojectionError= 2.0; 
  double  minInliersCount = 0.94; 
  cv::Mat inliers ; 
  int flags              = cv::SOLVEPNP_ITERATIVE;

  std::cout <<"cv::SOLVEPNP_ITERATIVE: " << cv::SOLVEPNP_ITERATIVE << std::endl; 
  flags = 1; 

  cv::solvePnPRansac(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec, useExtrinsicGuess,iterationsCount,reprojectionError,minInliersCount,inliers,flags);

  //cv::solvePnP(objPointsLocal, imgPointsLocal, cameraMatrix, distCoeffs, rvec, tvec);
 
  std::cout << "rvec: " << rvec << std::endl;
  std::cout << "tvec: " << tvec << std::endl;
 
  std::vector<cv::Point2d> projectedPoints; 
  cv::projectPoints(objPointsLocal, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);
 
  int rows = 2160; 
  int cols = 4096; 
  
  cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
 
  for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  {
    std::cout << "Image point: " << imagePoints[i].ptxr<<"," << imagePoints[i].ptyr << " Projected to " << projectedPoints[i] << std::endl;
     circle(matShow, cv::Point(imagePoints[i].ptx, imagePoints[i].pty), 10,  Scalar(255,0,0), 2, 8); 
     circle(matShow, projectedPoints[i], 8,  Scalar(0,0,255), 2, 8); 
  }
  
  namedWindow("show", 2); 
  imshow("show", matShow);
  imwrite("firstRight.jpg", matShow);
  waitKey(0); 

  return true; 
}
/*
    struct DataUV
    {
        double ptx; 
        double pty;
        double ptxr; 
        double ptyr;   
        double XGlobal; 
        double YGlobal; 
        double ZGlobal; 
        std::string name; 
        double XShift; 
        double YShift; 
        double ZShift;
        double Roll;
        double Pitch; 
        double Yaw; 
    };
*/
bool IBD::loadUVleft(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    //std::cout<<"path: " << path << std::endl; 

    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }
        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        IBD::DataUVXYZ point; 
        point.ptx   = atof(vecLine[0].c_str() ); 
        point.pty   = atof(vecLine[1].c_str() ); 
        point.ptxr  = atof(vecLine[2].c_str() ); 
        point.ptyr  = atof(vecLine[3].c_str() ); 
        point.XGlobal= atof(vecLine[4].c_str() ); 
        point.YGlobal= atof(vecLine[5].c_str() ); 
        point.ZGlobal= atof(vecLine[6].c_str() ); 
        point.name   = vecLine[7];
        point.XShift = pose[0];
        point.YShift = pose[1];
        point.ZShift = pose[2];
        point.Roll   = pose[3];
        point.Pitch  = pose[4];
        point.Yaw    = pose[5];

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}

bool IBD::loadUVleft(std::string path, std::vector<IBD::DataUV>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        IBD::DataUV point; 
        point.ptx   = atof(vecLine[0].c_str() ); 
        point.pty   = atof(vecLine[1].c_str() ); 
        point.Roll  = atof(vecLine[2].c_str() ); 
        point.Pitch  = atof(vecLine[3].c_str() ); 
        point.Yaw    = atof(vecLine[4].c_str() ); 
        point.XShift = atof(vecLine[5].c_str() ); 
        point.YShift = atof(vecLine[6].c_str() ); 
        point.ZShift = atof(vecLine[7].c_str() );

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}


bool IBD::loadUVleft(std::string path, std::vector<cv::Point2d>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        cv::Point2d point; 
        point.x = atof(vecLine[0].c_str() ); 
        point.y = atof(vecLine[1].c_str() ); 

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}

bool IBD::loadXYZleft(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints,  std::vector<cv::Point3d>& objectPoints)
{
    for(int i = 0 ; i < imagePoints.size() ; i++)
    {
        cv::Point3d pt3d(imagePoints[i].XGlobal, imagePoints[i].YGlobal, imagePoints[i].ZGlobal); 
        objectPoints.push_back(pt3d); 
    }
}
bool IBD::loadXYZleft(std::string path, std::vector<cv::Point3d>& objectPoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        cv::Point3d point3d; 
        point3d.x = atof(vecLine[0].c_str() ); 
        point3d.y = atof(vecLine[1].c_str() ); 
        point3d.z = atof(vecLine[2].c_str() ); 

        objectPoints.push_back(point3d); 
        //std::cout<< point3d.x << "," << point3d.y << "," <<point3d.z << std::endl; 
        vecLine.clear();          
    }
}

bool IBD::loadUVright(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        IBD::DataUVXYZ point; 
        point.ptx   = atof(vecLine[0].c_str() ); 
        point.pty   = atof(vecLine[1].c_str() ); 
        point.ptxr  = atof(vecLine[2].c_str() ); 
        point.ptyr  = atof(vecLine[3].c_str() ); 
        point.XGlobal= atof(vecLine[4].c_str() ); 
        point.YGlobal= atof(vecLine[5].c_str() ); 
        point.ZGlobal= atof(vecLine[6].c_str() ); 
        point.name   = vecLine[7];
        point.XShift = pose[0];
        point.YShift = pose[1];
        point.ZShift = pose[2];
        point.Roll   = pose[3];
        point.Pitch  = pose[4];
        point.Yaw    = pose[5];

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}

bool IBD::loadUVright(std::string path, std::vector<IBD::DataUV>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        IBD::DataUV point; 
        point.ptx   = atof(vecLine[0].c_str() ); 
        point.pty   = atof(vecLine[1].c_str() ); 
        point.Roll  = atof(vecLine[2].c_str() ); 
        point.Pitch  = atof(vecLine[3].c_str() ); 
        point.Yaw    = atof(vecLine[4].c_str() ); 
        point.XShift = atof(vecLine[5].c_str() ); 
        point.YShift = atof(vecLine[6].c_str() ); 
        point.ZShift = atof(vecLine[7].c_str() );

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}

bool IBD::loadUVright(std::string path, std::vector<cv::Point2d>& imagePoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        cv::Point2d point; 
        point.x = atof(vecLine[0].c_str() ); 
        point.y = atof(vecLine[1].c_str() ); 

        imagePoints.push_back(point); 
        //std::cout<< point.x << "," << point.y << std::endl; 
        vecLine.clear();          
    }
    
}

bool IBD::loadXYZright(std::string path, std::vector<cv::Point3d>& objectPoints)
{
    std::fstream file(path.c_str(), std::ios::in);
    if(!file.is_open() )
    {
        std::cout<<"file open failed!" << std::endl;
        return false; 
    }
    std::string strLine           = " "; 
    std::vector<std::string> vecLine; 

    while(1)
    {
        if( !getline(file, strLine) )
        {
            std::cout<< "file:" <<path << "  read finished!" << std::endl;
            break;  
        }

        boost::split(vecLine, strLine, boost::is_any_of(","));
      
        cv::Point3d point3d; 
        point3d.x = atof(vecLine[0].c_str() ); 
        point3d.y = atof(vecLine[1].c_str() ); 
        point3d.z = atof(vecLine[2].c_str() ); 

        objectPoints.push_back(point3d); 
        //std::cout<< point3d.x << "," << point3d.y << "," <<point3d.z << std::endl; 
        vecLine.clear();          
    }
}

bool IBD::loadXYZright(std::string path, std::vector<IBD::DataUVXYZ>& imagePoints, std::vector<cv::Point3d>& objectPoints)
{
    for(int i = 0 ; i < imagePoints.size() ; i++)
    {
        cv::Point3d pt3d(imagePoints[i].XGlobal, imagePoints[i].YGlobal, imagePoints[i].ZGlobal); 
        objectPoints.push_back(pt3d); 
    }
}

#endif 
