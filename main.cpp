#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
 
#include <iostream>
#include <fstream>
#include <string>
#include <DataBase.h>
#include <iomanip>
#include "CeresBa.h"
#include "shellCommand.h" 

std::vector<cv::Point2f> Generate2DPoints();
std::vector<cv::Point3f> Generate3DPoints();

int main( int argc, char* argv[])
{

  //[0]Load Parameters left
  cv::Mat cameraMatrixLeft(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixLeft);
  std::string common_patho = "../cameras_params/w/";
  std::string pathM1 = common_patho + "_M1.xml";

  cv::FileStorage fs(pathM1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M1"]>> cameraMatrixLeft;
      fs.release();
  } 
  std::cout << "cameraMatrixLeft: " << cameraMatrixLeft << std::endl;
 
  cv::Mat distCoeffsLeft;

  std::string pathD1 = common_patho + "_D1.xml";
  fs.open(pathD1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D1"]>> distCoeffsLeft;
      fs.release();
  } 
  std::cout << "distCoeffsLeft: " << distCoeffsLeft << std::endl;

  //[1]Load Parameters right
  cv::Mat cameraMatrixRight(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixRight);

  std::string pathM2 = common_patho + "_M2.xml";

  fs.open(pathM2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M2"]>> cameraMatrixRight;
      fs.release();
  } 
  std::cout << "cameraMatrixRight: " << cameraMatrixRight << std::endl;
 
  cv::Mat distCoeffsRight;

  std::string pathD2 = common_patho + "_D2.xml";
  fs.open(pathD2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D2"]>> distCoeffsRight;
      fs.release();
  } 
  std::cout << "distCoeffsRight: " << distCoeffsRight << std::endl;

  //[2] Read points 
  IBD_SD::CShellCommand shellCommandObj; 
  
  std::string path    = "../data/points/*.txt"; 
  std::string cmd_str = "ls " + path; 
 
  std::vector<std::string> veclists;
  shellCommandObj.getformatList(cmd_str, veclists);

  Problem problem;
  const int     numParameters  = 12; 
  double ratationleft[3] = {0};
  
  double Rotation[3]     = {0};
  double Translation[3] = {0};
  double Rotationl[3]   = {0}; 
  double Translationl[3]= {0}; 

  static int  numberFrame      = veclists.size(); 

  double Rotationlarray[numberFrame][3]={0};
  double Translationlarray[numberFrame][3]={0};

  //[3]calc R T 
  ///* 
  cv::Mat rvecl(3,1,cv::DataType<double>::type);
  cv::Mat tvecl(3,1,cv::DataType<double>::type);
  cv::Mat rvecr(3,1,cv::DataType<double>::type);
  cv::Mat tvecr(3,1,cv::DataType<double>::type);

  rvecl.at<double>(0,0)  = 0.00276471; 
  rvecl.at<double>(1,0)  = 0.0385392;
  rvecl.at<double>(2,0)  = -0.00296829;
  rvecr.at<double>(0,0) = -0.00127287; 
  rvecr.at<double>(1,0) = 0.0407736;
  rvecr.at<double>(2,0) = -0.00598665;

  tvecl.at<double>(0,0)  = 0.568852; 
  tvecl.at<double>(1,0)  = -0.288769;
  tvecl.at<double>(2,0)  = -0.0438256;
  tvecr.at<double>(0,0) = -0.532085; 
  tvecr.at<double>(1,0) = -0.280421;
  tvecr.at<double>(2,0) = -0.0285846;

  //*/

  //[3-1] R-T betweern two camera 
  Mat RRl, TTl, RRr,TTr, RRlinv, RR, TT;

  Rodrigues(rvecl, RRl);
  Rodrigues(rvecr, RRr);

  TTr = tvecr; 
  TTl = tvecl; 

  std::cout<<"Rodrigues over!" << std::endl; 

  RRlinv = RRl.inv(DECOMP_LU); 
  RR = RRr * RRlinv; 
  TT = TTr - RR * TTl; 

  std::cout<<"TT: " << TT << std::endl;  
  std::cout<<"RR: " << RR << std::endl; 

  double sy = sqrt(RR.at<double>(2,1) * RR.at<double>(2,1) + RR.at<double>(2,2) * RR.at<double>(2,2) );

  bool singular = (sy - 0.0000001 < 0); 
  if(!singular)
  {
      Rotation[0] = atan2(RR.at<double>(2,1), RR.at<double>(2,2) );
      Rotation[1] = atan2(-RR.at<double>(2,0), sy );
      Rotation[2] = atan2(RR.at<double>(1,0), RR.at<double>(0,0) );
  }
  else 
  {
      Rotation[0] = atan2(-RR.at<double>(1,2), RR.at<double>(1,1) );
      Rotation[1] = atan2(-RR.at<double>(2,0), sy );
      Rotation[2] = 0;      
  }
  Translation[0] = TT.at<double>(0,0);
  Translation[1] = TT.at<double>(1,0);
  Translation[2] = TT.at<double>(2,0);

  //[3-2]angle to rotation
    double theta[3] = {0};
    theta[0] = Rotation[0]; 
    theta[1] = Rotation[1]; 
    theta[2] = Rotation[2]; 

    // 计算旋转矩阵的X分量
    Mat R_x = (Mat_<double>(3,3) <<
               1,       0,              0,
               0,       cos(theta[0]),   -sin(theta[0]),
               0,       sin(theta[0]),   cos(theta[0])
               );
 
 
    // 计算旋转矩阵的Y分量
    Mat R_y = (Mat_<double>(3,3) <<
               cos(theta[1]),    0,      sin(theta[1]),
               0,               1,      0,
               -sin(theta[1]),   0,      cos(theta[1])
               );
 
 
    // 计算旋转矩阵的Z分量
    Mat R_z = (Mat_<double>(3,3) <<
               cos(theta[2]),    -sin(theta[2]),     0,
               sin(theta[2]),    cos(theta[2]),    0,
               0,               0,                  1);
 
 
    // 合并 
    Mat R = R_y * R_x ;
    R = R_z * R; 

    std::cout << "R: " << R << std::endl;  


    //[3-else]
    Mat RO, TO, RV; 

    std::string pathRO = common_patho + "_RO.xml";
    fs.open(pathRO, cv::FileStorage::READ);
    if(fs.isOpened())
    {
        fs["_R"]>> RO;
        fs.release();
    } 
    std::cout << "RO: " << RO << std::endl;   

    std::string pathTO = common_patho + "_TO.xml";
    fs.open(pathTO, cv::FileStorage::READ);
    if(fs.isOpened())
    {
        fs["_T"]>> TO;
        fs.release();
    } 
    std::cout << "TO: " << TO << std::endl; 

    Rodrigues(RO,RV);
    Translation[0] = TO.at<double>(0,0);
    Translation[1] = TO.at<double>(1,0);
    Translation[2] = TO.at<double>(2,0);
 
    Rotation[0] = RV.at<double>(0,0);
    Rotation[1] = RV.at<double>(1,0);
    Rotation[2] = RV.at<double>(2,0);

    std::vector<cv::Point3d> objPointsRightout;
    std::vector<cv::Point3d> objPointsLeftout;
    std::vector<cv::Point2d> imageRightout; 
    std::vector<cv::Point2d> imageLeftout;

    std::vector< std::vector<cv::Point3d> > vecobjPointsRightout;
    std::vector< std::vector<cv::Point3d> > vecobjPointsLeftout;
    std::vector< std::vector<cv::Point2d> > vecimageRightout; 
    std::vector< std::vector<cv::Point2d> > vecimageLeftout;

    int rowsleft = 2160; 
    int colsleft = 4096; 

    cv::Mat matLeftShow(rowsleft, colsleft, CV_8UC3, Scalar::all(100) );
    //matLeftShow.setTo(0);  

//////////////////////////////////////////////loop first 
for(int i = 0; i < numberFrame; i++)
{

  //[4]calc Rvec / Tvec 
  cv::Mat rvecleft(3,1,cv::DataType<double>::type);
  cv::Mat tvecleft(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsLeft;
  std::vector<cv::Point2d> imageLeft; 

  objPointsLeft.resize(0);
  imageLeft.resize(0);
 
  
  IBD::calcLeftRT(cameraMatrixLeft, distCoeffsLeft, veclists[i] ,objPointsLeft, rvecleft, tvecleft, imageLeft);
  //return 0;
  std::cout<<imageLeft[0].x << "," << imageLeft[0].y << "," << objPointsLeft[0].x << "," << objPointsLeft[0].y<< "," << objPointsLeft[0].z << std::endl; 

  for(int k = 0; k < imageLeft.size() ; k++)
  {
      cv::Point2d pt(imageLeft[k].x, imageLeft[k].y); 
      cv::circle(matLeftShow, pt, 8, Scalar(0,255,0), 2, 8); 
  }

  cv::Mat rvecright(3,1,cv::DataType<double>::type);
  cv::Mat tvecright(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsRight;
  std::vector<cv::Point2d> imageRight; 

  objPointsRight.resize(0);
  imageRight.resize(0);
  imageLeftout.resize(0); 
  imageRightout.resize(0);
  objPointsRightout.resize(0);
  objPointsLeftout.resize(0); 

 
  IBD::calcRightRT(cameraMatrixRight, distCoeffsRight, veclists[i] ,objPointsRight, rvecright, tvecright, imageRight);
  std::cout<<imageRight[0].x << "," << imageRight[0].y << "," << objPointsRight[0].x << "," << objPointsRight[0].y<< "," << objPointsRight[0].z << std::endl; 

  std::cout<<"list: " << veclists[i] << std::endl;
  imageLeftout      = imageLeft;  
  imageRightout     = imageRight; 
  objPointsRightout = objPointsRight;
  objPointsLeftout  = objPointsLeft; 

  vecimageLeftout.push_back(imageLeftout); 
  vecimageRightout.push_back(imageRightout);
  vecobjPointsRightout.push_back(objPointsRightout);
  vecobjPointsLeftout.push_back(objPointsLeftout);
 
   
  std::cout<<"objPointsLeft.size() : " << objPointsLeft.size() << std::endl; 
  std::cout<<"objPointsRight.size() : " << objPointsRight.size() << std::endl;
  assert(objPointsLeft.size() == objPointsRight.size()); 


  //return 0; 
  //[5]calc R & T 
/* 
  rvecleft.at<double>(0,0)  = 0.00276471; 
  rvecleft.at<double>(0,1)  = 0.0385392;
  rvecleft.at<double>(0,2)  = -0.00296829;
  rvecright.at<double>(0,0) = -0.00127287; 
  rvecright.at<double>(0,1) = 0.0407736;
  rvecright.at<double>(0,2) = -0.00598665;

  tvecleft.at<double>(0,0)  = 0.568852; 
  tvecleft.at<double>(1,0)  = -0.288769;
  tvecleft.at<double>(2,0)  = -0.0438256;
  tvecright.at<double>(0,0) = -0.532085; 
  tvecright.at<double>(1,0) = -0.280421;
  tvecright.at<double>(2,0) = -0.0285846;
*/
  Mat Rl, Tl, Rr,Tr, Rlinv, R, T;

  Rodrigues(rvecleft, Rl);
  Rodrigues(rvecright, Rr);

  Tr = tvecright; 
  Tl = tvecleft; 

  std::cout<<"Rodrigues over!" << std::endl; 

  Rlinv = Rl.inv(DECOMP_LU); 
  R = Rr * Rlinv; 
  T = Tr - R * Tl; 

  std::cout<<"R: " << R << "," << std::endl; 
  std::cout<< "T: " << T << std::endl; 
  std::cout<<"Rr: " << Rr << std::endl; 
  std::cout<< "Tr: " << Tr << std::endl;
  std::cout<<"Rl: " << Rl << std::endl; 
  std::cout<< "Tl: " << Tl << std::endl; 

  double syl = sqrt(Rl.at<double>(2,1) * Rl.at<double>(2,1) + Rl.at<double>(2,2) * Rl.at<double>(2,2) );

  bool singularl = (syl - 0.0000001 < 0); 
  if(!singularl)
  {
      Rotationl[0] = atan2(Rl.at<double>(2,1), Rl.at<double>(2,2) );
      Rotationl[1] = atan2(-Rl.at<double>(2,0), syl );
      Rotationl[2] = atan2(Rl.at<double>(1,0), Rl.at<double>(0,0) );
  }
  else 
  {
      Rotationl[0] = atan2(-Rl.at<double>(1,2), Rl.at<double>(1,1) );
      Rotationl[1] = atan2(-Rl.at<double>(2,0), syl );
      Rotationl[2] = 0;      
  }

/*-----R-T-FOR-every-frame-------
  double sys = sqrt(R.at<double>(2,1) * R.at<double>(2,1) + R.at<double>(2,2) * R.at<double>(2,2) );

  bool singulars = (sys - 0.0000001 < 0); 
  if(!singulars)
  {
      Rotation[0] = atan2(R.at<double>(2,1), R.at<double>(2,2) );
      Rotation[1] = atan2(-R.at<double>(2,0), sys );
      Rotation[2] = atan2(R.at<double>(1,0), R.at<double>(0,0) );
  }
  else 
  {
      Rotation[0] = atan2(-R.at<double>(1,2), R.at<double>(1,1) );
      Rotation[1] = atan2(-R.at<double>(2,0), sys );
      Rotation[2] = 0;      
  }
  Translation[0] = T.at<double>(0,0);
  Translation[1] = T.at<double>(1,0);
  Translation[2] = T.at<double>(2,0);
*/

/*-----------------------------make sure-----------------
    //angle to rotation
    double theta[3] = {0};
    theta[0] = Rotationl[0]; 
    theta[1] = Rotationl[1]; 
    theta[2] = Rotationl[2]; 

    // 计算旋转矩阵的X分量
    Mat R_x = (Mat_<double>(3,3) <<
               1,       0,              0,
               0,       cos(theta[0]),   -sin(theta[0]),
               0,       sin(theta[0]),   cos(theta[0])
               );
 
 
    // 计算旋转矩阵的Y分量
    Mat R_y = (Mat_<double>(3,3) <<
               cos(theta[1]),    0,      sin(theta[1]),
               0,               1,      0,
               -sin(theta[1]),   0,      cos(theta[1])
               );
 
 
    // 计算旋转矩阵的Z分量
    Mat R_z = (Mat_<double>(3,3) <<
               cos(theta[2]),    -sin(theta[2]),     0,
               sin(theta[2]),    cos(theta[2]),    0,
               0,               0,                  1);
 
 
    // 合并 
    Mat Rll = R_z * R_y * R_x;
    std::cout << "Rll: " << Rll << std::endl;  
//-------------------------make sure end----------------*/
    
  std::string pathRR = common_patho + "_RR.xml";
  cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite << "_R" << R ;
      fswrite.release();
  } 

  std::string pathTT = common_patho + "_TT.xml";
  fswrite.open(pathTT, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite <<"_T" << T ;
      fswrite.release();
  }

  //[6]calc new points 
  std::vector<cv::Point3d> objPointsLeftCam;
  std::vector<cv::Point3d> objPointsRightCam;
  ////////loop second////////////
  for(int i = 0 ; i < 1; i++) //objPointsLeft.size(); i++)
  {
      cv::Point3d ptleftworld = objPointsLeft[i]; 
      Mat PtWorldl = (Mat_<double>(3,1)<<ptleftworld.x, ptleftworld.y, ptleftworld.z);
      
      PtWorldl     = Rl * PtWorldl + Tl;  
      ptleftworld.x= PtWorldl.at<double>(0,0); 
      ptleftworld.y= PtWorldl.at<double>(1,0); 
      ptleftworld.z= PtWorldl.at<double>(2,0); 

      cv::Point3d ptrightworld = objPointsRight[i]; 
      Mat PtWorldr = (Mat_<double>(3,1)<<ptrightworld.x, ptrightworld.y, ptrightworld.z);
      
      PtWorldr      = Rr * PtWorldr + Tr;   
      ptrightworld.x= PtWorldr.at<double>(0,0); 
      ptrightworld.y= PtWorldr.at<double>(1,0); 
      ptrightworld.z= PtWorldr.at<double>(2,0);

      std::cout<< ptleftworld.x << "," << ptleftworld.y << "," << ptleftworld.z <<"," << ptrightworld.x << "," << ptrightworld.y << "," << ptrightworld.z << std::endl; 

      cv::Point3d ptrightworldnew; 
      PtWorldl = R * PtWorldl + T; 
      ptrightworldnew.x = PtWorldl.at<double>(0,0);
      ptrightworldnew.y = PtWorldl.at<double>(1,0);
      ptrightworldnew.z = PtWorldl.at<double>(2,0);

      double diffx = (ptrightworldnew.x- ptrightworld.x); 
      double diffy = (ptrightworldnew.y- ptrightworld.y);
      double diffz = (ptrightworldnew.z- ptrightworld.z);
      double dist  = diffx * diffx + diffy * diffy + diffz * diffz; 
      dist  = sqrt(dist);
      std::cout<<ptrightworldnew.x << "," << ptrightworldnew.y << "," << ptrightworldnew.z << "," << ptrightworld.x << "," << ptrightworld.y <<"," << ptrightworld.z << "," << dist<< std::endl; 

  }
    //return 0; 
    //[7] ceres 
    std::cout<<"R: " <<  R.type() << std::endl; 

    Rotationlarray[i][0] = Rotationl[0];
    Rotationlarray[i][1] = Rotationl[1];
    Rotationlarray[i][2] = Rotationl[2];

    Translationlarray[i][0] = Tl.at<double>(0,0);
    Translationlarray[i][1] = Tl.at<double>(1,0);
    Translationlarray[i][2] = Tl.at<double>(2,0);

    int indexm = i; 
    ////////loop third///////////
    for(int i = 0; i < objPointsLeft.size(); i++)
    {
        StereoCalibrationResiduals* pResidualX = new StereoCalibrationResiduals(
                                                    imageLeft[i].x, 
                                                    imageLeft[i].y,
                                                    objPointsLeft[i].x, 
                                                    objPointsLeft[i].y, 
                                                    objPointsLeft[i].z, 
                                                    imageRight[i].x, 
                                                    imageRight[i].y, 
                                                    objPointsRight[i].x, 
                                                    objPointsRight[i].y, 
                                                    objPointsRight[i].z, 
                                                    cameraMatrixLeft.at<double>(0,0),
                                                    cameraMatrixLeft.at<double>(0,2),
                                                    cameraMatrixLeft.at<double>(1,1), 
                                                    cameraMatrixLeft.at<double>(1,2),
                                                    distCoeffsLeft.at<double>(0,0), 
                                                    distCoeffsLeft.at<double>(0,1), 
                                                    distCoeffsLeft.at<double>(0,2), 
                                                    distCoeffsLeft.at<double>(0,3),
                                                    distCoeffsLeft.at<double>(0,4), 
                                                    cameraMatrixRight.at<double>(0,0),
                                                    cameraMatrixRight.at<double>(0,2),
                                                    cameraMatrixRight.at<double>(1,1),
                                                    cameraMatrixRight.at<double>(1,2),
                                                    distCoeffsRight.at<double>(0,0),
                                                    distCoeffsRight.at<double>(0,1),
                                                    distCoeffsRight.at<double>(0,2),
                                                    distCoeffsRight.at<double>(0,3),
                                                    distCoeffsRight.at<double>(0,4));

        CostFunction* cost_function = new AutoDiffCostFunction<StereoCalibrationResiduals, 1, 3 , 3, 3 ,3>(pResidualX);
        problem.AddResidualBlock(cost_function, new CauchyLoss(5), Rotation, Translation, Rotationlarray[indexm], Translationlarray[indexm]);
         
    }
}
    //return 0; 
    Solver::Options  options;
    options.max_num_iterations = 10000; 
    options.linear_solver_type = ceres::DENSE_SCHUR; //ceres::DENSE_QR; 
    options.minimizer_progress_to_stdout = true; 
 
    Solver::Summary summary; 
    Solve(options, &problem, &summary);

    std::cout<<"summary: " << summary.BriefReport() << std::endl; 

    for(int i = 0 ; i < 3 ; i++)
    {
         std::cout<<"Translation: " << Translation[i] << std::endl; 
    }
    for(int i = 0 ; i < 3 ; i++)
    {
        std::cout<<"Rotation: " << Rotation[i] << std::endl;
    }

    //R angle to rotation
    double thetas[3] = {0};
    thetas[0] = Rotation[0]; 
    thetas[1] = Rotation[1]; 
    thetas[2] = Rotation[2]; 

    //R X分量
    Mat Rx = (Mat_<double>(3,3) <<
               1,       0,              0,
               0,       cos(thetas[0]),   -sin(thetas[0]),
               0,       sin(thetas[0]),   cos(thetas[0])
               );
 
 
    //R Y分量
    Mat Ry = (Mat_<double>(3,3) <<
               cos(thetas[1]),    0,      sin(thetas[1]),
               0,               1,      0,
               -sin(thetas[1]),   0,      cos(thetas[1])
               );
 
 
    //R Z分量
    Mat Rz = (Mat_<double>(3,3) <<
               cos(thetas[2]),    -sin(thetas[2]),     0,
               sin(thetas[2]),    cos(thetas[2]),    0,
               0,               0,                  1);
 

    Mat matR = Rz * Ry * Rx;
    cv::Mat matT(3, 1, CV_64FC1, Translation); 

    std::cout << "matR: " << matR << std::endl; 
    std::string pathR = common_patho + "_R.xml";
    cv::FileStorage fswrite(pathR, cv::FileStorage::WRITE);
    if(fswrite.isOpened())
    {
        fswrite << "_R" << matR ;
        fswrite.release();
    } 

    std::string pathT = common_patho + "_T.xml";
    fswrite.open(pathT, cv::FileStorage::WRITE);
    if(fswrite.isOpened())
    {
        fswrite <<"_T" << matT ;
        fswrite.release();
    }

    cv::namedWindow("matShow", 1); 
    cv::imshow("matShow", matLeftShow); 
    imwrite("leftShow.jpg", matLeftShow);
    cv::waitKey(100);  

    //return 0; 
    //file storage 
    fstream file_obj("obj3d.dat",ios::out); 
    for(int k = 0; k < vecobjPointsRightout.size() ; k++)
    {
        //Rl angle to rotation Rl
        double thetal[3] = {0};
        thetal[0] = Rotationlarray[k][0]; 
        thetal[1] = Rotationlarray[k][1]; 
        thetal[2] = Rotationlarray[k][2]; 

        std::cout<< "Rotationlarray[0][0]: " << Rotationlarray[k][0]<< std::endl; 
        std::cout<< "Rotationlarray[0][1]: " << Rotationlarray[k][1]<< std::endl;
        std::cout<< "Rotationlarray[0][2]: " << Rotationlarray[k][2]<< std::endl;

        std::cout<< "Translationlarray[0][0]: " << Translationlarray[k][0]<< std::endl; 
        std::cout<< "Translationlarray[0][1]: " << Translationlarray[k][1]<< std::endl;
        std::cout<< "Translationlarray[0][2]: " << Translationlarray[k][2]<< std::endl;

        // RL X分量
        Mat Rxl = (Mat_<double>(3,3) <<
               1,       0,              0,
               0,       cos(thetal[0]),   -sin(thetal[0]),
               0,       sin(thetal[0]),   cos(thetal[0])
               );
 
 
        // RL Y分量
        Mat Ryl = (Mat_<double>(3,3) <<
               cos(thetal[1]),    0,      sin(thetal[1]),
               0,               1,      0,
               -sin(thetal[1]),   0,      cos(thetal[1])
               );
 
 
        // RL Z分量
        Mat Rzl = (Mat_<double>(3,3) <<
               cos(thetal[2]),    -sin(thetal[2]),     0,
               sin(thetal[2]),    cos(thetal[2]),    0,
               0,               0,                  1);
 
        //Rr
        Mat matRl = Rzl * Ryl * Rxl;
        Mat matTl(3, 1, CV_64FC1, Translationlarray[k]);

        std::cout << "matRl: " << matRl << std::endl; 
 
        //return 0; 
        //right show 
        int rows = 2160; 
        int cols = 4096; 
 
    
        Mat matTr;    
        Mat matRr = matR * matRl;

        matTr = matT + matR * matTl;
 
        Mat vecRr; 
        Rodrigues(matRr, vecRr); 


        objPointsRightout.clear(); 
        objPointsRightout = vecobjPointsRightout[k]; 
 
        cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
        std::vector<cv::Point2d> projectedPoints;

        cv::projectPoints(objPointsRightout, vecRr, matTr, cameraMatrixRight, distCoeffsRight, projectedPoints);
        
        imageRightout.clear(); 
        imageRightout = vecimageRightout[k];
        for(unsigned int i = 0; i < projectedPoints.size(); ++i)
        {
            file_obj<<fixed<<setprecision(10) <<  objPointsRightout[i].x << "," << objPointsRightout[i].y << "," << objPointsRightout[i].z << std::endl; 

            std::cout << "Image point: " << imageRightout[i].x<<"," << imageRightout[i].y << " Projected to " << projectedPoints[i] << std::endl;
            circle(matShow, cv::Point(imageRightout[i].x, imageRightout[i].y), 10,  Scalar(255,0,0), 2, 8); 
            circle(matShow, projectedPoints[i], 8,  Scalar(0,0,255), 2, 8); 
        }
  
        namedWindow("show", 2); 
        imwrite("rightfirst.jpg", matShow);
        imshow("show", matShow);
        waitKey(0); 

    //left show 
    Mat vecRl; 
    Rodrigues(matRl, vecRl);  

    cv::Mat matShowl(rows, cols, CV_8UC3, Scalar::all(100));
    std::vector<cv::Point2d> projectedPointsl;

    objPointsLeftout.clear(); 
    objPointsLeftout = vecobjPointsLeftout[k];

    cv::projectPoints(objPointsLeftout, vecRl, matTl, cameraMatrixLeft, distCoeffsLeft, projectedPointsl);
 

    imageLeftout.clear(); 
    imageLeftout = vecimageLeftout[k]; 
    for(unsigned int i = 0; i < projectedPoints.size(); ++i)
    {
        file_obj<<fixed<<setprecision(10) <<  objPointsLeftout[i].x << "," << objPointsLeftout[i].y << "," << objPointsLeftout[i].z << std::endl;

      std::cout << "Image point: " << imageLeftout[i].x<<"," << imageLeftout[i].y << " Projected to " << projectedPointsl[i] << std::endl;
       circle(matShowl, cv::Point(imageLeftout[i].x, imageLeftout[i].y), 10,  Scalar(255,0,0), 2, 8); 
       circle(matShowl, projectedPointsl[i], 8,  Scalar(0,0,255), 2, 8); 
    }
  
    namedWindow("show", 2); 
    imshow("show", matShowl);
    imwrite("leftfirst.jpg", matShowl);
    waitKey(0); 

    }

    return 0; 

}
 
int main_Rl_Tl_R_T( int argc, char* argv[])
{

  //[0]Load Parameters left
  cv::Mat cameraMatrixLeft(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixLeft);
  std::string common_patho = "../cameras_params/w/";
  std::string pathM1 = common_patho + "_M1.xml";

  cv::FileStorage fs(pathM1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M1"]>> cameraMatrixLeft;
      fs.release();
  } 
  std::cout << "cameraMatrixLeft: " << cameraMatrixLeft << std::endl;
 
  cv::Mat distCoeffsLeft;

  std::string pathD1 = common_patho + "_D1.xml";
  fs.open(pathD1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D1"]>> distCoeffsLeft;
      fs.release();
  } 
  std::cout << "distCoeffsLeft: " << distCoeffsLeft << std::endl;

  //[1]Load Parameters right
  cv::Mat cameraMatrixRight(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixRight);

  std::string pathM2 = common_patho + "_M2.xml";

  fs.open(pathM2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M2"]>> cameraMatrixRight;
      fs.release();
  } 
  std::cout << "cameraMatrixRight: " << cameraMatrixRight << std::endl;
 
  cv::Mat distCoeffsRight;

  std::string pathD2 = common_patho + "_D2.xml";
  fs.open(pathD2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D2"]>> distCoeffsRight;
      fs.release();
  } 
  std::cout << "distCoeffsRight: " << distCoeffsRight << std::endl;

  //[2] Read points 
  IBD_SD::CShellCommand shellCommandObj; 
  
  std::string path    = "/home/ibd01/DVISION/src/tools/calibration/StereoCalibration/data/points/*.txt"; 
  std::string cmd_str = "ls " + path; 
 
  std::vector<std::string> veclists;
  shellCommandObj.getformatList(cmd_str, veclists);

  Problem problem;
  const int     numParameters  = 12; 
  double ratationleft[3] = {0};
  
  std::vector< std::vector<double> > vvecStereoCalibrateParameters;
  vvecStereoCalibrateParameters.resize(veclists.size()); 

  double StereoCalibrateParameters[12] = {0};
  double Rotation[9]     = {0};
  double Translation[3] = {0};
  double Rotationl[9]   = {0}; 
  double Translationl[3]= {0}; 

  static int  numberFrame      = veclists.size(); 

  double Rotationlarray[numberFrame][9]={0};
  double Translationlarray[numberFrame][3]={0};

  //[3]calc R T 
  ///* 
  cv::Mat rvecl(3,1,cv::DataType<double>::type);
  cv::Mat tvecl(3,1,cv::DataType<double>::type);
  cv::Mat rvecr(3,1,cv::DataType<double>::type);
  cv::Mat tvecr(3,1,cv::DataType<double>::type);

  rvecl.at<double>(0,0)  = 0.00276471; 
  rvecl.at<double>(1,0)  = 0.0385392;
  rvecl.at<double>(2,0)  = -0.00296829;
  rvecr.at<double>(0,0) = -0.00127287; 
  rvecr.at<double>(1,0) = 0.0407736;
  rvecr.at<double>(2,0) = -0.00598665;

  tvecl.at<double>(0,0)  = 0.568852; 
  tvecl.at<double>(1,0)  = -0.288769;
  tvecl.at<double>(2,0)  = -0.0438256;
  tvecr.at<double>(0,0) = -0.532085; 
  tvecr.at<double>(1,0) = -0.280421;
  tvecr.at<double>(2,0) = -0.0285846;

  //*/
  Mat RRl, TTl, RRr,TTr, RRlinv, RR, TT;

  Rodrigues(rvecl, RRl);
  Rodrigues(rvecr, RRr);

  TTr = tvecr; 
  TTl = tvecl; 

  std::cout<<"Rodrigues over!" << std::endl; 

  RRlinv = RRl.inv(DECOMP_LU); 
  RR = RRr * RRlinv; 
  TT = TTr - RR * TTl; 
  std::cout<<"TT: " << TT << std::endl;  
  
  //return 0; 
for(int i = 0; i < numberFrame; i++)
{

  //[4]calc Rvec / Tvec 
  cv::Mat rvecleft(1,3,cv::DataType<double>::type);
  cv::Mat tvecleft(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsLeft;
  std::vector<cv::Point2d> imageLeft; 

  objPointsLeft.resize(0);
  imageLeft.resize(0);
 
  
  IBD::calcLeftRT(cameraMatrixLeft, distCoeffsLeft, veclists[i] ,objPointsLeft, rvecleft, tvecleft, imageLeft);
  //return 0;
  std::cout<<imageLeft[0].x << "," << imageLeft[0].y << "," << objPointsLeft[0].x << "," << objPointsLeft[0].y<< "," << objPointsLeft[0].z << std::endl; 

  cv::Mat rvecright(1,3,cv::DataType<double>::type);
  cv::Mat tvecright(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsRight;
  std::vector<cv::Point2d> imageRight; 

  objPointsRight.resize(0);
  imageRight.resize(0);
 
  IBD::calcRightRT(cameraMatrixRight, distCoeffsRight, veclists[i] ,objPointsRight, rvecright, tvecright, imageRight);
  std::cout<<imageRight[0].x << "," << imageRight[0].y << "," << objPointsRight[0].x << "," << objPointsRight[0].y<< "," << objPointsRight[0].z << std::endl; 

  std::cout<<"list: " << veclists[i] << std::endl; 
   
  std::cout<<"objPointsLeft.size() : " << objPointsLeft.size() << std::endl; 
  std::cout<<"objPointsRight.size() : " << objPointsRight.size() << std::endl;
  assert(objPointsLeft.size() == objPointsRight.size()); 


  //return 0; 
  //[5]calc R & T 
/* 
  rvecleft.at<double>(0,0)  = 0.00276471; 
  rvecleft.at<double>(0,1)  = 0.0385392;
  rvecleft.at<double>(0,2)  = -0.00296829;
  rvecright.at<double>(0,0) = -0.00127287; 
  rvecright.at<double>(0,1) = 0.0407736;
  rvecright.at<double>(0,2) = -0.00598665;

  tvecleft.at<double>(0,0)  = 0.568852; 
  tvecleft.at<double>(1,0)  = -0.288769;
  tvecleft.at<double>(2,0)  = -0.0438256;
  tvecright.at<double>(0,0) = -0.532085; 
  tvecright.at<double>(1,0) = -0.280421;
  tvecright.at<double>(2,0) = -0.0285846;
*/
  Mat Rl, Tl, Rr,Tr, Rlinv, R, T;

  Rodrigues(rvecleft, Rl);
  Rodrigues(rvecright, Rr);

  Tr = tvecright; 
  Tl = tvecleft; 

  std::cout<<"Rodrigues over!" << std::endl; 

  Rlinv = Rl.inv(DECOMP_LU); 
  R = Rr * Rlinv; 
  T = Tr - R * Tl; 

  std::cout<<"R: " << R << "," << std::endl; 
  std::cout<< "T: " << T << std::endl; 
  std::cout<<"Rr: " << Rr << std::endl; 
  std::cout<< "Tr: " << Tr << std::endl;
  std::cout<<"Rl: " << Rl << std::endl; 
  std::cout<< "Tl: " << Tl << std::endl; 
   
  std::string pathRR = common_patho + "_RR.xml";
  cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite << "_R" << R ;
      fswrite.release();
  } 

  std::string pathTT = common_patho + "_TT.xml";
  fswrite.open(pathTT, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite <<"_T" << T ;
      fswrite.release();
  }

  //[6]calc new points 
  std::vector<cv::Point3d> objPointsLeftCam;
  std::vector<cv::Point3d> objPointsRightCam;

  for(int i = 0 ; i < 1; i++) //objPointsLeft.size(); i++)
  {
      cv::Point3d ptleftworld = objPointsLeft[i]; 
      Mat PtWorldl = (Mat_<double>(3,1)<<ptleftworld.x, ptleftworld.y, ptleftworld.z);
      
      PtWorldl     = Rl * PtWorldl + Tl;  
      ptleftworld.x= PtWorldl.at<double>(0,0); 
      ptleftworld.y= PtWorldl.at<double>(1,0); 
      ptleftworld.z= PtWorldl.at<double>(2,0); 

      cv::Point3d ptrightworld = objPointsRight[i]; 
      Mat PtWorldr = (Mat_<double>(3,1)<<ptrightworld.x, ptrightworld.y, ptrightworld.z);
      
      PtWorldr      = Rr * PtWorldr + Tr;   
      ptrightworld.x= PtWorldr.at<double>(0,0); 
      ptrightworld.y= PtWorldr.at<double>(1,0); 
      ptrightworld.z= PtWorldr.at<double>(2,0);

      std::cout<< ptleftworld.x << "," << ptleftworld.y << "," << ptleftworld.z <<"," << ptrightworld.x << "," << ptrightworld.y << "," << ptrightworld.z << std::endl; 

      cv::Point3d ptrightworldnew; 
      PtWorldl = R * PtWorldl + T; 
      ptrightworldnew.x = PtWorldl.at<double>(0,0);
      ptrightworldnew.y = PtWorldl.at<double>(1,0);
      ptrightworldnew.z = PtWorldl.at<double>(2,0);

      double diffx = (ptrightworldnew.x- ptrightworld.x); 
      double diffy = (ptrightworldnew.y- ptrightworld.y);
      double diffz = (ptrightworldnew.z- ptrightworld.z);
      double dist  = diffx * diffx + diffy * diffy + diffz * diffz; 
      dist  = sqrt(dist);
      std::cout<<ptrightworldnew.x << "," << ptrightworldnew.y << "," << ptrightworldnew.z << "," << ptrightworld.x << "," << ptrightworld.y <<"," << ptrightworld.z << "," << dist<< std::endl; 

  }
    //return 0; 
    //[7] ceres 
    std::cout<<"R: " <<  R.type() << std::endl; 

    Rotation[0] = RR.at<double>(0,0);
    Rotation[1] = RR.at<double>(0,1);
    Rotation[2] = RR.at<double>(0,2);
    Rotation[3] = RR.at<double>(1,0);
    Rotation[4] = RR.at<double>(1,1);
    Rotation[5] = RR.at<double>(1,2);
    Rotation[6] = RR.at<double>(2,0);
    Rotation[7] = RR.at<double>(2,1);
    Rotation[8] = RR.at<double>(2,2);

    Translation[0] = TT.at<double>(0,0);
    Translation[1] = TT.at<double>(1,0);
    Translation[2] = TT.at<double>(2,0);
/*
    Rotationl[0] = Rl.at<double>(0,0);
    Rotationl[1] = Rl.at<double>(0,1);
    Rotationl[2] = Rl.at<double>(0,2);
    Rotationl[3] = Rl.at<double>(1,0);
    Rotationl[4] = Rl.at<double>(1,1);
    Rotationl[5] = Rl.at<double>(1,2);
    Rotationl[6] = Rl.at<double>(2,0);
    Rotationl[7] = Rl.at<double>(2,1);
    Rotationl[8] = Rl.at<double>(2,2);

    Translationl[0] = Tl.at<double>(0,0);
    Translationl[1] = Tl.at<double>(1,0);
    Translationl[2] = Tl.at<double>(2,0);
*/
    Rotationlarray[i][0] = Rl.at<double>(0,0);
    Rotationlarray[i][1] = Rl.at<double>(0,1);
    Rotationlarray[i][2] = Rl.at<double>(0,2);
    Rotationlarray[i][3] = Rl.at<double>(1,0);
    Rotationlarray[i][4] = Rl.at<double>(1,1);
    Rotationlarray[i][5] = Rl.at<double>(1,2);
    Rotationlarray[i][6] = Rl.at<double>(2,0);
    Rotationlarray[i][7] = Rl.at<double>(2,1);
    Rotationlarray[i][8] = Rl.at<double>(2,2);

    Translationlarray[i][0] = Tl.at<double>(0,0);
    Translationlarray[i][1] = Tl.at<double>(1,0);
    Translationlarray[i][2] = Tl.at<double>(2,0);

    int indexm = i; 

    for(int i = 0; i < objPointsLeft.size(); i++)
    {
        //std::cout<<imageLeft[i].x <<"," << imageLeft[i].y << "," <<imageRight[i].x <<"," << imageRight[i].y << std::endl;   
        //return 0; 
        StereoCalibrationResiduals* pResidualX = new StereoCalibrationResiduals(
                                                    imageLeft[i].x, 
                                                    imageLeft[i].y,
                                                    objPointsLeft[i].x, 
                                                    objPointsLeft[i].y, 
                                                    objPointsLeft[i].z, 
                                                    imageRight[i].x, 
                                                    imageRight[i].y, 
                                                    objPointsRight[i].x, 
                                                    objPointsRight[i].y, 
                                                    objPointsRight[i].z, 
                                                    cameraMatrixLeft.at<double>(0,0),
                                                    cameraMatrixLeft.at<double>(0,2),
                                                    cameraMatrixLeft.at<double>(1,1), 
                                                    cameraMatrixLeft.at<double>(1,2),
                                                    distCoeffsLeft.at<double>(0,0), 
                                                    distCoeffsLeft.at<double>(0,1), 
                                                    distCoeffsLeft.at<double>(0,2), 
                                                    distCoeffsLeft.at<double>(0,3),
                                                    distCoeffsLeft.at<double>(0,4), 
                                                    cameraMatrixRight.at<double>(0,0),
                                                    cameraMatrixRight.at<double>(0,2),
                                                    cameraMatrixRight.at<double>(1,1),
                                                    cameraMatrixRight.at<double>(1,2),
                                                    distCoeffsRight.at<double>(0,0),
                                                    distCoeffsRight.at<double>(0,1),
                                                    distCoeffsRight.at<double>(0,2),
                                                    distCoeffsRight.at<double>(0,3),
                                                    distCoeffsRight.at<double>(0,4));

        CostFunction* cost_function = new AutoDiffCostFunction<StereoCalibrationResiduals, 1, 9 , 3, 9,3>(pResidualX);
        problem.AddResidualBlock(cost_function, new CauchyLoss(5), Rotation, Translation, Rotationlarray[indexm], Translationlarray[indexm]);
         
    }
}
    //return 0; 
    Solver::Options  options;
    options.max_num_iterations = 1000; 
    options.linear_solver_type = ceres::DENSE_QR; 
    options.minimizer_progress_to_stdout = true; 
 
    Solver::Summary summary; 
    Solve(options, &problem, &summary);
    for(int i = 0 ; i < 3 ; i++)
    {
         std::cout<<"Translation: " << Translation[i] << std::endl; 
    }
    for(int i = 0 ; i < 3 ; i++)
    {
        for(int j = 0; j < 3 ; j++)
        {
            std::cout<<"Rotation: " << Rotation[i * 3 + j] << std::endl;
        }
    }

    cv::Mat matR(3, 3, CV_64FC1, Rotation); 
    cv::Mat matT(3, 1, CV_64FC1, Translation); 

    std::string pathR = common_patho + "_R.xml";
    cv::FileStorage fswrite(pathR, cv::FileStorage::WRITE);
    if(fswrite.isOpened())
    {
        fswrite << "_R" << matR ;
        fswrite.release();
    } 

    std::string pathT = common_patho + "_T.xml";
    fswrite.open(pathT, cv::FileStorage::WRITE);
    if(fswrite.isOpened())
    {
        fswrite <<"_T" << matT ;
        fswrite.release();
    }
    return 0; 

}


////////////////////////////////////////////////////////////////////////////////////////////////
int main_single( int argc, char* argv[])
{

  //[0]Load Parameters left
  cv::Mat cameraMatrixLeft(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixLeft);
  std::string common_patho = "../cameras_params/w/";
  std::string pathM1 = common_patho + "_M1.xml";

  cv::FileStorage fs(pathM1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M1"]>> cameraMatrixLeft;
      fs.release();
  } 
  std::cout << "cameraMatrixLeft: " << cameraMatrixLeft << std::endl;
 
  cv::Mat distCoeffsLeft;

  std::string pathD1 = common_patho + "_D1.xml";
  fs.open(pathD1, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D1"]>> distCoeffsLeft;
      fs.release();
  } 
  std::cout << "distCoeffsLeft: " << distCoeffsLeft << std::endl;

  //[1]Load Parameters right
  cv::Mat cameraMatrixRight(3,3,cv::DataType<double>::type);
  cv::setIdentity(cameraMatrixRight);

  std::string pathM2 = common_patho + "_M2.xml";

  fs.open(pathM2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_M2"]>> cameraMatrixRight;
      fs.release();
  } 
  std::cout << "cameraMatrixRight: " << cameraMatrixRight << std::endl;
 
  cv::Mat distCoeffsRight;

  std::string pathD2 = common_patho + "_D2.xml";
  fs.open(pathD2, cv::FileStorage::READ);
  if(fs.isOpened())
  {
      fs["_D2"]>> distCoeffsRight;
      fs.release();
  } 
  std::cout << "distCoeffsRight: " << distCoeffsRight << std::endl;

  //[2] Read points left
  IBD_SD::CShellCommand shellCommandObj; 
  
  std::string path    = "/home/ibd01/DVISION/src/tools/calibration/StereoCalibration/data/points/*.txt"; 
  std::string cmd_str = "ls " + path; 
 
  std::vector<std::string> veclists;
  shellCommandObj.getformatList(cmd_str, veclists);

  Problem problem;
  const int     numParameters  = 12; 
  double StereoCalibrateParameters[12] = {0};
  
  //std::vector< std::vector<double> > vvecStereoCalibrateParameters;
  //vvecStereoCalibrateParameters.resize(veclists.size()); 

for(int i = 0; i < veclists.size(); i++)
{

/*
  std::string path_uvleft  = "/home/ibd01/DVISION/src/tools/calibration/SolvePNP-master/data/imu2camera/imu2leftCamera/LoadUV.txt"; 
  std::string path_xyzleft = "/home/ibd01/DVISION/src/tools/calibration/SolvePNP-master/data/imu2camera/imu2leftCamera/LoadXYZ.txt"; 
  //[3] Read points right
  std::string path_uvright  = "/home/ibd01/DVISION/src/tools/calibration/SolvePNP-master/data/imu2camera/imu2rightCamera/LoadUV.txt"; 
  std::string path_xyzright = "/home/ibd01/DVISION/src/tools/calibration/SolvePNP-master/data/imu2camera/imu2rightCamera/LoadXYZ.txt"; 
*/
  //vvecStereoCalibrateParameters[i].resize(numParameters); 

  //[4]calc Rvec / Tvec 
  cv::Mat rvecleft(1,3,cv::DataType<double>::type);
  cv::Mat tvecleft(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsLeft;
  std::vector<cv::Point2d> imageLeft; 

  objPointsLeft.resize(0);
  imageLeft.resize(0);
 
  
  IBD::calcLeftRT(cameraMatrixLeft, distCoeffsLeft, veclists[i] ,objPointsLeft, rvecleft, tvecleft, imageLeft);
  //return 0;
  std::cout<<imageLeft[0].x << "," << imageLeft[0].y << "," << objPointsLeft[0].x << "," << objPointsLeft[0].y<< "," << objPointsLeft[0].z << std::endl; 

  cv::Mat rvecright(1,3,cv::DataType<double>::type);
  cv::Mat tvecright(3,1,cv::DataType<double>::type);
  std::vector<cv::Point3d> objPointsRight;
  std::vector<cv::Point2d> imageRight; 

  objPointsRight.resize(0);
  imageRight.resize(0);
 
  IBD::calcRightRT(cameraMatrixRight, distCoeffsRight, veclists[i] ,objPointsRight, rvecright, tvecright, imageRight);
  std::cout<<imageRight[0].x << "," << imageRight[0].y << "," << objPointsRight[0].x << "," << objPointsRight[0].y<< "," << objPointsRight[0].z << std::endl; 

  //return 0; 
  std::cout<<"objPointsLeft.size() : " << objPointsLeft.size() << std::endl; 
  std::cout<<"objPointsRight.size() : " << objPointsRight.size() << std::endl;
  assert(objPointsLeft.size() == objPointsRight.size()); 


  //return 0; 
  //[5]calc R & T 
/* 
  rvecleft.at<double>(0,0)  = 0.00276471; 
  rvecleft.at<double>(0,1)  = 0.0385392;
  rvecleft.at<double>(0,2)  = -0.00296829;
  rvecright.at<double>(0,0) = -0.00127287; 
  rvecright.at<double>(0,1) = 0.0407736;
  rvecright.at<double>(0,2) = -0.00598665;

  tvecleft.at<double>(0,0)  = 0.568852; 
  tvecleft.at<double>(1,0)  = -0.288769;
  tvecleft.at<double>(2,0)  = -0.0438256;
  tvecright.at<double>(0,0) = -0.532085; 
  tvecright.at<double>(1,0) = -0.280421;
  tvecright.at<double>(2,0) = -0.0285846;
*/
  Mat Rl, Tl, Rr,Tr, Rlinv, R, T;

  Rodrigues(rvecleft, Rl);
  Rodrigues(rvecright, Rr);

  Tr = tvecright; 
  Tl = tvecleft; 
/*
  Rodrigues(rvecleft, Rl);
  Rodrigues(rvecright, Rr);

  Tl = tvecleft; 
  Tr = tvecright; 
*/

  std::cout<<"Rodrigues over!" << std::endl; 
  //Rlinv = Rl; 
  //transpose(Rl, Rlinv); 
  Rlinv = Rl.inv(DECOMP_LU); 
  R = Rr * Rlinv; 
  T = Tr - R * Tl; 

  std::cout<<"R: " << R << "," << "T: " << T << std::endl; 

  //return 0; 
  std::string pathRR = common_patho + "_RR.xml";
  cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite << "_R" << R ;
      fswrite.release();
  } 

  std::string pathTT = common_patho + "_TT.xml";
  fswrite.open(pathTT, cv::FileStorage::WRITE);
  if(fswrite.isOpened())
  {
      fswrite <<"_T" << T ;
      fswrite.release();
  }
//return 0; 
  //[6]calc new points 
  std::vector<cv::Point3d> objPointsLeftCam;
  std::vector<cv::Point3d> objPointsRightCam;

  for(int i = 0 ; i < 1; i++) //i < 1; i++) //objPointsLeft.size(); i++)
  {
      cv::Point3d ptleftworld = objPointsLeft[i]; 
      Mat PtWorldl = (Mat_<double>(3,1)<<ptleftworld.x, ptleftworld.y, ptleftworld.z);
      
      PtWorldl     = Rl * PtWorldl + Tl;  
      ptleftworld.x= PtWorldl.at<double>(0,0); 
      ptleftworld.y= PtWorldl.at<double>(1,0); 
      ptleftworld.z= PtWorldl.at<double>(2,0); 

      cv::Point3d ptrightworld = objPointsRight[i]; 
      Mat PtWorldr = (Mat_<double>(3,1)<<ptrightworld.x, ptrightworld.y, ptrightworld.z);
      
      PtWorldr      = Rr * PtWorldr + Tr;   
      ptrightworld.x= PtWorldr.at<double>(0,0); 
      ptrightworld.y= PtWorldr.at<double>(1,0); 
      ptrightworld.z= PtWorldr.at<double>(2,0);

      //objPointsLeftCam.push_back(ptleftworld);
      //objPointsRightCam.push_back(ptrightworld); 

      std::cout<< ptleftworld.x << "," << ptleftworld.y << "," << ptleftworld.z <<"," << ptrightworld.x << "," << ptrightworld.y << "," << ptrightworld.z << std::endl; 

      cv::Point3d ptrightworldnew; 
      PtWorldl = R * PtWorldl + T; 
      ptrightworldnew.x = PtWorldl.at<double>(0,0);
      ptrightworldnew.y = PtWorldl.at<double>(1,0);
      ptrightworldnew.z = PtWorldl.at<double>(2,0);

      double diffx = (ptrightworldnew.x- ptrightworld.x); 
      double diffy = (ptrightworldnew.y- ptrightworld.y);
      double diffz = (ptrightworldnew.z- ptrightworld.z);
      double dist  = diffx * diffx + diffy * diffy + diffz * diffz; 
      dist  = sqrt(dist);
      std::cout<<ptrightworldnew.x << "," << ptrightworldnew.y << "," << ptrightworldnew.z << "," << ptrightworld.x << "," << ptrightworld.y <<"," << ptrightworld.z << "," << dist<< std::endl; 

/*      
  std::vector<cv::Point2d> projectedPoints; 
  cv::projectPoints(objPointsLeft, rvecleft, tvecleft, cameraMatrixLeft, distCoeffsLeft, projectedPoints);
 
  int rows = 2160; 
  int cols = 4096; 
  
  cv::Mat matShow(rows, cols, CV_8UC3, Scalar::all(100));
 
  for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  {
    std::cout << "Image point: " << imageLeft[i].x<<"," << imageLeft[i].y << " Projected to " << projectedPoints[i] << std::endl;
     circle(matShow, cv::Point(imageLeft[i].x, imageLeft[i].y), 15,  Scalar(255,0,0), 2, 8); 
     circle(matShow, projectedPoints[i], 15,  Scalar(0,0,255), 2, 8); 
  }
*/

  }
    //return 0; 
    //[7] ceres 
    //double StereoCalibrateParameters[12] = {0};

    StereoCalibrateParameters[0] =   rvecleft.at<double>(0,0); 
    StereoCalibrateParameters[1] =   rvecleft.at<double>(1,0);
    StereoCalibrateParameters[2] =   rvecleft.at<double>(2,0);
    StereoCalibrateParameters[3] =   rvecright.at<double>(0,0); 
    StereoCalibrateParameters[4] =   rvecright.at<double>(1,0);
    StereoCalibrateParameters[5] =   rvecright.at<double>(2,0);

    StereoCalibrateParameters[6]  = tvecleft.at<double>(0,0); 
    StereoCalibrateParameters[7]  = tvecleft.at<double>(1,0);
    StereoCalibrateParameters[8]  = tvecleft.at<double>(2,0);
    StereoCalibrateParameters[9]  = tvecright.at<double>(0,0); 
    StereoCalibrateParameters[10] = tvecright.at<double>(1,0);
    StereoCalibrateParameters[11] = tvecright.at<double>(2,0);

/*
    vvecStereoCalibrateParameters[i][0] =   rvecleft.at<double>(0,0); 
    vvecStereoCalibrateParameters[i][1] =   rvecleft.at<double>(1,0);
    vvecStereoCalibrateParameters[i][2] =   rvecleft.at<double>(2,0);
    vvecStereoCalibrateParameters[i][3] =   rvecright.at<double>(0,0); 
    vvecStereoCalibrateParameters[i][4] =   rvecright.at<double>(1,0);
    vvecStereoCalibrateParameters[i][5] =   rvecright.at<double>(2,0);

    vvecStereoCalibrateParameters[i][6]  = tvecleft.at<double>(0,0); 
    vvecStereoCalibrateParameters[i][7]  = tvecleft.at<double>(1,0);
    vvecStereoCalibrateParameters[i][8]  = tvecleft.at<double>(2,0);
    vvecStereoCalibrateParameters[i][9]  = tvecright.at<double>(0,0); 
    vvecStereoCalibrateParameters[i][10] = tvecright.at<double>(1,0);
    vvecStereoCalibrateParameters[i][11] = tvecright.at<double>(2,0);
*/

    for(int i = 0; i < 1; i++)//objPointsLeftCam.size() ; i++)//100; i++) //objPointsLeftCam.size() ; i++)
    {
        //std::cout<<imageLeft[i].x <<"," << imageLeft[i].y << "," <<imageRight[i].x <<"," << imageRight[i].y << std::endl;   
        //return 0; 
        StereoCalibrationResidual* pResidualX = new StereoCalibrationResidual(
                                                    imageLeft[i].x, 
                                                    imageLeft[i].y,
                                                    objPointsLeft[i].x, 
                                                    objPointsLeft[i].y, 
                                                    objPointsLeft[i].z, 
                                                    imageRight[i].x, 
                                                    imageRight[i].y, 
                                                    objPointsRight[i].x, 
                                                    objPointsRight[i].y, 
                                                    objPointsRight[i].z, 
                                                    cameraMatrixLeft.at<double>(0,0),
                                                    cameraMatrixLeft.at<double>(0,2),
                                                    cameraMatrixLeft.at<double>(1,1), 
                                                    cameraMatrixLeft.at<double>(1,2),
                                                    distCoeffsLeft.at<double>(0,0), 
                                                    distCoeffsLeft.at<double>(0,1), 
                                                    distCoeffsLeft.at<double>(0,2), 
                                                    distCoeffsLeft.at<double>(0,3),
                                                    distCoeffsLeft.at<double>(0,4), 
                                                    cameraMatrixRight.at<double>(0,0),
                                                    cameraMatrixRight.at<double>(0,2),
                                                    cameraMatrixRight.at<double>(1,1),
                                                    cameraMatrixRight.at<double>(1,2),
                                                    distCoeffsRight.at<double>(0,0),
                                                    distCoeffsRight.at<double>(0,1),
                                                    distCoeffsRight.at<double>(0,2),
                                                    distCoeffsRight.at<double>(0,3),
                                                    distCoeffsRight.at<double>(0,4));

        CostFunction* cost_function = new NumericDiffCostFunction<StereoCalibrationResidual,ceres::CENTRAL, 1, 28>(pResidualX);
        problem.AddResidualBlock(cost_function, new CauchyLoss(5), StereoCalibrateParameters);
    }
}
    Solver::Options  options;
    options.max_num_iterations = 1000; 
    options.linear_solver_type = ceres::DENSE_QR; 
    options.minimizer_progress_to_stdout = true; 
 
    Solver::Summary summary; 
    Solve(options, &problem, &summary);

    StereoCalibrationResidual stereoCalibrationObj; 
    stereoCalibrationObj.GetCalibrationParameters(StereoCalibrateParameters);

    return 0; 
/*
    fstream file("/home/ibd01/DVISION/src/Binocular_Ceres_Optimize/data/camera.txt",ios::in);
    if(!file.is_open())
    {
         std::cout<< "open error !"<< std::endl; 
         return 0; 
    }

    int      nCount = 4; 
    double * pData  = new double[5 * nCount];

    string       strLine;  

    int      count  = 0; 
    while(1)
    {
        if(!getline(file,strLine,','))
        {
            std::cout<<"finished !" << std::endl; 
            break; 
        }
        //std::cout<<strLine << std::endl; 

        stringstream ss;
        ss << strLine; 
        ss >> pData[count];
        
        ss.clear();
        std::cout << pData[count] << std::endl;  
        count = count + 1;
        if(count == 20)
        {
            break; 
        }        
    }

    double dBackCrossParameters[6] = {0};
    int    mCount  = 5; 
    for(int i = 0 ; i < nCount ; i = i + 1)
    {
        std::cout<<pData[mCount*i]<<","<<pData[mCount*i + 1]<<","<<pData[mCount*i + 2]<<","<<pData[mCount*i + 3]<<","<<pData[mCount*i + 4]<<std::endl; 
        dBackCrossParameters[0] += pData[mCount*i + 2]; 
        dBackCrossParameters[1] += pData[mCount*i + 3];
    }//data read

    file.close(); 
   
    double df               = 153.24; //mm
    dBackCrossParameters[0] = dBackCrossParameters[0]/4;
    dBackCrossParameters[1] = dBackCrossParameters[1]/4; 
    std::cout<<dBackCrossParameters[0]<<","<<dBackCrossParameters[1]<<std::endl;

    dBackCrossParameters[2] = 50 * df ; //Z
    
    Problem problem; 
    for(int i = 0; i < nCount ; i++)
    {
        double* nPoint = pData + 5 * i; 
        
        BackCrossResidual* pResidualX = new BackCrossResidual(nPoint[0]/1000, nPoint[1]/1000, df/1000, nPoint[2], nPoint[3], nPoint[4]);

        CostFunction* cost_function = new AutoDiffCostFunction<BackCrossResidual, 2, 6>(pResidualX);
        problem.AddResidualBlock(cost_function, new CauchyLoss(5), dBackCrossParameters);
    }

    Solver::Options  options;
    options.max_num_iterations = 50; 
    options.linear_solver_type = ceres::DENSE_QR; 
    options.minimizer_progress_to_stdout = true; 
 
    Solver::Summary summary; 
    Solve(options, &problem, &summary);

    std::cout<< summary.BriefReport() << std::endl; 
    std::cout<<  dBackCrossParameters[0] << "," << dBackCrossParameters[1] << "," << dBackCrossParameters[2] << std::endl; 
    std::cout<<  dBackCrossParameters[3] << "," << dBackCrossParameters[4] << "," << dBackCrossParameters[5] << std::endl;   

    BackCrossResidual pResidualX; 
    pResidualX.RotationTrans(dBackCrossParameters);

*/

}

