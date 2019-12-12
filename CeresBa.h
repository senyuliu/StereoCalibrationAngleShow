#ifndef __CERESBA_H__
#define __CERESBA_H__

#include <iostream>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv2/calib3d/calib3d.hpp"

#include <vector>
#include "ceres/ceres.h"
#include "glog/logging.h"

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>
#include <opencv2/core/eigen.hpp>

using namespace cv; 

using namespace ceres; 
//using namespace Eigen;

struct PointDBLXYZ {
	double x;
	double y;
	double z;
};
struct PointDBLXYZT {
	double x;
	double y;
	double z;
        double time; 
};
struct PointDBLXYZK {
	double x;
	double y;
	double z;
        double k; 
        double c; 
        double cost; 
        double time; 
};

class CeresBa
{
public: 
    CeresBa(void);
    ~CeresBa(void);

public: 
    PointDBLXYZ ShiftingCoordinateValue(std::vector<PointDBLXYZ> & pointClouds,
	pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud);
    PointDBLXYZT ShiftingCoordinateValueT(std::vector<PointDBLXYZT> & pointClouds,
	std::vector<PointDBLXYZT> & cloud);
public: 
    //offset variables
    double maxX;
    double maxY;
    double maxZ;
    PointDBLXYZT               slaneData; 
    std::vector<PointDBLXYZT>  glaneData;
}; 

struct CostFunctor
{
    template<typename T>
    bool operator() (const T* const x, T* residual) const
    {
        residual[0] = T(5.0) - x[0]; 
        return true; 
    }
};

void  TwoMatrixsMultiplyLane(double returnMatrix[3][3],const double A[3][3],const double B[3][3]);

void matrixMultiVectorLane(double coordinate[3], double R[3][3],double T[3]);

void matrixMultiVector(double coordinate[3], double R[3][3],double T[3]);

void matrixInv(double matrixin[3][3], double matrixinv[3][3]);

/*----------------------------exp(mx+c)---------------------
struct ExponentialResidual
{
    ExponentialResidual(double x, double y) : x_(x), y_(y)
    {

    }
    ExponentialResidual()
    {

    }
    template<typename T>
    bool operator() (const T* const m, const T* const c, T* residual) const 
    {
        //residual[0] = T(y_) - exp( m[0] * T(x_) + c[0] );  
        //residual[0] = T(y_) - exp( ( m[0] * T(x_) ) + c[0] );   
        //residual[0] = T(y_) - ( m[0]*T(x_) + c[0] )  ;  //exp( ( m[0] * T(x_) ) + c[0] );  
        
        return true; 
    }

    template<typename T>
    double getPoints(const T* m, const T* c, T x)
    {
       //double y = exp( m[0] * x + c[0] );
       //double y  = m[0] * x + c[0];
       return y; 
    }
    
    private: //observ point 
    //const double x_; 
    //const double y_;
    double x_; 
    double y_;  
};
//--------------------------------------------------------*/

/*----------------Aexp(mx + c)----------------------------
struct ExponentialResidual
{
    ExponentialResidual(double x, double y) : x_(x), y_(y)
    {

    }
    ExponentialResidual()
    {

    }
    template<typename T>
    bool operator() (const T* const m, const T* const c, const T* a, T* residual) const 
    {
        //residual[0] = T(y_) - exp( m[0] * T(x_) + c[0] );  
        //residual[0] = T(y_) - exp( ( m[0] * T(x_) ) + c[0] );   
        //residual[0] = T(y_) - ( m[0]*T(x_) + c[0] )  ;  //exp( ( m[0] * T(x_) ) + c[0] );  
        residual[0]   = T(y_) - ( a[0] + exp( ( (-m[0]) * T(x_) ) + c[0] ) );
        return true; 
    }

    template<typename T>
    double getPoints(const T* m, const T* c, const T* a, T x)
    {
       //double y = exp( m[0] * x + c[0] );
       //double y  = m[0] * x + c[0];
       std::cout<<"param: " << ((-m[0]) * x ) + c[0] << std::endl; 
       double y = a[0] + exp( ( (-m[0]) * x ) + c[0] );
       return y; 
    }
    
    private: //observ point 
    //const double x_; 
    //const double y_;
    double x_; 
    double y_;  
};
//-------------------------------------------------------*/

/*----------------Aexp(mx + c)----------------------------
struct ExponentialResidual
{
    ExponentialResidual(double x, double y) : x_(x), y_(y)
    {

    }
    ExponentialResidual()
    {

    }
    template<typename T>
    bool operator() (const T* const m, const T* const c, const T* a,const T* const m1, const T* const c1, const T* a1 ,T* residual) const 
    {
        //residual[0] = T(y_) - exp( m[0] * T(x_) + c[0] );  
        //residual[0] = T(y_) - exp( ( m[0] * T(x_) ) + c[0] );   
        //residual[0] = T(y_) - ( m[0]*T(x_) + c[0] )  ;  //exp( ( m[0] * T(x_) ) + c[0] );  
        //residual[0]   = T(y_) - ( a[0] + exp( ( (m[0]) * T(x_) ) + c[0] ) ) - ( a1[0] + exp( ( (m1[0]) * T(x_) ) + c1[0] ) );
        residual[0]   = T(y_) - ( a[0] * exp( m[0] * T(x_) + c[0] ) + m1[0] );
        //residual[0]  =  T(y_) - ( m[0]*(T(x_)*T(x_)*T(x_)) + c[0]* (T(x_)*T(x_)) + a[0] * T(x_) + m1[0]);
        return true; 
    }

    template<typename T>
    double getPoints(const T* m, const T* c, const T* a, const T* const m1, const T* const c1, const T* a1 ,T x)
    {
       //double y = exp( m[0] * x + c[0] );
       //double y  = m[0] * x + c[0];
       double y =  a[0] * exp( m[0] * x + c[0] ) + m1[0]; //+ a[0] ; 
       //double y = ( a[0] + exp( ( (m[0]) * T(x_) ) + c[0] ) ) - ( a1[0] + exp( ( (m1[0]) * T(x_) ) + c1[0] ) );
       //double y = m[0]*(T(x_)*T(x_)*T(x_)) + c[0]* (T(x_)*T(x_)) + a[0] * T(x_) + m1[0]; 
       return y; 
    }
    
    private: //observ point 
    //const double x_; 
    //const double y_;
    double x_; 
    double y_;  
};
//-------------------------------------------------------*/

///*----------------Aexp(mx + c)----------------------------
struct ExponentialResidual
{
    ExponentialResidual(double x, double y) : x_(x), y_(y)
    {

    }
    ExponentialResidual()
    {

    }
    template<typename T>
    bool operator() (const T* const m, const T* const c, T* residual) const 
    {
        residual[0]   = T(y_) - ( m[0] * T(x_) + c[0] ) ;

        return true; 
    }

    template<typename T>
    double getPoints(const T* m, const T* c ,T x)
    {

       double y =   m[0] * x + c[0] ; //+ a[0] ; 
       return y; 
    }
    
    private: //observ point 
    //const double x_; 
    //const double y_;
    double x_; 
    double y_;  
};
//-------------------------------------------------------*/

/*------
StereoCalibrationResidual(double tlX, double tlY, double tlZ, double trX, double trY, double trZ, double rlx, double rly, double rlz, double rrx, double rry, double rrz)
        :_tlx(tlX), _tly(tlY), _tlz(tlZ), _trx(trX), _try(trY), _trz(trZ), _rlx(rlx), _rly(rly), _rlz(rlz), _rrx(rrx), _rry(rry), _rrz(rrz)
cv::Mat cameraMatrixl, cv::Mat distCoeffsl, cv::Mat cameraMatrixr, cv::Mat distCoeffsr
---------*/
struct StereoCalibrationResiduals
{
    StereoCalibrationResiduals(
                              double imageu, 
                              double imagev, 
                              double objectx, 
                              double objecty, 
                              double objectz, 
                              double imageur, 
                              double imagevr, 
                              double objectxr, 
                              double objectyr, 
                              double objectzr, 
                              double m00, 
                              double m02, 
                              double m11,
                              double m12, 
                              double d00, 
                              double d10, 
                              double d20, 
                              double d30, 
                              double d40, 
                              double mr00, 
                              double mr02,
                              double mr11,  
                              double mr12, 
                              double dr00, 
                              double dr10, 
                              double dr20, 
                              double dr30, 
                              double dr40  ): _imageu(imageu), _imagev(imagev), _objectx(objectx), _objecty(objecty), _objectz(objectz),_imageur(imageur), _imagevr(imagevr), _objectxr(objectxr), _objectyr(objectyr), _objectzr(objectzr), _m00(m00), _m02(m02), _m11(m11),_m12(m12), _d00(d00), _d10(d10), _d20(d20), _d30(d30),_d40(d40), _mr00(mr00), _mr02(mr02), _mr11(mr11), _mr12(mr12), _dr00(dr00), _dr10(dr10),_dr20(dr20),_dr30(dr30),_dr40(dr40)
    {
    
    }

    StereoCalibrationResiduals()
    {

    }
    template<typename T>
    bool operator() (const T* const Rotation, const T* const Translation, const T* const Rotationl, const T* const Translationl,T* residual) const 
    {

        //[0]get data input  
        T pointsl[2] = {T(_imageu), T(_imagev)}; 
        T pointsr[2] = {T(_imageur), T(_imagevr)}; 
        T objectl[3] = {T(_objectx), T(_objecty), T(_objectz)}; 
        T objectr[3] = {T(_objectxr),T(_objectyr),T(_objectzr)}; 

        //[1]get parameters 
        T cameral[4]  = {T(_m00),  T(_m02), T(_m11), T(_m12)}; 
        T camerar[4]  = {T(_mr00), T(_mr02),T(_mr11),T(_mr12)};
        T distCoffl[5]= {T(_d00), T(_d10), T(_d20), T(_d30), T(_d40)};
        T distCoffr[5]= {T(_dr00), T(_dr10), T(_dr20), T(_dr30), T(_dr40)};

        //[2]rotaion R
        T theta[3]; 
        theta[0] = Rotation[0];
        theta[1] = Rotation[1];
        theta[2] = Rotation[2]; 

        T R_x[3][3] = {T(1), T(0),                T(0), 
                       T(0), T( cos(theta[0]) ) , T( -sin(theta[0]) ),
                       T(0), T( sin(theta[0]) ) , T(  cos(theta[0]) ) };

        T R_y[3][3] = {T( cos(theta[1]) ), T( 0 ), T( sin(theta[1]) ), 
                       T(0),               T( 1 ) , T(0 ),
                       T( -sin(theta[1]) ),T( 0 ) , T(  cos(theta[1]) ) };

        T R_z[3][3] = {T( cos(theta[2]) ), T( -sin(theta[2]) ), T( 0 ), 
                       T( sin(theta[2]) ), T(  cos(theta[2]) ) , T( 0 ),
                       T(0),               T( 0 ) ,              T( 1 )};     

        T Temp[3][3]; 

        Temp[0][0] = R_z[0][0] * R_y[0][0] + R_z[0][1] * R_y[1][0] + R_z[0][2] * R_y[2][0];
        Temp[0][1] = R_z[0][0] * R_y[0][1] + R_z[0][1] * R_y[1][1] + R_z[0][2] * R_y[2][1];
        Temp[0][2] = R_z[0][0] * R_y[0][2] + R_z[0][1] * R_y[1][2] + R_z[0][2] * R_y[2][2];

        Temp[1][0] = R_z[1][0] * R_y[0][0] + R_z[1][1] * R_y[1][0] + R_z[1][2] * R_y[2][0];
        Temp[1][1] = R_z[1][0] * R_y[0][1] + R_z[1][1] * R_y[1][1] + R_z[1][2] * R_y[2][1];
        Temp[1][2] = R_z[1][0] * R_y[0][2] + R_z[1][1] * R_y[1][2] + R_z[1][2] * R_y[2][2];

        Temp[2][0] = R_z[2][0] * R_y[0][0] + R_z[2][1] * R_y[1][0] + R_z[2][2] * R_y[2][0];
        Temp[2][1] = R_z[2][0] * R_y[0][1] + R_z[2][1] * R_y[1][1] + R_z[2][2] * R_y[2][1];
        Temp[2][2] = R_z[2][0] * R_y[0][2] + R_z[2][1] * R_y[1][2] + R_z[2][2] * R_y[2][2];

        T A[3][3];

        A[0][0] = Temp[0][0] * R_x[0][0] + Temp[0][1] * R_x[1][0] + Temp[0][2] * R_x[2][0];
        A[0][1] = Temp[0][0] * R_x[0][1] + Temp[0][1] * R_x[1][1] + Temp[0][2] * R_x[2][1];
        A[0][2] = Temp[0][0] * R_x[0][2] + Temp[0][1] * R_x[1][2] + Temp[0][2] * R_x[2][2];

        A[1][0] = Temp[1][0] * R_x[0][0] + Temp[1][1] * R_x[1][0] + Temp[1][2] * R_x[2][0];
        A[1][1] = Temp[1][0] * R_x[0][1] + Temp[1][1] * R_x[1][1] + Temp[1][2] * R_x[2][1];
        A[1][2] = Temp[1][0] * R_x[0][2] + Temp[1][1] * R_x[1][2] + Temp[1][2] * R_x[2][2];

        A[2][0] = Temp[2][0] * R_x[0][0] + Temp[2][1] * R_x[1][0] + Temp[2][2] * R_x[2][0];
        A[2][1] = Temp[2][0] * R_x[0][1] + Temp[2][1] * R_x[1][1] + Temp[2][2] * R_x[2][1];
        A[2][2] = Temp[2][0] * R_x[0][2] + Temp[2][1] * R_x[1][2] + Temp[2][2] * R_x[2][2];

        //[2]rotaion Rl 
        T thetal[3]; 
        thetal[0] = Rotationl[0];
        thetal[1] = Rotationl[1];
        thetal[2] = Rotationl[2]; 

        T R_xl[3][3] = {T(1), T(0),                T(0), 
                       T(0), T( cos(thetal[0]) ) , T( -sin(thetal[0]) ),
                       T(0), T( sin(thetal[0]) ) , T(  cos(thetal[0]) ) };

        T R_yl[3][3] = {T( cos(thetal[1]) ), T( 0 ), T( sin(thetal[1]) ), 
                       T(0),               T( 1 ) , T(0 ),
                       T( -sin(thetal[1]) ),T( 0 ) , T(  cos(thetal[1]) ) };

        T R_zl[3][3] = {T( cos(thetal[2]) ), T( -sin(thetal[2]) ), T( 0 ), 
                       T( sin(thetal[2]) ), T(  cos(thetal[2]) ) , T( 0 ),
                       T(0),               T( 0 ) ,              T( 1 )};     

        T Templ[3][3]; 

        Templ[0][0] = R_zl[0][0] * R_yl[0][0] + R_zl[0][1] * R_yl[1][0] + R_zl[0][2] * R_yl[2][0];
        Templ[0][1] = R_zl[0][0] * R_yl[0][1] + R_zl[0][1] * R_yl[1][1] + R_zl[0][2] * R_yl[2][1];
        Templ[0][2] = R_zl[0][0] * R_yl[0][2] + R_zl[0][1] * R_yl[1][2] + R_zl[0][2] * R_yl[2][2];

        Templ[1][0] = R_zl[1][0] * R_yl[0][0] + R_zl[1][1] * R_yl[1][0] + R_zl[1][2] * R_yl[2][0];
        Templ[1][1] = R_zl[1][0] * R_yl[0][1] + R_zl[1][1] * R_yl[1][1] + R_zl[1][2] * R_yl[2][1];
        Templ[1][2] = R_zl[1][0] * R_yl[0][2] + R_zl[1][1] * R_yl[1][2] + R_zl[1][2] * R_yl[2][2];

        Templ[2][0] = R_zl[2][0] * R_yl[0][0] + R_zl[2][1] * R_yl[1][0] + R_zl[2][2] * R_yl[2][0];
        Templ[2][1] = R_zl[2][0] * R_yl[0][1] + R_zl[2][1] * R_yl[1][1] + R_zl[2][2] * R_yl[2][1];
        Templ[2][2] = R_zl[2][0] * R_yl[0][2] + R_zl[2][1] * R_yl[1][2] + R_zl[2][2] * R_yl[2][2];

        T B[3][3];

        B[0][0] = Templ[0][0] * R_xl[0][0] + Templ[0][1] * R_xl[1][0] + Templ[0][2] * R_xl[2][0];
        B[0][1] = Templ[0][0] * R_xl[0][1] + Templ[0][1] * R_xl[1][1] + Templ[0][2] * R_xl[2][1];
        B[0][2] = Templ[0][0] * R_xl[0][2] + Templ[0][1] * R_xl[1][2] + Templ[0][2] * R_xl[2][2];

        B[1][0] = Templ[1][0] * R_xl[0][0] + Templ[1][1] * R_xl[1][0] + Templ[1][2] * R_xl[2][0];
        B[1][1] = Templ[1][0] * R_xl[0][1] + Templ[1][1] * R_xl[1][1] + Templ[1][2] * R_xl[2][1];
        B[1][2] = Templ[1][0] * R_xl[0][2] + Templ[1][1] * R_xl[1][2] + Templ[1][2] * R_xl[2][2];

        B[2][0] = Templ[2][0] * R_xl[0][0] + Templ[2][1] * R_xl[1][0] + Templ[2][2] * R_xl[2][0];
        B[2][1] = Templ[2][0] * R_xl[0][1] + Templ[2][1] * R_xl[1][1] + Templ[2][2] * R_xl[2][1];
        B[2][2] = Templ[2][0] * R_xl[0][2] + Templ[2][1] * R_xl[1][2] + Templ[2][2] * R_xl[2][2];
        
        T Rr[3][3];

        Rr[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
        Rr[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
        Rr[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

        Rr[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
        Rr[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
        Rr[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

        Rr[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
        Rr[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
        Rr[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
        

        //[3]translation Tr
        T Tr[3]; 
        Tr[0] = A[0][0] * Translationl[0] + A[0][1] * Translationl[1] + A[0][2] * Translationl[2];
        Tr[1] = A[1][0] * Translationl[0] + A[1][1] * Translationl[1] + A[1][2] * Translationl[2];
        Tr[2] = A[2][0] * Translationl[0] + A[2][1] * Translationl[1] + A[2][2] * Translationl[2];

        Tr[0] += Translation[0];
        Tr[1] += Translation[1];
        Tr[2] += Translation[2];

        //[4]projectpoints left 
/*
        std::cout<<"Translationlx: " << Translationl[0] << std::endl;
        std::cout<<"Translationly: " << Translationl[1] << std::endl; 
        std::cout<<"Translationlz: " << Translationl[2] << std::endl;  
        std::cout<<"B00: " << B[0][0] << std::endl; 
        std::cout<<"B01: " << B[0][1] << std::endl; 
        std::cout<<"B02: " << B[0][2] << std::endl; 
        std::cout<<"B10: " << B[1][0] << std::endl; 
        std::cout<<"B11: " << B[1][1] << std::endl; 
        std::cout<<"B12: " << B[1][2] << std::endl; 
        std::cout<<"B20: " << B[2][0] << std::endl; 
        std::cout<<"B21: " << B[2][1] << std::endl; 
        std::cout<<"B22: " << B[2][2] << std::endl; 
*/
        T objectsl[3]; 

        objectsl[0] = B[0][0] * objectl[0] + B[0][1] * objectl[1] + B[0][2] * objectl[2] ;
        objectsl[1] = B[1][0] * objectl[0] + B[1][1] * objectl[1] + B[1][2] * objectl[2] ;
        objectsl[2] = B[2][0] * objectl[0] + B[2][1] * objectl[1] + B[2][2] * objectl[2] ;

        objectsl[0] = objectsl[0] + Translationl[0];
        objectsl[1] = objectsl[1] + Translationl[1];
        objectsl[2] = objectsl[2] + Translationl[2];

        T xl = objectsl[0]/objectsl[2];
        T yl = objectsl[1]/objectsl[2]; 

        T x2 = xl*xl, y2 = yl*yl;
        T r2 = x2 + y2, _2xy = T(2)*xl*yl;

        T kl  = (T(1) + distCoffl[4]*r2*r2*r2 + distCoffl[1]*r2*r2 + distCoffl[0]*r2 ); 
        T Plx = distCoffl[2]*_2xy + distCoffl[3]*(r2 + T(2)*x2);
        T ply = distCoffl[2]*(r2 + T(2) * y2) + distCoffl[3] * _2xy; 

        T ul = (xl * kl + Plx) * cameral[0] + cameral[1]; 
        T vl = (yl * kl + ply) * cameral[2] + cameral[3];
/*
        std::cout<<"ul: " << ul << std::endl; 
        std::cout<<"vl: " << vl << std::endl; 
        std::cout<<"_imageu: " << _imageu << std::endl; 
        std::cout<<"_imagev: " << _imagev << std::endl; 
        std::cout<<"objectlx: " << objectl[0] << std::endl; 
        std::cout<<"objectly: " << objectl[1] << std::endl; 
        std::cout<<"objectlz: " << objectl[2] << std::endl; 
*/

        //[5]projectpoinst right 
        T objectsr[3];
        objectsr[0] = Rr[0][0] * objectr[0] + Rr[0][1] * objectr[1] + Rr[0][2] * objectr[2] ;
        objectsr[1] = Rr[1][0] * objectr[0] + Rr[1][1] * objectr[1] + Rr[1][2] * objectr[2] ;
        objectsr[2] = Rr[2][0] * objectr[0] + Rr[2][1] * objectr[1] + Rr[2][2] * objectr[2] ;

        objectsr[0] = objectsr[0] + Tr[0];
        objectsr[1] = objectsr[1] + Tr[1];
        objectsr[2] = objectsr[2] + Tr[2];

        T xr = objectsr[0]/objectsr[2];
        T yr = objectsr[1]/objectsr[2]; 

        T xr2 = xr*xr, yr2 = yr*yr;
        T rr2 = xr2 + yr2, _2xyr = T(2)*xr*yr;

        T kr  = (T(1) + distCoffr[4]*rr2*rr2*rr2 + distCoffr[1]*rr2*rr2 + distCoffr[0]*rr2 ); 
        T Prx = distCoffr[2]*_2xyr + distCoffr[3]*(rr2 + T(2)*xr2);
        T pry = distCoffr[2]*(rr2 + T(2) * yr2) + distCoffr[3] * _2xyr; 

        T ur = (xr * kr + Prx) * camerar[0] + camerar[1]; 
        T vr = (yr * kr + pry) * camerar[2] + camerar[3];
/*
        std::cout<<"ur: " << ur << std::endl; 
        std::cout<<"vr: " << vr << std::endl; 
        std::cout<<"_imageur: " << _imageur << std::endl; 
        std::cout<<"_imagevr: " << _imagevr << std::endl; 
        std::cout<<"objectrx: " << objectr[0] << std::endl; 
        std::cout<<"objectry: " << objectr[1] << std::endl; 
        std::cout<<"objectrz: " << objectr[2] << std::endl; 
*/ 
        //return false; 
        //[6]result 
        residual[0] = T( sqrt ( (ul-pointsl[0]) * (ul-pointsl[0]) + (vl - pointsl[1])*(vl - pointsl[1]) + (ur - pointsr[0]) * (ur - pointsr[0]) + (vr - pointsr[1]) * (vr - pointsr[1]) ) );
        //std::cout<<"residual[0]: " << residual[0] << std::endl; 
        return true; 
    }

public: 
    bool GetCalibrationParameters(const double* const StereoCalibrateParameters)
    {
        double rlx =  StereoCalibrateParameters[0];
        double rly =  StereoCalibrateParameters[1];
        double rlz =  StereoCalibrateParameters[2];
        double tlx =  StereoCalibrateParameters[3];
        double tly =  StereoCalibrateParameters[4];
        double tlz =  StereoCalibrateParameters[5];
        double rrx =  StereoCalibrateParameters[6];
        double rry =  StereoCalibrateParameters[7];
        double rrz =  StereoCalibrateParameters[8];
        double trx =  StereoCalibrateParameters[9];
        double trys = StereoCalibrateParameters[10];
        double trz =  StereoCalibrateParameters[11];

        cv::Mat rvecl = (cv::Mat_<double>(1,3)<< rlx, rly, rlz);
        cv::Mat tvecl = (cv::Mat_<double>(3,1)<< rrx, rry, rrz); 
        cv::Mat rvecr = (cv::Mat_<double>(1,3)<< tlx, tly, tlz);
        cv::Mat tvecr = (cv::Mat_<double>(3,1)<< trx, trys, trz);     
    
        std::cout<< "rvecl:" << rvecl << std::endl; 
        std::cout<< "tvecl:" << tvecl << std::endl;
        std::cout<< "rvecr:" << rvecr << std::endl; 
        std::cout<< "tvecr:" << tvecr << std::endl;           
       
        cv::Mat Rl, Rr; 
        Rodrigues(rvecl, Rl);
        Rodrigues(rvecr, Rr);

        cv::Mat Tl, Tr; 
        Tl = tvecl; 
        Tr = tvecr; 

        cv::Mat R, T; 
        Mat Rl_inv = Rl.inv(DECOMP_LU); 
        std::cout<<"Rr: " <<  Rr << std::endl; 
        std::cout<<"Rl: " <<  Rl << std::endl; 
        R = Rr * Rl_inv; 
        T = Tr - R * Tl; 

        std::cout<<"R: " << R << std::endl;
        std::cout<< "T: " << T << std::endl;

        std::string common_patho = "../cameras_params/w/";
        std::string pathRR = common_patho + "_R.xml";
        cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite << "_R" << R ;
            fswrite.release();
        } 
      
        std::string pathTT = common_patho + "_T.xml";
        fswrite.open(pathTT, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite <<"_T" << T ;
            fswrite.release();
        }
    }

private:
    double _imageu; 
    double _imagev; 
    double _objectx; 
    double _objecty; 
    double _objectz; 
    double _imageur; 
    double _imagevr; 
    double _objectxr; 
    double _objectyr; 
    double _objectzr;
    double _m00; 
    double _m02;
    double _m11;
    double _m12; 
    double _mr00;
    double _mr02;
    double _mr11;
    double _mr12; 
    double _d00;
    double _d10;
    double _d20;
    double _d30;
    double _d40;
    double _dr00;
    double _dr10;
    double _dr20;
    double _dr30;
    double _dr40;

    //cv::Mat _cameraMatrixl; 
    //cv::Mat _distCoeffsl; 
    //cv::Mat _cameraMatrixr; 
    //cv::Mat _distCoeffsr; 
    
};

struct StereoCalibrationResidualsr
{
    StereoCalibrationResidualsr(
                              double imageu, 
                              double imagev, 
                              double objectx, 
                              double objecty, 
                              double objectz, 
                              double imageur, 
                              double imagevr, 
                              double objectxr, 
                              double objectyr, 
                              double objectzr, 
                              double m00, 
                              double m02, 
                              double m11,
                              double m12, 
                              double d00, 
                              double d10, 
                              double d20, 
                              double d30, 
                              double d40, 
                              double mr00, 
                              double mr02,
                              double mr11,  
                              double mr12, 
                              double dr00, 
                              double dr10, 
                              double dr20, 
                              double dr30, 
                              double dr40  ): _imageu(imageu), _imagev(imagev), _objectx(objectx), _objecty(objecty), _objectz(objectz),_imageur(imageur), _imagevr(imagevr), _objectxr(objectxr), _objectyr(objectyr), _objectzr(objectzr), _m00(m00), _m02(m02), _m11(m11),_m12(m12), _d00(d00), _d10(d10), _d20(d20), _d30(d30),_d40(d40), _mr00(mr00), _mr02(mr02), _mr11(mr11), _mr12(mr12), _dr00(dr00), _dr10(dr10),_dr20(dr20),_dr30(dr30),_dr40(dr40)
    {
    
    }

    StereoCalibrationResidualsr()
    {

    }
    template<typename T>
    bool operator() (const T* const Rotation, const T* const Translation, const T* const Rotationl, const T* const Translationl,T* residual) const 
    {

        //[0]get data input  
        T pointsl[2] = {T(_imageu), T(_imagev)}; 
        T pointsr[2] = {T(_imageur), T(_imagevr)}; 
        T objectl[3] = {T(_objectx), T(_objecty), T(_objectz)}; 
        T objectr[3] = {T(_objectxr),T(_objectyr),T(_objectzr)}; 

        //[1]get parameters 
        T cameral[4]  = {T(_m00),  T(_m02), T(_m11), T(_m12)}; 
        T camerar[4]  = {T(_mr00), T(_mr02),T(_mr11),T(_mr12)};
        T distCoffl[5]= {T(_d00), T(_d10), T(_d20), T(_d30), T(_d40)};
        T distCoffr[5]= {T(_dr00), T(_dr10), T(_dr20), T(_dr30), T(_dr40)};

        //[2]rotaion Rr
        T A[3][3];
        T B[3][3];
        T Rr[3][3];

        for(int i = 0 ; i < 3 ; i++)
        {
            for(int j = 0 ; j < 3; j++)
            {
                A[i][j] = Rotation[i * 3 + j];
            }
        }
        for(int i = 0 ; i < 3 ; i++)
        {
            for(int j = 0 ; j < 3; j++)
            {
                B[i][j] = Rotationl[i * 3 + j];
            }
        }
        
        Rr[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
        Rr[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
        Rr[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

        Rr[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
        Rr[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
        Rr[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

        Rr[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
        Rr[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
        Rr[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
        

        //[3]translation Tr
        T Tr[3]; 
        Tr[0] = A[0][0] * Translationl[0] + A[0][1] * Translationl[1] + A[0][2] * Translationl[2];
        Tr[1] = A[1][0] * Translationl[0] + A[1][1] * Translationl[1] + A[1][2] * Translationl[2];
        Tr[2] = A[2][0] * Translationl[0] + A[2][1] * Translationl[1] + A[2][2] * Translationl[2];

        Tr[0] += Translation[0];
        Tr[1] += Translation[1];
        Tr[2] += Translation[2];

        //[4]projectpoints left 
        T objectsl[3]; 

        objectsl[0] = B[0][0] * objectl[0] + B[0][1] * objectl[1] + B[0][2] * objectl[2] ;
        objectsl[1] = B[1][0] * objectl[0] + B[1][1] * objectl[1] + B[1][2] * objectl[2] ;
        objectsl[2] = B[2][0] * objectl[0] + B[2][1] * objectl[1] + B[2][2] * objectl[2] ;

        objectsl[0] = objectsl[0] + Translationl[0];
        objectsl[1] = objectsl[1] + Translationl[1];
        objectsl[2] = objectsl[2] + Translationl[2];

        T xl = objectsl[0]/objectsl[2];
        T yl = objectsl[1]/objectsl[2]; 

        T x2 = xl*xl, y2 = yl*yl;
        T r2 = x2 + y2, _2xy = T(2)*xl*yl;

        T kl  = (T(1) + distCoffl[4]*r2*r2*r2 + distCoffl[1]*r2*r2 + distCoffl[0]*r2 ); 
        T Plx = distCoffl[2]*_2xy + distCoffl[3]*(r2 + T(2)*x2);
        T ply = distCoffl[2]*(r2 + T(2) * y2) + distCoffl[3] * _2xy; 

        T ul = (xl * kl + Plx) * cameral[0] + cameral[1]; 
        T vl = (yl * kl + ply) * cameral[2] + cameral[3];

        //[5]projectpoinst right 
        T objectsr[3];
        objectsr[0] = Rr[0][0] * objectr[0] + Rr[0][1] * objectr[1] + Rr[0][2] * objectr[2] ;
        objectsr[1] = Rr[1][0] * objectr[0] + Rr[1][1] * objectr[1] + Rr[1][2] * objectr[2] ;
        objectsr[2] = Rr[2][0] * objectr[0] + Rr[2][1] * objectr[1] + Rr[2][2] * objectr[2] ;

        objectsr[0] = objectsr[0] + Tr[0];
        objectsr[1] = objectsr[1] + Tr[1];
        objectsr[2] = objectsr[2] + Tr[2];

        T xr = objectsr[0]/objectsr[2];
        T yr = objectsr[1]/objectsr[2]; 

        T xr2 = xr*xr, yr2 = yr*yr;
        T rr2 = xr2 + yr2, _2xyr = T(2)*xr*yr;

        T kr  = (T(1) + distCoffr[4]*rr2*rr2*rr2 + distCoffr[1]*rr2*rr2 + distCoffr[0]*rr2 ); 
        T Prx = distCoffr[2]*_2xyr + distCoffr[3]*(rr2 + T(2)*xr2);
        T pry = distCoffr[2]*(rr2 + T(2) * yr2) + distCoffr[3] * _2xyr; 

        T ur = (xr * kr + Prx) * camerar[0] + camerar[1]; 
        T vr = (yr * kr + pry) * camerar[2] + camerar[3];
 
        //[6]result 
        residual[0] = (ul-pointsl[0]) * (ul-pointsl[0]) + (vl - pointsl[1])*(vl - pointsl[1]) + (ur - pointsr[0]) * (ur - pointsr[0]) + (ur - pointsr[1]) * (ur - pointsr[1]);
        
        return true; 
    }

public: 
    bool GetCalibrationParameters(const double* const StereoCalibrateParameters)
    {
        double rlx =  StereoCalibrateParameters[0];
        double rly =  StereoCalibrateParameters[1];
        double rlz =  StereoCalibrateParameters[2];
        double tlx =  StereoCalibrateParameters[3];
        double tly =  StereoCalibrateParameters[4];
        double tlz =  StereoCalibrateParameters[5];
        double rrx =  StereoCalibrateParameters[6];
        double rry =  StereoCalibrateParameters[7];
        double rrz =  StereoCalibrateParameters[8];
        double trx =  StereoCalibrateParameters[9];
        double trys = StereoCalibrateParameters[10];
        double trz =  StereoCalibrateParameters[11];

        cv::Mat rvecl = (cv::Mat_<double>(1,3)<< rlx, rly, rlz);
        cv::Mat tvecl = (cv::Mat_<double>(3,1)<< rrx, rry, rrz); 
        cv::Mat rvecr = (cv::Mat_<double>(1,3)<< tlx, tly, tlz);
        cv::Mat tvecr = (cv::Mat_<double>(3,1)<< trx, trys, trz);     
    
        std::cout<< "rvecl:" << rvecl << std::endl; 
        std::cout<< "tvecl:" << tvecl << std::endl;
        std::cout<< "rvecr:" << rvecr << std::endl; 
        std::cout<< "tvecr:" << tvecr << std::endl;           
       
        cv::Mat Rl, Rr; 
        Rodrigues(rvecl, Rl);
        Rodrigues(rvecr, Rr);

        cv::Mat Tl, Tr; 
        Tl = tvecl; 
        Tr = tvecr; 

        cv::Mat R, T; 
        Mat Rl_inv = Rl.inv(DECOMP_LU); 
        std::cout<<"Rr: " <<  Rr << std::endl; 
        std::cout<<"Rl: " <<  Rl << std::endl; 
        R = Rr * Rl_inv; 
        T = Tr - R * Tl; 

        std::cout<<"R: " << R << std::endl;
        std::cout<< "T: " << T << std::endl;

        std::string common_patho = "../cameras_params/w/";
        std::string pathRR = common_patho + "_R.xml";
        cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite << "_R" << R ;
            fswrite.release();
        } 
      
        std::string pathTT = common_patho + "_T.xml";
        fswrite.open(pathTT, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite <<"_T" << T ;
            fswrite.release();
        }
    }

private:
    double _imageu; 
    double _imagev; 
    double _objectx; 
    double _objecty; 
    double _objectz; 
    double _imageur; 
    double _imagevr; 
    double _objectxr; 
    double _objectyr; 
    double _objectzr;
    double _m00; 
    double _m02;
    double _m11;
    double _m12; 
    double _mr00;
    double _mr02;
    double _mr11;
    double _mr12; 
    double _d00;
    double _d10;
    double _d20;
    double _d30;
    double _d40;
    double _dr00;
    double _dr10;
    double _dr20;
    double _dr30;
    double _dr40;

    //cv::Mat _cameraMatrixl; 
    //cv::Mat _distCoeffsl; 
    //cv::Mat _cameraMatrixr; 
    //cv::Mat _distCoeffsr; 
    
};

struct StereoCalibrationResidual
{
    StereoCalibrationResidual(
                              double imageu, 
                              double imagev, 
                              double objectx, 
                              double objecty, 
                              double objectz, 
                              double imageur, 
                              double imagevr, 
                              double objectxr, 
                              double objectyr, 
                              double objectzr, 
                              double m00, 
                              double m02, 
                              double m11,
                              double m12, 
                              double d00, 
                              double d10, 
                              double d20, 
                              double d30, 
                              double d40, 
                              double mr00, 
                              double mr02,
                              double mr11,  
                              double mr12, 
                              double dr00, 
                              double dr10, 
                              double dr20, 
                              double dr30, 
                              double dr40  ): _imageu(imageu), _imagev(imagev), _objectx(objectx), _objecty(objecty), _objectz(objectz),_imageur(imageur), _imagevr(imagevr), _objectxr(objectxr), _objectyr(objectyr), _objectzr(objectzr), _m00(m00), _m02(m02), _m11(m11),_m12(m12), _d00(d00), _d10(d10), _d20(d20), _d30(d30),_d40(d40), _mr00(mr00), _mr02(mr02), _mr11(mr11), _mr12(mr12), _dr00(dr00), _dr10(dr10),_dr20(dr20),_dr30(dr30),_dr40(dr40)
    {
    
    }

    StereoCalibrationResidual()
    {

    }

    bool operator() (const double* const StereoCalibrateParameters, double* residual) const 
    {
        double rlx =  StereoCalibrateParameters[0];
        double rly =  StereoCalibrateParameters[1];
        double rlz =  StereoCalibrateParameters[2];
        double tlx =  StereoCalibrateParameters[3];
        double tly =  StereoCalibrateParameters[4];
        double tlz =  StereoCalibrateParameters[5];
        double rrx =  StereoCalibrateParameters[6];
        double rry =  StereoCalibrateParameters[7];
        double rrz =  StereoCalibrateParameters[8];
        double trx =  StereoCalibrateParameters[9];
        double trys = StereoCalibrateParameters[10];
        double trz =  StereoCalibrateParameters[11];

        cv::Mat rvecl = (cv::Mat_<double>(1,3)<< rlx, rly, rlz);
        cv::Mat tvecl = (cv::Mat_<double>(1,3)<< rrx, rry, rrz); 
        cv::Mat rvecr = (cv::Mat_<double>(1,3)<< tlx, tly, tlz);
        cv::Mat tvecr = (cv::Mat_<double>(1,3)<< trx, trys, trz);         
       
        cv::Mat Rl, Rr; 
        Rodrigues(rvecl, Rl);
        Rodrigues(rvecr, Rr);

        double _objectzs  = _objectz; 
        double _objectzrs = _objectzr;  
        cv::Point2d imagepointsl(_imageu, _imagev);
        cv::Point3d objpointsl(_objectx, _objecty, _objectzs);  
        cv::Point2d imagepointsr(_imageur, _imagevr);
        cv::Point3d objpointsr(_objectxr,_objectyr,  _objectzrs);

        std::vector<cv::Point2d> projectedPointsl; 
        std::vector<cv::Point3d> objectPointsl;
        std::vector<cv::Point2d> projectedPointsr; 
        std::vector<cv::Point3d> objectPointsr; 

        objectPointsl.push_back(objpointsl);
        objectPointsr.push_back(objpointsr); 

        cv::Mat _cameraMatrixl(3, 3, CV_32FC1); 
        cv::Mat _distCoeffsl(1, 5, CV_32FC1); 
        cv::Mat _cameraMatrixr(3, 3, CV_32FC1); 
        cv::Mat _distCoeffsr(1, 5, CV_32FC1); 

        _cameraMatrixl.setTo(0); 
        _cameraMatrixr.setTo(0); 
        _distCoeffsl.setTo(0);
        _distCoeffsr.setTo(0); 
        
        float* cameraDatal = (float*) _cameraMatrixl.data; 
        int    cameraStepl = _cameraMatrixl.step; 
        float* cameraDatar = (float*) _cameraMatrixr.data; 
        int    cameraStepr = _cameraMatrixr.step;

        //std::cout<<"_m11: " << _m11 <<",cameraStepl: " << cameraStepl<< std::endl;   
        cameraDatal[0 * (cameraStepl/4) + 0]   = _m00; 
        cameraDatal[0 * (cameraStepl/4) + 2]   = _m02;
        cameraDatal[1 * (cameraStepl/4) + 1]   = _m11;
        cameraDatal[1 * (cameraStepl/4) + 2]   = _m12; 
        cameraDatar[0 * (cameraStepr/4) + 0]   = _mr00; 
        cameraDatar[0 * (cameraStepr/4) + 2]   = _mr02;
        cameraDatar[1 * (cameraStepr/4) + 1]   = _mr11;
        cameraDatar[1 * (cameraStepr/4) + 2]   = _mr12; 

        //std::cout<<"_d00: " << _d00 << "_d10: " << _d10 << "_d20: " <<_d20 << "_d30: " << _d30 << std::endl; 
        _distCoeffsl.at<float>(0,0) = _d00;
        _distCoeffsl.at<float>(0,1) = _d10;
        _distCoeffsl.at<float>(0,2) = _d20;
        _distCoeffsl.at<float>(0,3) = _d30;
        _distCoeffsl.at<float>(0,4) = _d40;
        _distCoeffsr.at<float>(0,0) = _dr00;
        _distCoeffsr.at<float>(0,1) = _dr10;
        _distCoeffsr.at<float>(0,2) = _dr20;
        _distCoeffsr.at<float>(0,3) = _dr30;
        _distCoeffsr.at<float>(0,4) = _dr40;
     
        //std::cout<< "rvecl:" << rvecl << std::endl; 
        //std::cout<< "tvecl:" << tvecl << std::endl;
        //std::cout<< "rvecr:" << rvecr << std::endl; 
        //std::cout<< "tvecr:" << tvecr << std::endl;

        //std::cout<< "_cameraMatrixl: " << _cameraMatrixl << std::endl; 
        //std::cout<< "_distCoeffsl: " << _distCoeffsl << std::endl;

        //std::cout<< "_cameraMatrixr: " << _cameraMatrixr << std::endl; 
        //std::cout<< "_distCoeffsr: " << _distCoeffsr << std::endl;

        cv::projectPoints(objectPointsl, rvecl, tvecl, _cameraMatrixl, _distCoeffsl, projectedPointsl);
        cv::projectPoints(objectPointsr, rvecr, tvecr, _cameraMatrixr, _distCoeffsr, projectedPointsr);
       
        //std::cout <<"objectPointsl.x: " << objectPointsl[0].x << "objectPointsl.y: " << objectPointsl[0].y << "objectPointsl.z: " << objectPointsl[0].z << "image.x: " <<  imagepointsl.x <<"," << imagepointsl.y << "," << projectedPointsl[0].x << "," << projectedPointsl[0].y << std::endl; 
        //std::cout <<"objectPointsr.x: " << objectPointsr[0].x << "objectPointsr.y: " << objectPointsr[0].y << "objectPointsr.z: " << objectPointsr[0].z << "image.x: " <<  imagepointsr.x <<"," << imagepointsr.y << "," << projectedPointsr[0].x << "," << projectedPointsr[0].y << std::endl; 

        //return false; 
        double distancelx = (imagepointsl.x - projectedPointsl[0].x) * (imagepointsl.x - projectedPointsl[0].x);
        double distancely = (imagepointsl.y - projectedPointsl[0].y) * (imagepointsl.y - projectedPointsl[0].y); 
        double distancerx = (imagepointsr.x - projectedPointsr[0].x) * (imagepointsr.x - projectedPointsr[0].x); 
        double distancery = (imagepointsr.y - projectedPointsr[0].y) * (imagepointsr.y - projectedPointsr[0].y);
 
        residual[0] = sqrt(distancelx + distancely) + sqrt(distancerx + distancery); 
        //std::cout<<"rlx: " << rlx << std::endl; 
        //std::cout<<"residual[0] : " << residual[0] << std::endl; 
        return true; 
    }

public: 
    bool GetCalibrationParameters(const double* const StereoCalibrateParameters)
    {
        double rlx =  StereoCalibrateParameters[0];
        double rly =  StereoCalibrateParameters[1];
        double rlz =  StereoCalibrateParameters[2];
        double tlx =  StereoCalibrateParameters[3];
        double tly =  StereoCalibrateParameters[4];
        double tlz =  StereoCalibrateParameters[5];
        double rrx =  StereoCalibrateParameters[6];
        double rry =  StereoCalibrateParameters[7];
        double rrz =  StereoCalibrateParameters[8];
        double trx =  StereoCalibrateParameters[9];
        double trys = StereoCalibrateParameters[10];
        double trz =  StereoCalibrateParameters[11];

        cv::Mat rvecl = (cv::Mat_<double>(1,3)<< rlx, rly, rlz);
        cv::Mat tvecl = (cv::Mat_<double>(3,1)<< rrx, rry, rrz); 
        cv::Mat rvecr = (cv::Mat_<double>(1,3)<< tlx, tly, tlz);
        cv::Mat tvecr = (cv::Mat_<double>(3,1)<< trx, trys, trz);     
    
        std::cout<< "rvecl:" << rvecl << std::endl; 
        std::cout<< "tvecl:" << tvecl << std::endl;
        std::cout<< "rvecr:" << rvecr << std::endl; 
        std::cout<< "tvecr:" << tvecr << std::endl;           
       
        cv::Mat Rl, Rr; 
        Rodrigues(rvecl, Rl);
        Rodrigues(rvecr, Rr);

        cv::Mat Tl, Tr; 
        Tl = tvecl; 
        Tr = tvecr; 

        cv::Mat R, T; 
        Mat Rl_inv = Rl.inv(DECOMP_LU); 
        std::cout<<"Rr: " <<  Rr << std::endl; 
        std::cout<<"Rl: " <<  Rl << std::endl; 
        R = Rr * Rl_inv; 
        T = Tr - R * Tl; 

        std::cout<<"R: " << R << std::endl;
        std::cout<< "T: " << T << std::endl;

        std::string common_patho = "../cameras_params/w/";
        std::string pathRR = common_patho + "_R.xml";
        cv::FileStorage fswrite(pathRR, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite << "_R" << R ;
            fswrite.release();
        } 
      
        std::string pathTT = common_patho + "_T.xml";
        fswrite.open(pathTT, cv::FileStorage::WRITE);
        if(fswrite.isOpened())
        {
            fswrite <<"_T" << T ;
            fswrite.release();
        }
    }

private:
    double _imageu; 
    double _imagev; 
    double _objectx; 
    double _objecty; 
    double _objectz; 
    double _imageur; 
    double _imagevr; 
    double _objectxr; 
    double _objectyr; 
    double _objectzr;
    double _m00; 
    double _m02;
    double _m11;
    double _m12; 
    double _mr00;
    double _mr02;
    double _mr11;
    double _mr12; 
    double _d00;
    double _d10;
    double _d20;
    double _d30;
    double _d40;
    double _dr00;
    double _dr10;
    double _dr20;
    double _dr30;
    double _dr40;

    //cv::Mat _cameraMatrixl; 
    //cv::Mat _distCoeffsl; 
    //cv::Mat _cameraMatrixr; 
    //cv::Mat _distCoeffsr; 
    
};

struct BackCrossResidual 
{
    BackCrossResidual(double x, double y, double f, double X, double Y, double Z)
        :_x(x),_y(y),_f(f),_X(X),_Y(Y),_Z(Z)
    {

    }
   
    BackCrossResidual()
    {

    }

    template<typename T>
    bool operator() (const T* const pBackCrossParameters, T* residual) const
    {
        T dXs = pBackCrossParameters[0];
        T dYs = pBackCrossParameters[1]; 
        T dZs = pBackCrossParameters[2];
        T dPhi   = pBackCrossParameters[3];
        T dOmega = pBackCrossParameters[4];
        T dKappa = pBackCrossParameters[5];

        T a1  = cos(dPhi)*cos(dKappa) - sin(dPhi)*sin(dOmega)*sin(dKappa); 
        T a2  = -cos(dPhi)*sin(dKappa) - sin(dPhi)*sin(dOmega)*cos(dKappa); 
        T a3  = -sin(dPhi)*cos(dOmega); 
        T b1  = cos(dOmega)*sin(dKappa);
        T b2  = cos(dOmega)*cos(dKappa);
        T b3  = -sin(dOmega);
        T c1  = sin(dPhi)*cos(dKappa)  + cos(dPhi)*sin(dOmega)*sin(dKappa);
        T c2  = -sin(dPhi)*sin(dKappa) + cos(dPhi)*sin(dOmega)*cos(dKappa);
        T c3  = cos(dPhi)*cos(dOmega);

        residual[0] = T(_x) + T(_f) * T( ( a1*(_X-dXs) + b1*(_Y-dYs) + c1*(_Z-dZs) )/( a3*(_X-dXs) + b3*(_Y-dYs) + c3*(_Z-dZs) ) );

        residual[1] = T(_y) + T(_f) * T( ( a2*(_X-dXs) + b2*(_Y-dYs) + c2*(_Z-dZs) )/( a3*(_X-dXs) + b3*(_Y-dYs) + c3*(_Z-dZs) ) ); 

        return true; 

    }
    template<typename T>
    bool RotationTrans(T* pBackCrossParameters)
    {
        T dXs = pBackCrossParameters[0];
        T dYs = pBackCrossParameters[1]; 
        T dZs = pBackCrossParameters[2];
        T dPhi   = pBackCrossParameters[3];
        T dOmega = pBackCrossParameters[4];
        T dKappa = pBackCrossParameters[5];

        T a1  = cos(dPhi)*cos(dKappa) - sin(dPhi)*sin(dOmega)*sin(dKappa); 
        T a2  = -cos(dPhi)*sin(dKappa) - sin(dPhi)*sin(dOmega)*cos(dKappa); 
        T a3  = -sin(dPhi)*cos(dOmega); 
        T b1  = cos(dOmega)*sin(dKappa);
        T b2  = cos(dOmega)*cos(dKappa);
        T b3  = -sin(dOmega);
        T c1  = sin(dPhi)*cos(dKappa)  + cos(dPhi)*sin(dOmega)*sin(dKappa);
        T c2  = -sin(dPhi)*sin(dKappa) + cos(dPhi)*sin(dOmega)*cos(dKappa);
        T c3  = cos(dPhi)*cos(dOmega);

        std::cout<< a1 << ","<< a2 << ","<< a3 << "," << std::endl; 
        std::cout<< b1 << ","<< b2 << ","<< b3 << "," << std::endl;
        std::cout<< c1 << ","<< c2 << ","<< c3 << "," << std::endl;
        std::cout<< dXs << ","<< dYs << ","<< dZs << "," << std::endl;

        return true; 
    }

private :
    double _x; 
    double _y; 
    double _f; 
    double _X; 
    double _Y; 
    double _Z; 
};


class QuadraticCostFunction : public ceres::SizedCostFunction<1,/*number of residuals*/
                                                       1/*size of first parameters*/>
{
public: 
     virtual ~QuadraticCostFunction() {}

     virtual bool Evaluate(double const* const* parameters, 
                           double * residuals, 
                           double ** jacobians) const
     {
         double x     = parameters[0][0];
          
         residuals[0] = (x - 4) * (x - 8);//(5 - x)*(5 - x); 

         if(jacobians != NULL && jacobians[0] != NULL)
         {
             jacobians[0][0] = -1; 
         }
          
         return true; 
     }

}; 

#endif 
