#include "CeresBa.h"

CeresBa::CeresBa(void)
{

}

CeresBa::~CeresBa(void)
{

}

PointDBLXYZ CeresBa::ShiftingCoordinateValue(std::vector<PointDBLXYZ> & pointClouds,
	pcl::PointCloud<pcl::PointXYZ>::Ptr & cloud) {

	//std::cout << "input : " << std::endl; 
	//std::cout << "ptsize: " << pointClouds.size() << std::endl;
        //std::cout << "width : " << cloud->width << std::endl; 

	cloud->width = pointClouds.size();
	cloud->height = 1;
	cloud->is_dense = false;
	cloud->points.resize(cloud->width*cloud->height);

        //std::cout << "cloud->points.size : " << cloud->points.size() << std::endl; 
	//count max value
	for (int i = 0; i != pointClouds.size(); i++) {
		if (pointClouds[i].x>maxX)
			maxX = pointClouds[i].x;
		if (pointClouds[i].y>maxY)
			maxY = pointClouds[i].y;
		if (pointClouds[i].z>maxZ)
			maxZ = pointClouds[i].z;
	}

        //new a output as the offset
	PointDBLXYZ Offset;
	Offset.x = maxX;
	Offset.y = maxY;
	//Offset.z = maxZ;

        //all elevation values should be positive for further compare,responded by XiaoQi
        Offset.z = 0.0;

	//shift coordinate value of each point
	for (int i = 0; i != pointClouds.size(); i++) {
		cloud->points[i].x = pointClouds[i].x - Offset.x;
		cloud->points[i].y = pointClouds[i].y - Offset.y;
	        cloud->points[i].z = pointClouds[i].z - Offset.z;
                //std::cout << cloud->points[i].x << " , " << cloud->points[i].y << " , " <<cloud->points[i].z << std::endl; 
	}

	return Offset;

}

PointDBLXYZT CeresBa::ShiftingCoordinateValueT(std::vector<PointDBLXYZT> & pointClouds,
	std::vector<PointDBLXYZT> & cloud) {

        cloud.resize( pointClouds.size() );
	//count max value
	for (int i = 0; i < pointClouds.size(); i++) {
		if (pointClouds[i].x>maxX)
			maxX = pointClouds[i].x;
		if (pointClouds[i].y>maxY)
			maxY = pointClouds[i].y;
		if (pointClouds[i].z>maxZ)
			maxZ = pointClouds[i].z;
	}

        //new a output as the offset
	PointDBLXYZT Offset;
	Offset.x = maxX;
	Offset.y = maxY;
	//Offset.z = maxZ;

        //all elevation values should be positive for further compare,responded by XiaoQi
        Offset.z = 0.0;
        Offset.time = 0.0; 

	//shift coordinate value of each point
	for (int i = 0; i < pointClouds.size(); i++) {
		cloud[i].x    = pointClouds[i].x - Offset.x;
		cloud[i].y    = pointClouds[i].y - Offset.y;
	        cloud[i].z    = pointClouds[i].z - Offset.z;
                cloud[i].time = pointClouds[i].time - Offset.time; 
                //std::cout << cloud->points[i].x << " , " << cloud->points[i].y << " , " <<cloud->points[i].z << std::endl; 
	}

	return Offset;

}

void  TwoMatrixsMultiplyLane(double returnMatrix[3][3], const  double A[3][3], const double B[3][3])
{
    returnMatrix[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
    returnMatrix[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
    returnMatrix[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

    returnMatrix[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
    returnMatrix[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
    returnMatrix[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

    returnMatrix[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
    returnMatrix[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
    returnMatrix[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

void matrixMultiVectorLane(double coordinate[3], double R[3][3],double T[3])
{
    double b[3];
    b[0] = coordinate[0];
    b[1] = coordinate[1];
    b[2] = coordinate[2];
    coordinate[0] = R[0][0] * b[0] + R[0][1] * b[1] + R[0][2] * b[2];
    coordinate[1] = R[1][0] * b[0] + R[1][1] * b[1] + R[1][2] * b[2];
    coordinate[2] = R[2][0] * b[0] + R[2][1] * b[1] + R[2][2] * b[2];
    coordinate[0] = coordinate[0] + T[0];
    coordinate[1] = coordinate[1] + T[1];
    coordinate[2] = coordinate[2] + T[2];
}

void matrixMultiVector(double coordinate[3], double R[3][3],double T[3])
{
    double b[3];
    b[0] = coordinate[0];
    b[1] = coordinate[1];
    b[2] = coordinate[2];
    coordinate[0] = R[0][0] * b[0] + R[0][1] * b[1] + R[0][2] * b[2];
    coordinate[1] = R[1][0] * b[0] + R[1][1] * b[1] + R[1][2] * b[2];
    coordinate[2] = R[2][0] * b[0] + R[2][1] * b[1] + R[2][2] * b[2];
    T[0] = coordinate[0]; 
    T[1] = coordinate[1];
    T[2] = coordinate[2];
}

void matrixInv(double matrixin[3][3], double matrixinv[3][3])
{
    matrixinv[0][0] =  matrixin[0][0]; 
    matrixinv[0][1] =  matrixin[1][0];
    matrixinv[0][2] =  matrixin[2][0];

    matrixinv[1][0] =  matrixin[0][1];
    matrixinv[1][1] =  matrixin[1][1];
    matrixinv[1][2] =  matrixin[2][1];

    matrixinv[2][0] =  matrixin[0][2];
    matrixinv[2][1] =  matrixin[1][2];
    matrixinv[2][2] =  matrixin[2][2];
}

