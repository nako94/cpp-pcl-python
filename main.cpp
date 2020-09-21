#include <iostream>
#include <pcl/common/transforms.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/registration/ndt.h>
#include <pcl/filters/voxel_grid.h>
#include <eigen3/Eigen/Dense>
#include "geofunc.h"
using namespace std;
using namespace Eigen;
int main() {
    cout<<"testest"<<endl;
    pcl::PointCloud<pcl::PointXYZ> scan,output;
    pcl::PointCloud<pcl::PointXYZ>::Ptr output_ptr (new pcl::PointCloud<pcl::PointXYZ>());


    pcl::PointXYZ org, testpoint;
    org.x = 0;
    org.y = 0;
    org.z = 0;
    scan.push_back(org);
    testpoint.x=1;
    testpoint.y=1;
    testpoint.z=1;
    scan.push_back(testpoint);

    Eigen::Matrix4f update_tf(Eigen::Matrix4f::Identity());
    Eigen::Translation3f tl_update(0, 0, 0);  // tl: translation roll,pitch, yaw -> n,e,d
    Eigen::AngleAxisf rot_x_update(0, Eigen::Vector3f::UnitX());  // rot: rotation
    Eigen::AngleAxisf rot_y_update(45*M_PI/180, Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf rot_z_update(90*M_PI/180, Eigen::Vector3f::UnitZ());
//        Eigen::Translation3f tl_update(DataBuf.refe, DataBuf.refn, DataBuf.refu);  // tl: translation
//        Eigen::AngleAxisf rot_x_update(DataBuf.refroll, Eigen::Vector3f::UnitX());  // rot: rotation
//        Eigen::AngleAxisf rot_y_update(DataBuf.refpitch, Eigen::Vector3f::UnitY());
//        Eigen::AngleAxisf rot_z_update(DataBuf.refyaw, Eigen::Vector3f::UnitZ());
    update_tf = (tl_update * rot_z_update * rot_y_update * rot_x_update).matrix();    //transformation matrix
    //=>in order of roll,pitch,yaw,tl
    pcl::transformPointCloud (scan, output, update_tf);
    for (pcl::PointCloud<pcl::PointXYZ>::const_iterator item = output.begin(); item != output.end(); item++)
    {
        pcl::PointXYZ pp;
        pp.x = (double) item->x;
        pp.y = (double) item->y;
        pp.z = (double) item->z;
cout<<"pp.x : "<<pp.x<<", pp.y : "<<pp.y<<", pp.z : "<<pp.z<<endl;
    }
    Vector4d attitude;
    Vector3d euler;
    attitude<<0.6426  ,  0.3532  ,  0.5809,   -0.3532;
    Matrix3d dcm;double r,pit,y;
//    attitude=-attitude;
////    attitude(0)=-attitude(0);
    qua2dcm(attitude,&dcm); //body to navi dcm
    getRPY(dcm,&r,&pit,&y);
    euler<<r,pit,y;
    cout<<r*180/M_PI<<endl;
    cout<<pit*180/M_PI<<endl;
    cout<<y*180/M_PI<<endl;
//    pcl::PointCloud<pcl::PointXYZ>::Ptr pcmapsaved_ptr (new pcl::PointCloud<pcl::PointXYZ>);
//    if (pcl::io::loadPCDFile<pcl::PointXYZ> ("pcmapsaved.pcd", *pcmapsaved_ptr) == -1) //* load the file
//    {
//        PCL_ERROR ("Couldn't read file pcmapsaved.pcd \n");
//    }
//    std::cout << "pcmap Loaded "
//              << pcmapsaved_ptr->width * pcmapsaved_ptr->height
//              << " data points from test_pcd.pcd with the following fields: "
//              << std::endl;
    return 0;
}


//#include <iostream>
//#include <eigen3/Eigen/Dense>
//#include <pcl/registration/ndt.h>
//using namespace std;
//using namespace Eigen;
//int main() {
//    cout<<"testest"<<endl;
//}