//
// Created by konaco on 20. 8. 7..
//

#ifndef LIOTEST_GEOFUNC_H
#define LIOTEST_GEOFUNC_H
#include <cmath>
#include <eigen3/Eigen/Dense>

extern "C" {
//Convert Earth-Centered-Earth-Fixed (ECEF) to lat, Lon, Altitude
//x, y, z in meters
//lat and lon in degrees, and altitude in meters
void ecef2llh(double x, double y, double z, double *lat, double *lon, double *h);

//Convert Lat, Lon, Altitude to Earth-Centered-Earth-Fixed (ECEF)
//lat, lon (degs) and alt (m)
//x, y, z in meters
void llh2ecef(double lat, double lon, double h, double *x, double *y, double *z);

void llh2local(double lat, double lon, double h,
               double latOrg, double lonOrg, double hOrg,
               double *n, double *e, double *d);

void local2llh(double n, double e, double d,
               double latOrg, double lonOrg, double hOrg,
               double *lat, double *lon, double *h);

void ecef2local(double x, double y, double z,
                double orgX, double orgY, double orgZ,
                double *n, double *e, double *d);

void local2ecef(double n, double e, double d, double orgX, double orgY, double orgZ, double *x, double *y, double *z);

// eulr (rads, rads, rads)
void eulr2qua(Eigen::Vector3d eulr, Eigen::Vector4d* q);

void qua2dcm(Eigen::Vector4d q, Eigen::Matrix3d* dcmbn);

//
void cross(Eigen::Vector3d aVec, Eigen::Vector3d bVec, Eigen::Vector3d* Vect);
void ang2qua(Eigen::Vector3d three_angle, Eigen::Vector4d *qua);
void qua_update(Eigen::Vector4d q_a, Eigen::Vector4d q_b, Eigen::Vector4d* q);

// degs
double gravityWgs84(double lat, double lon, double h);

void getEulerYPR(Eigen::Matrix3d m_el,double *yaw, double *pitch, double *roll);
void getRPY(Eigen::Matrix3d m_el,double *roll, double *pitch,double *yaw );
} // extern C
#endif //LIOTEST_GEOFUNC_H
