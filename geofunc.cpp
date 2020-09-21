//
// Created by konaco on 20. 8. 7..
//

#include "geofunc.h"
#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

const double a = 6378137.0;              //WGS-84 semi-major axis
const double e2 = 6.6943799901377997e-3;  //WGS-84 first eccentricity squared
const double a1 = 4.2697672707157535e+4;  //a1 = a*e2
const double a2 = 1.8230912546075455e+9;  //a2 = a1*a1
const double a3 = 1.4291722289812413e+2;  //a3 = a1*e2/2
const double a4 = 4.5577281365188637e+9;  //a4 = 2.5*a2
const double a5 = 4.2840589930055659e+4;  //a5 = a1+a3
const double a6 = 9.9330562000986220e-1;  //a6 = 1-e2

extern "C" {
//Convert Earth-Centered-Earth-Fixed (ECEF) to lat, Lon, Altitude
//x, y, z in meters
//lat and lon in degrees, and altitude in meters
void ecef2llh(double x, double y, double z, double *lat, double *lon, double *h){
    double zp,w2,w,r2,r,s2,c2,s,c,ss;
    double g,rg,rf,u,v,m,f,p;
    double latBuf, lonBuf, hBuf;
    zp = std::abs( z );
    w2 = x*x + y*y;
    w = sqrt( w2 );
    r2 = w2 + z*z;
    r = sqrt( r2 );
    lonBuf = atan2( y, x );       //Lon (final)
    s2 = z*z/r2;
    c2 = w2/r2;
    u = a2/r;
    v = a3 - a4/r;
    if( c2 > 0.3 ){
        s = ( zp/r )*( 1.0 + c2*( a1 + u + s2*v )/r );
        latBuf = asin( s );      //Lat
        ss = s*s;
        c = sqrt( 1.0 - ss );
    }
    else{
        c = ( w/r )*( 1.0 - s2*( a5 - u - c2*v )/r );
        latBuf = acos( c );      //Lat
        ss = 1.0 - c*c;
        s = sqrt( ss );
    }
    g = 1.0 - e2*ss;
    rg = a/sqrt( g );
    rf = a6*rg;
    u = w - rg*c;
    v = zp - rf*s;
    f = c*u + s*v;
    m = c*v - s*u;
    p = m/( rf/g + f );
    latBuf = latBuf + p;      //Lat
    hBuf = f + m*p/2.0;     //Altitude
    if( z < 0.0 ){
        latBuf *= -1.0;     //Lat
    }
    *lat = latBuf*180/M_PI;
    *lon = lonBuf*180/M_PI;
    *h = hBuf;
}

//Convert Lat, Lon, Altitude to Earth-Centered-Earth-Fixed (ECEF)
//lat, lon (degs) and alt (m)
//x, y, z in meters
void llh2ecef(double lat, double lon, double h, double *x, double *y, double *z) {
    double latBuf, lonBuf, n;
    latBuf = lat*M_PI/180;
    lonBuf = lon*M_PI/180;
    n = a/sqrt( 1 - e2*sin( latBuf )*sin( latBuf ) );
    *x = ( n + h )*cos( latBuf )*cos( lonBuf );    //ECEF x
    *y = ( n + h )*cos( latBuf )*sin( lonBuf );    //ECEF y
    *z = ( n*(1 - e2 ) + h )*sin( latBuf );          //ECEF z
}

void llh2local(double lat, double lon, double h,
               double latOrg, double lonOrg, double hOrg,
               double *n, double *e, double *d){
    Matrix3d c_e2n;
    Vector3d locBuf, locBufp;
    double latOrgRad, lonOrgRad;
    double x, y, z;
    double orgX, orgY, orgZ;
    llh2ecef(latOrg, lonOrg, hOrg, &orgX, &orgY, &orgZ);	// origin (llh) -> (ECEF)
    llh2ecef(lat, lon, h, &x, &y, &z);

    latOrgRad = latOrg*M_PI/180;
    lonOrgRad = lonOrg*M_PI/180;
    c_e2n << -sin(latOrgRad)*cos(lonOrgRad), -sin(latOrgRad)*sin(lonOrgRad), cos(latOrgRad),
            -sin(lonOrgRad), cos(lonOrgRad), 0,
            -cos(latOrgRad)*cos(lonOrgRad), -cos(latOrgRad)*sin(lonOrgRad), -sin(latOrgRad);

    locBuf << x - orgX, y - orgY, z - orgZ;
    locBufp = c_e2n*locBuf;
    *n = locBufp(0);
    *e = locBufp(1);
    *d = locBufp(2);
}

void local2llh(double n, double e, double d,
               double latOrg, double lonOrg, double hOrg,
               double *lat, double *lon, double *h){

    Matrix3d c_e2n, c_n2e;
    Vector3d locBuf, centerBuf, ecefBuf;

    double orgX, orgY, orgZ;
    double latOrgRad, lonOrgRad;
    double latBuf, lonBuf, hBuf;

    locBuf << n, e, d;
    llh2ecef(latOrg, lonOrg, hOrg, &orgX, &orgY, &orgZ);	// origin (llh) -> (ECEF)
    centerBuf << orgX, orgY, orgZ;
    latOrgRad = latOrg*M_PI/180;
    lonOrgRad = lonOrg*M_PI/180;
    c_e2n << -sin(latOrgRad)*cos(lonOrgRad), -sin(latOrgRad)*sin(lonOrgRad), cos(latOrgRad),
            -sin(lonOrgRad), cos(lonOrgRad), 0,
            -cos(latOrgRad)*cos(lonOrgRad), -cos(latOrgRad)*sin(lonOrgRad), -sin(latOrgRad);
    c_n2e = c_e2n.transpose();
    ecefBuf = centerBuf + c_n2e*locBuf;
    ecef2llh(ecefBuf(0), ecefBuf(1), ecefBuf(2), &latBuf, &lonBuf, &hBuf);
    *lat = latBuf;
    *lon = lonBuf;
    *h = hBuf;

}

void ecef2local(double x, double y, double z,
                double orgX, double orgY, double orgZ,
                double *n, double *e, double *d) {

    double latOrg, lonOrg, latOrgRad, lonOrgRad, h;
    ecef2llh(orgX, orgY, orgZ, &latOrg, &lonOrg, &h);

    Matrix3d c_e2n;
    Vector3d locBuf, locBufp;

    latOrgRad = latOrg*M_PI / 180;
    lonOrgRad = lonOrg*M_PI / 180;
    c_e2n << -sin(latOrgRad)*cos(lonOrgRad), -sin(latOrgRad)*sin(lonOrgRad), cos(latOrgRad),
            -sin(lonOrgRad), cos(lonOrgRad), 0,
            -cos(latOrgRad)*cos(lonOrgRad), -cos(latOrgRad)*sin(lonOrgRad), -sin(latOrgRad);

    locBuf << x - orgX, y - orgY, z - orgZ;
    locBufp = c_e2n*locBuf;
    *n = locBufp(0);
    *e = locBufp(1);
    *d = locBufp(2);
}

void local2ecef(double n, double e, double d, double orgX, double orgY, double orgZ, double *x, double *y, double *z)
{
    Matrix3d c_e2n, c_n2e;
    Vector3d locBuf, centerBuf, ecefBuf;

    double latOrg, lonOrg, hOrg;
    double latOrgRad, lonOrgRad;
    double latBuf, lonBuf, hBuf, nn;


    locBuf << n, e, d;
    ecef2llh( orgX, orgY, orgZ,&latOrg,&lonOrg, &hOrg);
    centerBuf << orgX, orgY, orgZ;
    latOrgRad = latOrg*M_PI/180;
    lonOrgRad = lonOrg*M_PI/180;
    c_e2n << -sin(latOrgRad)*cos(lonOrgRad), -sin(latOrgRad)*sin(lonOrgRad), cos(latOrgRad),
            -sin(lonOrgRad), cos(lonOrgRad), 0,
            -cos(latOrgRad)*cos(lonOrgRad), -cos(latOrgRad)*sin(lonOrgRad), -sin(latOrgRad);
    c_n2e = c_e2n.transpose();
    ecefBuf = centerBuf + c_n2e*locBuf;
    *x = ecefBuf(0);
    *y= ecefBuf(1);
    *z = ecefBuf(2);
}

void eulr2qua(Vector3d eulr, Vector4d* q){
    double phi, theta, psi;
    double cpsi2, spsi2, cthe2, sthe2, cphi2, sphi2;

    phi = eulr(0)*M_PI/180;
    theta = eulr(1)*M_PI/180;
    psi = eulr(2)*M_PI/180;
    cpsi2 = cos(psi/2);
    spsi2 = sin(psi/2);
    cthe2 = cos(theta/2);
    sthe2 = sin(theta/2);
    cphi2 = cos(phi/2);
    sphi2 = sin(phi/2);

    *q << cphi2*cthe2*cpsi2 + sphi2*sthe2*spsi2,
            sphi2*cthe2*cpsi2 - cphi2*sthe2*spsi2,
            cphi2*sthe2*cpsi2 + sphi2*cthe2*spsi2,
            cphi2*cthe2*spsi2 + sphi2*sthe2*cpsi2;
}

void qua2dcm(Vector4d q, Matrix3d* dcmbn){
    // from body to nav.
    double a, b, c, d;
    a = q(0);
    b = q(1);
    c = q(2);
    d = q(3);

    *dcmbn << a*a + b*b - c*c - d*d,			2*(b*c - a*d),			2*(b*d + a*c),
            2*(b*c + a*d),			a*a - b*b + c*c - d*d,			2*(c*d - a*b),
            2*(b*d - a*c),					2*(c*d + a*b),	a*a - b*b - c*c + d*d;
}

void cross(Vector3d aVec, Vector3d bVec, Vector3d* Vect){
    double c1, c2, c3;

    c1 = aVec(1) * bVec(2) - aVec(2) * bVec(1);
    c2 = aVec(2) * bVec(0) - aVec(0) * bVec(2);
    c3 = aVec(0) * bVec(1) - aVec(1) * bVec(0);

    *Vect << c1, c2, c3;
}


void qua_update(Vector4d q_a, Vector4d q_b, Vector4d* q){
    double a, b, c, d;
    Vector4d qua;

    a = q_a(0)*q_b(0) - q_a(1)*q_b(1) - q_a(2)*q_b(2) - q_a(3)*q_b(3);
    b = q_a(0)*q_b(1) + q_a(1)*q_b(0) + q_a(2)*q_b(3) - q_a(3)*q_b(2);
    c = q_a(0)*q_b(2) + q_a(2)*q_b(0) + q_a(3)*q_b(1) - q_a(1)*q_b(3);
    d = q_a(0)*q_b(3) + q_a(3)*q_b(0) + q_a(1)*q_b(2) - q_a(2)*q_b(1);

    qua << a, b, c, d;
    qua.normalize();

    *q << qua;
}

void ang2qua(Vector3d three_angle, Vector4d *qua){

    double ang, a, b, c, d;
    Vector4d q;

    ang = sqrt(three_angle(0)*three_angle(0) + three_angle(1)*three_angle(1) + three_angle(2)*three_angle(2));

    if (abs(ang) > 0.0001){
        a = cos(ang/2);
        b = sin(ang/2) * (three_angle(0)/ang);
        c = sin(ang/2) * (three_angle(0)/ang);
        d = sin(ang/2) * (three_angle(0)/ang);
    } else {
        a = cos(ang/2);
        b = three_angle(0)/2;
        c = three_angle(0)/2;
        d = three_angle(0)/2;
    }

    q << a, b, c, d;
    q.normalize();
    *qua << q;
}



double gravityWgs84(double lat, double lon, double h){
    double latRad, lonRad, g0, gBuf, RN, RE, R0;
    const double R = 6378137.0;  //m
    const double e = 0.0818191908426;
    latRad = lat*M_PI/180;
    lonRad = lon*M_PI/180;
    RN = R*(1 - e*e)/pow(1 - e*e*pow(sin(latRad),2),1.5);
    RE = R/sqrt(1 - e*e*pow(sin(latRad),2));
    R0 = sqrt(RN*RE);
    g0 = 9.780318*(1 + 5.3024*1e-3*pow(sin(latRad),2) - 5.9*1e-6*pow(sin(2*lonRad),2));
    gBuf = g0/pow(1 + h/R0,2);
    return gBuf;
}
void getEulerYPR(Matrix3d m_el,double *yaw, double *pitch, double *roll)
{
    struct Euler
    {
        double yaw;
        double pitch;
        double roll;
    };

    Euler euler_out;
    //get the pointer to the raw data

    // Check that pitch is not at a singularity
    // Check that pitch is not at a singularity
    if (abs(m_el(2,0)) >= 1)
    {
        euler_out.yaw = 0;

        // From difference of angles formula
        if (m_el(2,0) < 0)  //gimbal locked down
        {
            double delta = atan2(m_el(0,1),m_el(0,2));
            euler_out.pitch = M_PI / (2.0);
            euler_out.roll = delta;
        }
        else // gimbal locked up
        {
            double delta = atan2(-m_el(0,1),-m_el(0,2));
            euler_out.pitch = -M_PI / (2.0);
            euler_out.roll = delta;
        }
    }
    else
    {
        euler_out.pitch = - asin(m_el(2,0));

        euler_out.roll = atan2(m_el(2,1)/cos(euler_out.pitch), m_el(2,2)/cos(euler_out.pitch));

        euler_out.yaw = atan2(m_el(1,0)/cos(euler_out.pitch), m_el(0,0)/cos(euler_out.pitch));
    }
    *pitch = euler_out.pitch;
    *roll= euler_out.roll;
    *yaw = euler_out.yaw;

}

void getRPY(Matrix3d m_el,double *roll, double *pitch,double *yaw )
{
    getEulerYPR(m_el,yaw, pitch, roll );
}
} // extern C
