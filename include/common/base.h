#pragma once

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;
namespace Polycs{ class polycrystal; }

#define SQR2 1.41421356237309
#define RSQ2 0.70710678118654744
#define RSQ3 0.57735026918962584
#define RSQ6 0.40824829046386304
#define pi 3.14159265358979323846
#define k_boltzmann 1.380649e-23
#define eV_to_J 1.60217662e-19
#define MPa_to_Pa 1e6

#define Intn 10 //number of integral points in Eshelby calculation
#define Mtr 8 //number of Multithread

#define CRIT_SHAPE 25 //critical shape of ellipsoid

typedef Matrix<int, 6, 1> Vector6i;
typedef Matrix<int, 6, 6> Matrix6i;
typedef Matrix<int, 6, 1> Vector6i;
typedef Matrix<int, 6, 6> Matrix6i;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 10, 1> Vector10d;
typedef Matrix<double,3,Intn*Intn> Integralpoint3;
typedef Matrix<double,6,Intn*Intn> Integralpoint6;
typedef Matrix<double,Intn*Intn,1> Integralpoint1;

