#ifndef FUNC
#define FUNC

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <time.h> 
#include <math.h>
using namespace std;
using Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix, Eigen::MatrixXd, Eigen::all, Eigen::last;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

// [constants]
#define pi 3.14159265358979323846
#define k_boltzmann 1.380649e-23
#define eV_to_J 1.60217662e-19
#define MPa_to_Pa 1e6

// [functions]
void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx);
void params_convert_to_matrix(Matrix<double, 15, 1> &params, Vector6d &unknown_params, vector<int> &unknown_idx, Matrix3d &vel_grad_elas, Matrix3d &stress_incr);
void cut_precision(Matrix3d &mat, int prec);
int sign(double x);
int heaviside(double x);
double cal_cosine(Vector3d vec_i, Vector3d vec_j);
double set_precision(double num, int prec);
double calc_relative_error(Vector6d &v1, Vector6d &v2);
double calc_relative_error(double x, double y);
double calc_equivalent_value(Matrix3d mat);
Vector6d tensor_trans_order(Matrix3d tensor);
Vector6d get_vec_only_ith(Vector6d &vector_base, int i); 
Vector6d set_precision(Vector6d &num, int prec);
Matrix3d tensor_trans_order(Vector6d tensor);
Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor);
Matrix3d calc_stress(Matrix3d strain_elastic, Matrix6d elastic_modulus);
Matrix3d vel_bc_to_vel_grad(Matrix3d vel_bc_tensor);
Matrix3d tensor_rot_to_CryCoord(Matrix3d tensor, Matrix3d orientation);
Matrix3d tensor_rot_to_RefCoord(Matrix3d tensor, Matrix3d orientation);
Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix);
Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor);
Matrix<double,9,1> vel_to_dw(Matrix3d tensor);

#endif

