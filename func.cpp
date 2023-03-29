#include "func.h"
#include <Eigen/Eigenvalues>

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M);
Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M);

Vector6d tensor_trans_order(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 1*6 Vector: [00,11,22,12,02,01] 
     */
    return Vector6d{{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1)}};
}

Matrix<double,9,1> tensor_trans_order_9(Matrix3d tensor){
    /** 
     * 3*3 Matrix --> 9*1 Vector: [00,11,22,12,02,01,21,20,10] 
     */
    return Matrix<double,9,1> {{tensor(0,0), tensor(1,1), tensor(2,2), tensor(1,2), tensor(0,2), tensor(0,1), tensor(2,1), tensor(2,0), tensor(1,0)}};
}

Matrix<double,9,1> vel_to_dw(Matrix3d tensor){
    // This function is actually invalid;
    Matrix<double,9,9> vel_to_dw_matrix = Matrix<double,9,9>::Zero();
    return vel_to_dw_matrix * tensor_trans_order_9(tensor);
}

Matrix3d tensor_trans_order(Vector6d tensor){
    /** 
     * 1*6 Vector --> 3*3 Matrix: [[0,5,4],[5,1,3],[4,3,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(5), tensor(1), tensor(3)}, {tensor(4), tensor(3), tensor(2)}};
}

Matrix3d tensor_trans_order_9(Matrix<double,9,1> tensor){
    /** 
     * 9*1 Vector --> 3*3 Matrix: [[0,5,4],[8,1,3],[7,6,2]] 
     */
    return Matrix3d{{tensor(0), tensor(5), tensor(4)}, {tensor(8), tensor(1), tensor(3)}, {tensor(7), tensor(6), tensor(2)}};
}

void params_convert_to_matrix(Matrix<double, 15, 1> &params, Vector6d &unknown_params, vector<int> &unknown_idx, Matrix3d &vel_grad_elas, Matrix3d &stress_incr){
    int unk_par_idx = 0;
    for (auto &unk_idx : unknown_idx) {
	params(unk_idx) = unknown_params(unk_par_idx);
	++unk_par_idx;
    }
    Matrix<double,9,1> temp_vel_grad_elas = params(Eigen::seq(0,8));
    Vector6d temp_stress_incr = params(Eigen::seq(9,14));
    vel_grad_elas = tensor_trans_order_9(temp_vel_grad_elas);
    stress_incr = tensor_trans_order(temp_stress_incr);
}

Matrix3d calc_stress(Matrix3d strain_elastic, Matrix6d elastic_modulus){
    Vector6d eps_elastic, sig_vector;
    Matrix3d stress_tensor = Matrix3d::Identity();
    Matrix6d strain_modi_tensor = Matrix6d::Zero();
    strain_modi_tensor.diagonal() << 1,1,1,2,2,2;
    eps_elastic = tensor_trans_order(strain_elastic);
    sig_vector = elastic_modulus * strain_modi_tensor * eps_elastic;
    stress_tensor = tensor_trans_order(sig_vector);
    return stress_tensor;
}

Vector6d get_vec_only_ith(Vector6d &vector_base, int i){
    Vector6d vec_ith = Vector6d::Zero();
    vec_ith(i) = vector_base(i) + 1e-2;
    return vec_ith;
}

void flag_to_idx(Matrix<double, 15, 1> flag, vector<int> &known_idx, vector<int> &unknown_idx){
    int unknown_cur = 0, known_cur = 0;
    for (int idx = 0;idx != 15; ++idx) {
        if(flag(idx) == 0) {unknown_idx[unknown_cur] = idx; ++unknown_cur;}
        else {known_idx[known_cur] = idx; ++known_cur;}
    }
    if (unknown_cur != 6 || known_cur != 9) throw "ERROR: wrong number of unknown parameters in L or dsigma!";
}

Matrix6d rotate_6d_stiff_modu(Matrix6d modulus, Matrix3d rotate_matrix){
    Matrix6d M66 = cal_rotation_trans_6d_for_stiff(rotate_matrix);
    return M66*modulus*M66.transpose();
}

Matrix6d cal_rotation_trans_6d_for_stiff(Matrix3d M){
   Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = 2*xy*xz;
    M66(0,4) = 2*xx*xz;
    M66(0,5) = 2*xx*xy;

    M66(1,3) = 2*yy*yz;
    M66(1,4) = 2*yx*yz;
    M66(1,5) = 2*yx*yy;

    M66(2,3) = 2*zy*zz;
    M66(2,4) = 2*zx*zz;
    M66(2,5) = 2*zx*zy;

    M66(3,0) = yx*zx;
    M66(3,1) = yy*zy;
    M66(3,2) = yz*zz;

    M66(4,0) = xx*zx;
    M66(4,1) = xy*zy;
    M66(4,2) = xz*zz;

    M66(5,0) = xx*yx;
    M66(5,1) = xy*yy;
    M66(5,2) = xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

Matrix6d rotate_6d_compl_modu(Matrix6d modulus, Matrix3d rotate_matrix)
{
    Matrix6d M66 = cal_rotation_trans_6d_for_compl(rotate_matrix);
    return M66*modulus*M66.transpose();
}

Matrix6d cal_rotation_trans_6d_for_compl(Matrix3d M){
	Matrix6d M66;
    double xx,yy,zz,yz,xz,xy,yx,zx,zy;
    xx = M(0,0); 
    yy = M(1,1);
    zz = M(2,2);
    yz = M(1,2);
    xz = M(0,2);
    xy = M(0,1);
	yx = M(1,0);
	zx = M(2,0);
	zy = M(2,1);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
    M66(i,j) = M(i,j)*M(i,j);

    M66(0,3) = xy*xz;
    M66(0,4) = xx*xz;
    M66(0,5) = xx*xy;

    M66(1,3) = yy*yz;
    M66(1,4) = yx*yz;
    M66(1,5) = yx*yy;

    M66(2,3) = zy*zz;
    M66(2,4) = zx*zz;
    M66(2,5) = zx*zy;

    M66(3,0) = 2*yx*zx;
    M66(3,1) = 2*yy*zy;
    M66(3,2) = 2*yz*zz;

    M66(4,0) = 2*xx*zx;
    M66(4,1) = 2*xy*zy;
    M66(4,2) = 2*xz*zz;

    M66(5,0) = 2*xx*yx;
    M66(5,1) = 2*xy*yy;
    M66(5,2) = 2*xz*yz;

    M66(3,3) = yy*zz+yz*zy;
    M66(3,4) = yx*zz+yz*zx;
    M66(3,5) = yx*zy+yy*zx;

    M66(4,3) = xy*zz+xz*zy;
    M66(4,4) = xx*zz+xz*zx;
    M66(4,5) = xx*zy+xy*zx;

    M66(5,3) = xy*yz+xz*yy;
    M66(5,4) = xx*yz+xz*yx;
    M66(5,5) = xx*yy+xy*yx;

	return M66;
}

Matrix3d tensor_rot_to_CryCoord(Matrix3d tensor, Matrix3d orientation){
    return orientation * tensor * orientation.transpose();
}

Matrix3d tensor_rot_to_RefCoord(Matrix3d tensor, Matrix3d orientation){
    return orientation.transpose() * tensor * orientation;
}

int sign(double x){
    return (0 < x) - (x < 0);
}

int heaviside(double x){
    return (x > 0);
}

Matrix3d vel_bc_to_vel_grad(Matrix3d vel_bc_tensor){
    // this function is actually invalid
    Matrix<double,9,9> bc_modi_matrix; // need to be modified.
    Matrix<double,9,1> d_w_vector = bc_modi_matrix * tensor_trans_order_9(vel_bc_tensor);
    Matrix3d d_tensor = tensor_trans_order((Vector6d)d_w_vector(Eigen::seq(0,5)));
    Matrix3d w_tensor;
    w_tensor << 0,d_w_vector(8),d_w_vector(7),-d_w_vector(8),0,d_w_vector(6),-d_w_vector(7),-d_w_vector(6),0;
    return d_tensor + w_tensor;
}

double set_precision(double num, int prec){
    if (num == 0) return 0.0;
    double expo = floor(log10(abs(num)));
    double result = sign(num) * round(abs(num) * pow(10,-expo + prec -1))/pow(10,-expo + prec -1);
    return result;
}

double calc_equivalent_value(Matrix3d mat){
    Matrix3d dev_mat = mat - Matrix3d::Identity() * mat.trace();
    double sum = (dev_mat.cwiseProduct(dev_mat)).sum();
    return sqrt(2./3. * sum);
}

Vector6d set_precision(Vector6d &num, int prec){
    Vector6d result; result << num;
    for (auto &inum : result){
	if (inum != 0){
    	    double expo = floor(log10(abs(inum)));
	    inum = sign(inum) * round(abs(inum) * pow(10,-expo + prec -1))/pow(10,-expo + prec -1);
	}
    }
    return result;
}

void cut_precision(Matrix3d &mat, int prec){
    double ref_value = mat.cwiseAbs().maxCoeff();
    if (ref_value != 0){
	double expo = floor(log10(ref_value)) - prec;
	mat = round(mat.array() * pow(10,-expo)).matrix() * pow(10,expo);
    }
}

double calc_relative_error(Vector6d &v1, Vector6d &v2){
    Vector6d v_error = (v1 - v2).cwiseAbs();
    Vector6d inv_v1 = (0.5 * v1 + 0.5 * v2);
    for (auto &iv : inv_v1) if (iv != 0) iv = 1/iv;
    double result = (v_error.cwiseProduct(inv_v1)).norm();
    return result;
}

double calc_relative_error(double x, double y){
    if (x == 0 && y == 0) return 0;
    else return abs(x)-abs(y)/abs(x);
}

