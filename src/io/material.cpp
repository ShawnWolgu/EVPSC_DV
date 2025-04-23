#include "common/common.h"
#include "io/material.h"

void materialPhase::add_trans_miller(string crysym){
    for (int i = 0; i < crysym.size(); i++) crysym[i] = tolower(crysym[i]);
    vector<double> Mtemp; int n_miller; 
    if(!crysym.compare("hexag"))
        {
            n_miller = 4;
            Mtemp = {1, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 1};
        }
    else if(!crysym.compare("cubic")) 
        {
        n_miller = 3;
        Mtemp = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    }
    else if(!crysym.compare("ortho")) 
        {
        n_miller = 3;
        Mtemp = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    }
    else
    {
        logger.error("Error code 1: crystal symmetry is not supported.");
        throw std::runtime_error("Error code 1: crystal symmetry is not supported.");
        exit(1);
    }
    crystalSymmetry = crysym;
    Miller_n = n_miller;
    Trans_Miller.resize(3, Miller_n);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < Miller_n; ++j) {
            Trans_Miller(i, j) = Mtemp[i * Miller_n + j];
        }
    }
}

void materialPhase::add_lattice_const(const VectorXd& ccon){
    //add lattice constants to sx_json
    vector<double> c_dim, c_angle, Mtemp;
    for(int i = 0; i < ccon.size(); i++){
        if (i < 3) c_dim.push_back(ccon(i));
        else c_angle.push_back(ccon(i)/180*M_PI);
        latticeConstants.push_back(ccon(i));
    }
    Mabc.resize(3, 3);
    //calculate Mabc
    Mabc(0,0)=sin(c_angle[1]);
    Mabc(1,0)=0.;
    Mabc(2,0)=cos(c_angle[1]);
    Mabc(0,1)=(cos(c_angle[2])-cos(c_angle[0])*cos(c_angle[1]))/sin(c_angle[1]);
    Mabc(2,1)=cos(c_angle[0]);
    Mabc(1,1)=sqrt(1.0-pow(Mabc(0,1),2)-pow(Mabc(2,1),2));
    Mabc(0,2)=0.;
    Mabc(1,2)=0.;
    Mabc(2,2)=1.;
    for(int i = 0; i < 3; i++){ 
        for(int j = 0; j < 3; j++) 
            Mabc(i,j) = c_dim[j] * Mabc(i,j);
    }
}

MatrixXd materialPhase::cal_sn_info(MatrixXd &Min, int system_n){
    MatrixXd Min_s, Min_n;
    Min_n = Min(all,seq(0,Miller_n-1)) * Trans_Miller.transpose();
    Min_s = Min(all,seq(Miller_n,2*Miller_n-1)) * Trans_Miller.transpose(); 
    //calculate the coordinate in Cartesian system
    /* logger.debug("Initial Min_n = "); logger.debug(Min_n); */
    MatrixXd Mtemp = Min_n.cwiseAbs().array().max(0.001).matrix();
    Min_n = (Min_n.array().sign() * Mtemp.array()).matrix();
    Min_n = Min_n.cwiseInverse();
    MatrixX3d Min_nv1(Min_n.rows(),3), Min_nv2(Min_n.rows(),3);
    Min_nv1 << -Min_n.col(0), Min_n.col(1), VectorXd::Zero(Min_n.rows());
    Min_nv2 << VectorXd::Zero(Min_n.rows()), -Min_n.col(1), Min_n.col(2);
    Min_nv1 = Min_nv1 * Mabc.transpose();
    Min_nv2 = Min_nv2 * Mabc.transpose();
    /* logger.debug("Min_nv1 = "); logger.debug(Min_nv1); */
    /* logger.debug("Min_nv2 = "); logger.debug(Min_nv2); */
    for (int i = 0; i < Min_n.rows(); ++i) {
        Min_n.row(i) = Min_nv1.row(i).cross(Min_nv2.row(i));
    }
    Min_s = Min_s*Mabc.transpose();
    //normalization
    for(int i = 0; i < system_n; i++)  Min_n.row(i) = Min_n.row(i).normalized();
    /* logger.debug("Min_n = "); logger.debug(Min_n); */
    MatrixXd Min_ns(system_n,6);
    Min_ns.block(0,0,system_n,3) = Min_n;
    Min_ns.block(0,3,system_n,3) = Min_s;
    for(int i = 0; i < system_n; i++)
        for(int j = 0; j < 6; j++)
            if(abs(Min_ns(i,j)) <= 1e-3 ) Min_ns(i,j) = 0.0;
    return Min_ns;
}

void materialPhase::sx_info_postprocess() {
    int mode_id = 0;
    for (auto &isys : m_systems) {
        int system_n = isys.mode_n;
        for (int crt_mode = 0; crt_mode < system_n; crt_mode++) {
            ieMode this_mode;
            this_mode.id = mode_id++;
            this_mode.type = isys.type;
            this_mode.sn.clear();
            for (int i = 0; i < 6; i++){
                this_mode.sn.push_back(isys.sn_info(crt_mode, i));
            }
            this_mode.plane_norm << this_mode.sn[0], this_mode.sn[1], this_mode.sn[2];
            this_mode.character_vec << this_mode.sn[3], this_mode.sn[4], this_mode.sn[5];
            
            this_mode.strainRateSens = isys.strainRateSens;
            this_mode.control_params = isys.control_params;
            this_mode.latent_params = isys.latent_params;
            if (this_mode.type == 0 || this_mode.type == 1){
                this_mode.shearModulus = cal_shear_modulus(this_mode);
            }
            else{
                this_mode.shearModulus = elasticConstants(3,3);
            }
            info_modes.push_back(this_mode);
        }
    }
}

double materialPhase::cal_shear_modulus(ieMode &mode) {
    Matrix3d slip_rotation;
    Vector3d plane_norm, burgers_vec, trav_direc;
    double shear_modulus;
    plane_norm = mode.plane_norm;
    burgers_vec = mode.character_vec;
    trav_direc = burgers_vec.cross(plane_norm);
    slip_rotation << (burgers_vec/burgers_vec.norm()), plane_norm, trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elasticConstants, slip_rotation.transpose())(3,3);
    return shear_modulus;
}
