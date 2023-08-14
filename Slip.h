#ifndef MODES_H
#define MODES_H

#include "Toolbox.h"
#include "func.h"
#include "Grains.h"
#include <vector>

using namespace std;
class grain;

//machnical system from which slip&twinning system can be derived 
class Slip
    {
    protected:
        double ref_rate = 1e-3; 
        double temperature = 273; // need to be synchronized with the grain;
        double crss_incr = 0;
        Matrix3d Pij;
        Matrix3d Rij;

    private:
        int mtype = 1; //type of deformation modes (1: slip; 0: twin)
        double gamma_rate_abs_m = 0; //shear strain of one mode
        //New Contents 
        void cal_strain_rate(Matrix3d stress_tensor);
        void cal_strain_rate_pow(Matrix3d stress_tensor);
        void cal_strain_rate_disvel(Matrix3d stress_tensor);
        void cal_drate_dtau(Matrix3d stress_tensor);
        void cal_drate_dtau_pow(Matrix3d stress_tensor);
        void cal_drate_dtau_disvel(Matrix3d stress_tensor);
        void update_voce(Slip* slip_sys, MatrixXd lat_hard_mat, int nmode, double dtime);
        void update_disvel(Slip* slip_sys, MatrixXd lat_hard_mat, int bv, double nmode ,double dtime);

    public:
        Slip();
        Slip(json &j_slip);
        int num = -1, flag_harden;
        Vector3d burgers_vec,plane_norm;
        vector<double> harden_params, update_params, latent_params;
        double ref_strain_rate = 0.001, rate_sen, strain_rate_slip, drate_dtau, shear_modulus, disloc_density, crss, acc_strain, disloc_velocity, rho_sat = 0.0, lh_coeff = 1.0, rho_mov = 0.0;
        double t_wait = 0.0, t_run = 0.0;
        // Original Contents
        int ini_sn_mode(VectorXd matrix_input, int mode_type, int system_n);
        int check_sn_mode();
        int ini_hardening_mode(double nrsx_in, VectorXd hardens_in, VectorXd latents_in);
        int check_hardening_mode();
        double get_gamma0();
        double get_nrsx();
        double cal_rss(Matrix3d stress_tensor);
        double cal_relative_rss(Matrix3d stress_tensor);
        double update_shear_strain_m();
        Matrix3d cal_dijpmode(Matrix3d);
        Matrix3d cal_rotslip_m();
        Matrix6d get_Fgradm(Matrix3d);
        // New Contents
        void update_status(grain &grain, double dtime); //update the status of slip/twinning system
        void update_ssd(Matrix3d strain_rate, double dtime);
        void update_lhparams(Matrix3d strain_rate);
        double disl_velocity(double rss);
        vector<double> disl_velocity_grad(double rss, double crss, vector<double> harden_params, vector<double> update_params);
        void cal_shear_modulus(Matrix6d elastic_modulus);
        void print();
    };

#endif
