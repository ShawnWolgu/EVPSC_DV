#ifndef EVPSC_H
#define EVPSC_H

#include "Toolbox.h"
#include <vector>

using namespace std;

//[Classes]
class grain;
class PMode;
class Slip;
class Twin;
enum class mode_type {slip, twin, undefined};
enum class twin_status {inactive, growth, saturated, governed};

//[Class realization]
class grain
{
    private:
        Matrix6d Cij6_SA_g; //Elastic stiffness in sample Axes
        Matrix6d Mij6_J_g;  //Elastic compliance invovling Jaumann rate in sample Axes
        Matrix6d Metilde_g; //the (Me~)^-1 in elastic consistency
        Matrix6d Cij6_SA_g_old;
        Matrix6d Mij6_J_g_old;  //Old Elastic compliance invovling Jaumann rate in sample Axes
        Matrix6d Metilde_g_old; //Old the (Me~)^-1 in elastic consistency
        double RSinv_C[3][3][3][3];

        Matrix5d Mptilde_g;// the (M~) in elastic consistency
        Matrix5d Mpij6_g;  // M visco-plastic compliance of grain
        Matrix5d Mptilde_g_old;
        Matrix5d Mpij6_g_old;
        double RSinv_VP[3][3][3][3];

        Vector5d d0_g;

        Matrix3d Fij_g;  //the deformation tensor in grain
        Matrix3d Udot_g; // the velocity gradient in grain
        Matrix3d Dij_g; //the strain rate tensor in grain
        Matrix3d Dije_g; //the elastic strain rate tensor in grain
        Matrix3d Dijp_g; //the vp strain rate tensor in grain
        Matrix3d therm_expansion_g; //the thermal expansion tensor in grain (grain coordinate)
        Matrix3d therm_strain_g; //the thermal strain tensor in grain (global coordinate)
        Matrix3d Wij_g; //the rotation rate tensor in grain
        // Udot_g = Dij_g + Wij_g
        Matrix3d eps_g; //strain of grain
        Matrix3d sig_g; //stress of grain
        // Some variables for restoration
        Matrix3d sig_g_old; //stress in last step
        Matrix3d Dij_g_old; //strain rate in last step
        Matrix3d Dije_g_old; //elastic strain rate in last step
        Matrix3d Dijp_g_old; //vp strain rate in last step
        Matrix3d therm_strain_g_old; //the thermal strain tensor in grain in last step
        double RSinv_C_old[3][3][3][3];
        double RSinv_VP_old[3][3][3][3];
	
        ///////
        //the shape of ellipsoid
        //Vector3d ell_axis_o_g;
        Vector3d ell_axis_g; //axial length of ellipsoid/grain (_o means the original value)
        Matrix3d ell_axisb_g;
        double ell_crit_shape_g = 25;
        Vector3d ellip_ang_g; //the rotate angle of the ellipsoid
        bool Iflat_g = 0; //flag of stretch the ellipsoid (0: yes; 1: no)
        ///////

        //Vector3d euler;  //Euler angles and weight (/degree)
        Matrix3d Euler_M;
        double weight;

        double* gamma_delta_gmode = NULL;

        double gamma_total = 0;
        double gamma_delta = 0; //the increment of gamma

    public:
        int grain_i; // The Number
        int modes_num = 0;
        int if_stress = 0; //flag of stress calculation
        double child_frac = 0.0, weight_ref = 0.0, temperature = 0., temp_old = 0.;
        bool twin_term_flag;

        grain();
        grain(const grain& g);
        grain& operator=(const grain& g);
        //move from private to public
        PMode** gmode = NULL; // deformation modes

        //input the euler angle and weights
        void ini_euler_g(Vector4d);
        Vector3d get_euler_g();
        Vector3d get_euler_g(int mode_num);
        Matrix3d get_Euler_M_g();
        double get_weight_g();
        double get_weight_g_eff();
        double get_weight_g(int mode_num);
        void set_weight_g(double);
        //set thermal expansion tensor
        void therm_expansion_config(Vector6d therm);

        //input the number of deformation modes
        int ini_gmode_g(int);
        int ini_gmode_g(json &);
        int ini_gmode_g(grain &);
        int check_gmode_g();

        //input the normal and Burgers vector in ONE mode (several systems)
        //(MatrixXd) sn, (int) flag of twin or slip, (int) number of systems,(int) mode label
        int ini_sn_g(MatrixXd, int, int, int, Matrix6d);
        int check_sn_g();

        //input the hardening parameters
        //input parameters:
        //double nrsx_in; VectorXd CRSS_p_in; VectorXd hst_in; int modei
        int ini_hardening_g(double, VectorXd, VectorXd, int, int);
        int check_hardening_g();

        //calculate Eshelby tensor in ESC
        //Parameters::
        //double 
        //Vector3d axis_t, Matrix6d C66, //the axis length of ellipsoid
        //Integralpoint6 aa6, Integralpoint6 aaww6, Integralpoint3 alpha // integral points and weights 
        void Eshelby_E(double ESIM[3][3][3][3],double ESCR[3][3][3][3],Vector3d, Matrix6d, Integralpoint6, Integralpoint6, Integralpoint3);

        //calculate Eshelby tensor in VPSC
        //Parameters::
        //Vector3d axis_t, Matrix6d C66, //the axis length of ellipsoid
        //Integralpoint6 aa6, Integralpoint6 aaww6, Integralpoint3 alpha // integral points and weights 
        //Integralpoint3 aww, Integralpoint1 ww
        void Eshelby_P(double ESIM[3][3][3][3],double ESCR[3][3][3][3],Vector3d, Matrix6d, Integralpoint6, Integralpoint6, Integralpoint3, Integralpoint3, Integralpoint1);

        //get the stress of grain
        Matrix3d get_stress_g();
        void set_stress_g(Matrix3d);
        Matrix3d get_strain_g();

        Matrix3d get_Dije_g();
        Matrix3d get_Dijp_g();
        Matrix3d get_Dij_g();
        Matrix3d get_Udot_g();
        Matrix3d get_Wij_g();

        //calculate the stress in grains with Newton-Rapthon iteration
        //Parameters:
        //double Tincr,
        //Matrix3d Wij_m, Matrix3d Dij_AV, Matrix3d Dije_AV, Matrix3d Dijp_AV,
        //Matrix3d Sig_m, Matrix3d Sig_m_old, Matrix3d thermal_expansion_ave
        void grain_stress(double, Matrix3d, Matrix3d, Matrix3d, Matrix3d, Matrix3d, Matrix3d, Matrix3d);

        //calculate the symmetric components
        Matrix3d cal_Dijp(Matrix3d);

        //calculate the skew symmetric components
        Matrix3d cal_rotslip();

        Matrix5d cal_Fgrad(Matrix3d);

        double cal_RSSxmax(Matrix3d); //Calculate the maxinum RSS/CRSS
        double cal_RSSxlim(Matrix3d); //Calculate the limit of RSS/CRSS

        //Elastic consistent
        void save_status_g();
        void restore_status_g();
        void restore_Mij6_J_g();
        void Update_Mij6_J_g(Matrix6d);
        void Update_Cij6_SA_g(Matrix6d);
        void Update_Metilde_g(Matrix6d);
        Matrix6d get_Metilde_g();
        Matrix6d get_Mij6_J_g();

        //return alpha matrix in system coordinate
        Matrix3d get_therm_expansion();
        void Update_RSinv_C_g(double A[3][3][3][3]);

        //Visco-plastic consistent
        void Update_Mpij6_g();
        void Update_Mptilde_g(Matrix5d);
        Matrix5d get_Mpij6_g(); 
        Vector5d get_d0_g();
        void Update_RSinv_VP_g(double A[3][3][3][3]);
        void save_RSinv_g();

        Vector3d get_ell_axis_g();
        Matrix3d get_ell_axisb_g();

        //if the Ishape = 1(in class polycrys)
        // activate these function:

        //according the Fij_g to update the shape;
        void Update_shape_g(); 

        //according the Udot_g to update the Fij_g;
        //parameter:
        //double Tincr: the time increment
        void Update_Fij_g(double);

        //update the accumulate shear strain in all deformation modes
        void Update_shear_strain(double);

        //update the grain orientation
        //parameters:
        //double Tincr, Matrix3d Wij_m
        //Matrix3d Dije_AV, Matrix3d Dijp_AV
        void update_orientation(double, Matrix3d, Matrix3d, Matrix3d);

        void update_strain(double);
        //update the CRSS in the deformation systems
        //parameter:
        //double Tincr
        void update_modes(double);

        //add from SXCpp
        /* MatrixXd lat_hard_mat; */
        /* std::shared_ptr<Eigen::MatrixXd> lat_hard_mat; */
        vector<vector<double>> lat_hard_mat;
        void set_lat_hard_mat();
        void print_latent_matrix();
        // add temperature
        void update_temperature(double Tincr);
    };

class PMode
    {
    protected:
        double temperature = 0.; // need to be synchronized with the grain;
        double ref_strain_rate = 0.001;
        Matrix3d Pij;
        Matrix3d Rij;
        double shear_rate_old, drate_dtau_old, disloc_density_old, crss_old, acc_strain_old, rss_old, velocity_old, rho_init_old, rho_H_old;

    public:
        PMode();
        PMode(json & j_mode);
        PMode(PMode* t_mode, bool a);
        mode_type type = mode_type::undefined;
        int num = -1;
        Vector3d burgers_vec,plane_norm;
        void cal_shear_modulus(Matrix6d elastic_modulus);
        /* 
         * [Slip parameters : disvel model]
         * 0. SSD_density, 1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent 
         * 6. saturated speed, 7. drag coefficient, 8. forest hardening coefficient, 9. nucleation coefficient 
         * 10. multiplication coefficient, 11. drag stress D, 12. reference strain rate, 13. c/g 
         * update_params: 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
         * [Slip parameters : Voce model]
         * * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1
         * [Twin parameters]
         * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. twin_strain, 5. A1 6. A2, 7. ref_rate
         */
        vector<double> harden_params, update_params, latent_params;
        double rate_sen, shear_rate, drate_dtau, shear_modulus, disloc_density, crss, acc_strain, rss = 0.0, velocity = 0.0;
        double rho_init = 0.0, rho_H=0.0;
        // Original Contents
        double get_gamma0();
        double get_nrsx();
        double cal_rss(Matrix3d stress_tensor);
        double cal_relative_rss(Matrix3d stress_tensor);
        double update_shear_strain_m();
        void update_temperature(double temp_in);
        Matrix3d cal_dijpmode(Matrix3d);
        Matrix3d cal_rot_mode();
        Matrix6d get_Fgradm(Matrix3d);
        void save_status();
        void restore_status();
        // Virtual funcs
        virtual void check_hardening_mode(){};
        virtual void check_sn_mode() {};
        virtual void update_status(grain &grain, double dtime) {}; //update the status of slip/twinning system
        virtual void update_ssd(Matrix3d strain_rate, double dtime) {};
        virtual void update_ssd_coplanar_reaction(int modes_num, PMode** sys, double time_incr) {};
        virtual void print();
        virtual void cal_strain_rate(Matrix3d stress_tensor) {};
        virtual void cal_drate_dtau(Matrix3d stress_tensor) {};
    };

class Slip : public PMode
    {
    protected:
        int flag_harden;
        double disloc_velocity = 0.0;

    private:
        double lh_coeff = 1.0, rho_mov = 0.0;
        //New Contents 
        void cal_strain_rate_pow(Matrix3d stress_tensor);
        void cal_strain_rate_disvel(Matrix3d stress_tensor);
        void cal_drate_dtau_pow(Matrix3d stress_tensor);
        void cal_drate_dtau_disvel(Matrix3d stress_tensor);
        void update_voce(PMode** mode_s, vector<vector<double>> lat_hard_mat, int nmode, double dtime);
        void update_disvel(PMode** mode_s, vector<vector<double>> lat_hard_mat, double bv, double nmode ,double dtime);
        double disl_velocity(double rss);
        vector<double> disl_velocity_grad(double rss);

    public:
        Slip();
        Slip(json &j_slip);
        Slip(Slip* t_mode, bool a);
        double t_wait = 0.0, t_run = 0.0, rho_sat = 0.0;
        // Override funcs
        void check_hardening_mode() override;
        void check_sn_mode() override;
        void update_status(grain &grain, double dtime) override; //update the status of slip/twinning system
        void update_ssd(Matrix3d strain_rate, double dtime) override;
        void update_ssd_coplanar_reaction(int modes_num, PMode** sys, double time_incr) override;
        void cal_strain_rate(Matrix3d stress_tensor) override;
        void cal_drate_dtau(Matrix3d stress_tensor) override;
        void print() override;
        // Unused Contents
        void update_lhparams(Matrix3d strain_rate);
    };


class Twin : public PMode
    {
    protected:
        int flag_harden;
        double disloc_velocity = 0.0;

    private:
        twin_status status = twin_status::inactive;
        Twin* link_variant = nullptr;

    public:
        Twin();
        Twin(json &j_twin);
        Twin(Twin* t_mode, bool a);
        int grain_link = -1;
        double t_wait = 0.0, t_run = 0.0, rho_sat = 0.0, child_frac = 0.0;
        void set_parent(int parent_id);
        // Override funcs
        void check_hardening_mode() override;
        void check_sn_mode() override;
        void update_status(grain &grain, double dtime) override; //update the status of slip/twinning system
        void update_ssd(Matrix3d strain_rate, double dtime) override;
        void cal_strain_rate(Matrix3d stress_tensor) override;
        void cal_drate_dtau(Matrix3d stress_tensor) override;
        void print() override;
        void set_status(twin_status s);
        twin_status get_status() { return status; };
    };


class TwinG : public PMode
{
    protected:
        int flag_harden;
        double disloc_velocity = 0.0;

    private:
        twin_status status = twin_status::inactive;
        double equivalent_frac = 0.0;

    public:
        TwinG();
        TwinG(json &j_twin);
        TwinG(TwinG* t_mode, bool a);
        int grain_link = -1;
        double t_wait = 0.0, t_run = 0.0, rho_sat = 0.0, child_frac = 0.0;
        Vector4d euler_twin;
        // Override funcs
        void check_hardening_mode() override;
        void check_sn_mode() override;
        void update_status(grain &grain, double dtime) override; //update the status of slip/twinning system
        void update_ssd(Matrix3d strain_rate, double dtime) override;
        void cal_strain_rate(Matrix3d stress_tensor) override;
        void cal_drate_dtau(Matrix3d stress_tensor) override;
        void print() override;
        void set_status(twin_status s);
        twin_status get_status() { return status; };
    };

#endif

