#pragma once

#include "common/base.h"
#include <locale>

class PMode;
struct materialPhase;

//[Class realization]
class grain
{
    private:
        Matrix6d Cij6; //Elastic stiffness in grain Axes
        double   Cijkl[3][3][3][3]; // Elastic modulus in 4th tensor.
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
        Vector5d d0_g_old; //d0 in last step
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
        int phase_id; // The Phase Number
        int modes_num = 0;
        int if_stress = 0; //flag of stress calculation
        double size = 0.0; //the size of the grain
        double child_frac = 0.0; // the fraction of the child grain
        double weight_ref = 0.0; // the reference weight of the grain
        double temperature = 0., temp_old = 0.;
        bool twin_term_flag;
        materialPhase* mat; // pointer to the materialPhase

        grain();
        grain(const grain& g);
        grain(grain&& other) noexcept;
        grain& operator=(const grain& other);
        grain& operator=(grain&& g) noexcept;
        ~grain();
        //move from private to public
        PMode** gmode = NULL; // deformation modes

        //input the euler angle and weights
        void initialization(int id, int phase_id, Vector4d euler, materialPhase* mat);
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
        int ini_gmode_g(const materialPhase &);
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
        Matrix6d get_Cij6_g();

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
        Matrix5d cal_M_secant(Matrix3d);

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

        // [Visco-plastic consistent Function]
        //
        // Calculate the Viscoplastic Modulus using different linearization method.
        // 1 = Affine; 2 = Tangent; 3 = Secant
        void Update_Mpij6_g(int);
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
        void update_jslip();
    };

