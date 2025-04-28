#pragma once

#include "common/base.h"

//[Classes]
struct ieMode;
struct modeSys;
//
class grain;
enum class mode_type {slip, twin, undefined};
enum class twin_status {inactive, growth, saturated, governed};

//[Class realization]
class PMode
    {
    protected:
        double temperature = 0.; // need to be synchronized with the grain;
        double ref_strain_rate = 0.001;
        Matrix3d Rij;
        double shear_rate_old, drate_dtau_old, disloc_density_old, crss_old, acc_strain_old, rss_old, velocity_old, rho_init_old, rho_H_old, rho_debri_old;

    public:
        PMode();
        PMode(const ieMode &);
        PMode(PMode* t_mode, bool a);
        PMode(const PMode& other);
        PMode(PMode&& other) noexcept;
        PMode& operator=(const PMode& other);
        PMode& operator=(PMode&& other) noexcept;
        virtual ~PMode() = default;
        virtual PMode* clone() const;
        Matrix3d Pij;
        Matrix3d JPij;//J schmid 
        mode_type type = mode_type::undefined;
        int num = -1;
        Vector3d burgers_vec, plane_norm;
        void cal_shear_modulus(Matrix6d elastic_modulus);
        /* 
         * [Slip parameters : disvel model]
         *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent 
         *  6. saturated speed, 7. drag coefficient
         *
         * [hardening parameters] 
         *  8. forest hardening coefficient
         *
         * [DD evolution parameters] 
         *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
         *  12. drag stress D, 13. reference strain rate, 14. c/g, 15. coplanar reaction coefficient, 
         *  16. HP stress, 17. debris_control_param
         *
         * update_params: 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
         *
         * [Slip parameters : Voce model]
         *  0. tau_0, 1. tau_1, 2. h_0, 3. h_1
         *
         * [Twin parameters]
         * 0. tau_0, 1. tau_1, 2. h_0, 3. h_1, 4. twin_strain, 5. A1 6. A2, 7. ref_rate, 8. tau_nuc, 9. tau_HP
         */
        vector<double> harden_params, update_params, latent_params;
        double rate_sen, shear_rate, drate_dtau, shear_modulus, disloc_density, crss, acc_strain, rss = 0.0, velocity = 0.0;
        double J_slipsystem = 0.0;
        double rho_init = 0.0, rho_H=0.0, rho_debri=0.0;
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
        Matrix6d get_M_secant(Matrix3d);
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
        Slip(const ieMode &);
        Slip(Slip* t_mode, bool a);
        ~Slip() override = default;
        Slip(const Slip& other);
        Slip& operator=(const Slip& other);
        Slip(Slip&& other) noexcept;
        Slip& operator=(Slip&& other) noexcept;
        Slip* clone() const override;

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
        Twin(const ieMode &);
        Twin(Twin* t_mode, bool a);

        ~Twin() override;
        Twin(const Twin& other);
        Twin& operator=(const Twin& other);
        Twin(Twin&& other) noexcept;
        Twin& operator=(Twin&& other) noexcept;
        Twin* clone() const override;

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
        TwinG(const ieMode &);
        TwinG(TwinG* t_mode, bool a);

        ~TwinG() override;
        TwinG(const TwinG& other);
        TwinG& operator=(const TwinG& other);
        TwinG(TwinG&& other) noexcept;
        TwinG& operator=(TwinG&& other) noexcept;
        TwinG* clone() const override;

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

