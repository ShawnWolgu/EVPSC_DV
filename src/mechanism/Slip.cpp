#include "mechanism/PMode_common.h"

Slip::Slip() {};

Slip::Slip(json &j_slip)
{
    /* logger.info("Initializing slip system..."); */
    num = j_slip["id"];  shear_modulus = j_slip["G"];
    type = mode_type::slip;
    /* logger.info("Slip system: " + to_string(j_slip["id"])); */
    /* logger.info("Shear modulus: " + to_string(shear_modulus)); */
    shear_rate = 0; drate_dtau = 0; acc_strain = 0; disloc_velocity = 0; 

    VectorXd info_vec = to_vector(j_slip, "sn", 6);
    plane_norm = info_vec(seq(0,2)); //normal 
    burgers_vec = info_vec(seq(3,5)); //Burgers
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());

    rate_sen = j_slip["nrsx"];
    int len = j_slip["CRSS_p"].size();
    if (len == 4) flag_harden = 0; else flag_harden = 1;
    for (int i = 0; i != len; ++i) harden_params.push_back(j_slip["CRSS_p"][i]); 
    for (int i = 0; i != j_slip["hst"].size(); ++i) latent_params.push_back(j_slip["hst"][i]);
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    if (flag_harden == 0) crss = harden_params[0];
    else {
        disloc_density = j_slip["CRSS_p"][0];
        rho_mov = disloc_density;
        rho_init = disloc_density;
        rho_H = disloc_density;
        update_params[0] = burgers_vec.norm() * 1e-10;
    }
}

Slip::Slip(Slip* t_slip, bool a)
{
    /* logger.info("Initializing slip system..."); */
    num = t_slip->num;  shear_modulus = t_slip->shear_modulus;
    type = mode_type::slip;
    shear_rate = 0; drate_dtau = 0; acc_strain = 0; disloc_velocity = 0; 
    plane_norm = t_slip->plane_norm;
    burgers_vec = t_slip->burgers_vec;
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());

    rate_sen = t_slip->rate_sen;
    int len = t_slip->harden_params.size();
    if (len == 4) flag_harden = 0; else flag_harden = 1;
    harden_params = t_slip->harden_params;
    latent_params = t_slip->latent_params;
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    if (flag_harden == 0) crss = harden_params[0];
    else {
        disloc_density = harden_params[0];
        rho_mov = disloc_density;
        rho_init = disloc_density;
        rho_H = disloc_density;
        update_params[0] = burgers_vec.norm() * 1e-10;
    }
}

void Slip::check_sn_mode()
{
    logger.info("Slip system: " + to_string(num));
    logger.info("Normal: ");
    logger.info(plane_norm.transpose());
    logger.info("Burgers: ");
    logger.info(burgers_vec.transpose());
    logger.info("Pij: ");
    logger.info(Pij);
}

void Slip::check_hardening_mode()
{
    logger.info("Slip system: " + to_string(num));
    string str = "";
    for (int i = 0; i != harden_params.size(); ++i) str += std::to_string(harden_params[i]) + "\t";
    logger.info("Hardening parameters: " + str);
}

// add function from SXCpp
void Slip::cal_strain_rate(Matrix3d stress_tensor){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Velocity Model;
     * Power model will be used in case 0, while Velocity model used in case 1;
     */
    switch (flag_harden)
    {
        case 0:
            cal_strain_rate_pow(stress_tensor);
            break;
        case 1:
            cal_strain_rate_disvel(stress_tensor);
            break;
        default:
            cal_strain_rate_disvel(stress_tensor);
            break;
    }
}

void Slip::cal_strain_rate_pow(Matrix3d stress_tensor){
    double rss_slip = cal_rss(stress_tensor);       
    if(abs(rss_slip) > 0.5 * crss){
        shear_rate = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen)* sign(rss_slip); 
    }
    rss = rss_slip;
}

// void Slip::cal_strain_rate_disvel(Matrix3d stress_tensor){
//     double burgers = update_params[0];
//     double rss_slip = cal_rss(stress_tensor);
//     //double rss_j = 0.0; //the stress influenced by the 
//     //custom_vars[1] = rss_slip; //original slip system
//     if (J_slipsystem >= 0 and sign(J_slipsystem*rss_slip)>=0){
//         rss_slip = rss_slip+J_slipsystem*K_ew;
//     }else{
//         rss_slip = 0.0;
//     }
//     //rss_slip = rss_slip+J_slipsystem*K_ew;
//     custom_vars[1] = rss_slip; //new rss
//     custom_vars[2] = J_slipsystem*K_ew;//to record the J_slip force
//     velocity = disl_velocity(rss_slip);
//     shear_rate = abs(rho_mov * burgers * velocity) * sign(rss_slip);
//     rss = rss_slip; //用来算形核率的。
// } 1.26

void Slip::cal_strain_rate_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    double rss_j = 0.0; //the stress influenced by the 
    custom_vars[1] = rss_slip;
    velocity = disl_velocity(rss_slip);
    shear_rate = abs(rho_mov * burgers * velocity) * sign(rss_slip);
    rss = rss_slip;
}

void Slip::cal_drate_dtau(Matrix3d stress_tensor){
    /* Select different model for the gradient of shear strain rate by rss, controlled by flag_harden.
     * Note the slip rate will also be update in this function.
     * 0 : Voce Hardening; 1 : Dislocation Velocity Model;
     * Power model will be used in case 0, while Velocity model used in case 1;
     */
    switch (flag_harden)
    {
        case 0:
            cal_drate_dtau_pow(stress_tensor);
            break;
        case 1:
            cal_drate_dtau_disvel(stress_tensor);
            break;
        default:
            cal_drate_dtau_disvel(stress_tensor);
            break;
    }
}

void Slip::cal_drate_dtau_pow(Matrix3d stress_tensor){
    double rss_slip = cal_rss(stress_tensor);       
    if(abs(rss_slip) > 0.5 * crss){
        drate_dtau = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen-1) * sign(rss_slip) / rate_sen / crss * sign(rss_slip); 
        shear_rate = ref_strain_rate * pow(abs(rss_slip / crss), 1/rate_sen) * sign(rss_slip); 
    }
    rss = rss_slip;
}   

void Slip::cal_drate_dtau_disvel(Matrix3d stress_tensor){
    double burgers = update_params[0];
    double rss_slip = cal_rss(stress_tensor);
    // rss_slip = rss_slip+J_slipsystem*K_ew;
    // if (J_slipsystem >= 0 and sign(J_slipsystem*rss_slip)>=0){
    //     rss_slip = rss_slip+J_slipsystem*K_ew;
    // }else{
    //     rss_slip = 0.0;
    // }
    vector<double> dvel_and_vel = disl_velocity_grad(rss_slip);
    drate_dtau = rho_mov * burgers * sign(rss_slip) * dvel_and_vel[0];
    shear_rate = rho_mov * burgers * dvel_and_vel[1] * sign(rss_slip);
    rss = rss_slip;//
    velocity = dvel_and_vel[1];
}

void Slip::update_status(grain &gr_, double dtime){
    /* Select different model for shear strain rate, controlled by flag_harden.
     * 0 : Voce Hardening; 1 : Dislocation Velocity Model;
     */
    //Update Schmidt here. But not realized yet, maybe pack the update into a function.
    Vector3d update_bv = burgers_vec;
    Vector3d update_nv = plane_norm;
    /* Pij = 0.5 * (update_bv / update_bv.norm() * update_nv.transpose() + update_nv * update_bv.transpose()/update_bv.norm()); */
    /* Rij = 0.5 * (update_bv / update_bv.norm() * update_nv.transpose() - update_nv * update_bv.transpose()/update_bv.norm()); */
    temperature = gr_.temperature;

    switch (flag_harden)
    {
        case 0:
            update_voce(gr_.gmode, gr_.lat_hard_mat, gr_.modes_num, dtime);
            break;
        case 1:
            update_disvel(gr_.gmode, gr_.lat_hard_mat, update_bv.norm(), gr_.modes_num, dtime);
            break;
        default:
            update_disvel(gr_.gmode, gr_.lat_hard_mat, update_bv.norm(), gr_.modes_num, dtime);
            break;
    }
}

void Slip::update_voce(PMode** gmode, vector<vector<double>> lat_hard_mat, int modes_num, double dtime){
    /*
     * Update crss and acc_strain.
     */
    double Gamma = 0;
    for(int i = 0; i < modes_num; i++){
        Gamma += abs(gmode[i]->acc_strain);
    }
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3];
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for (int i = 0; i < modes_num; i++){
        crss += abs(gmode[i]->shear_rate) * dtime * lat_hard_mat[num][i] * dtau_by_dGamma;
    }
}

void Slip::update_disvel(PMode** slip_sys, vector<vector<double>> lat_hard_mat, double bv, double nmode ,double dtime){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
     *  12. drag stress D, 13. reference strain rate, 14. c/g, 15. coplanar reaction coefficient, 16. HP stress
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
     */
    double c_mfp = harden_params[1], resistance_slip = harden_params[4], c_forest = harden_params[8], HP_stress = harden_params[16];
    resistance_slip = resistance_slip/factor(Current_intensity, ref_current_intensity_0);//the renewed resistence by the current pulsing
    double burgers, disl_density_for, disl_density_resist, joint_density, debris_density, forest_stress, boundary_resistance, mfp;
    disl_density_for = disl_density_resist = joint_density = debris_density = 0;
    for(int i = 0; i < nmode; i++){
        disl_density_for += slip_sys[i]->disloc_density;
        disl_density_resist += slip_sys[i]->rho_H * (lat_hard_mat[num][i]);
        if (slip_sys[i]->num != num) {
            joint_density += lat_hard_mat[num][i] * sqrt(slip_sys[i]->rho_H - slip_sys[i]->rho_init) * sqrt(rho_H-rho_init);
        }
        debris_density += slip_sys[i]->rho_debri;
    }
    burgers = bv * 1e-10;

    boundary_resistance = HP_stress + c_forest * shear_modulus * burgers * sqrt(debris_density);
    forest_stress = c_forest * shear_modulus * burgers * sqrt(disl_density_resist + 0.707*joint_density) + boundary_resistance;
    crss = forest_stress + resistance_slip;
    mfp = c_mfp / sqrt(disl_density_for);
    update_params[0] = burgers, update_params[1] = mfp, update_params[2] = disl_density_resist, update_params[3] = forest_stress;
    if (num == 0){
        custom_vars[0] = max(custom_vars[0], disl_density_for);
    }
}

void Slip::update_ssd(Matrix3d strain_rate, double dtime){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
     *  12. drag stress D, 13. reference strain rate, 14. c/g, 15. coplanar reaction coefficient, 16. HP stress, 17. debris_control_param
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
     */
    acc_strain += abs(shear_rate) * dtime;
    if (flag_harden == 0){ return; }
    if (flag_harden == 1){ 
        double c_forest = harden_params[8], c_nuc = harden_params[9], tau_nuc = harden_params[10],\
               c_multi = harden_params[11]/factor(Current_intensity, ref_current_intensity_0), c_annih = 0., c_debri = harden_params[17], \
               D = harden_params[12] * 1e6, ref_srate = harden_params[13], gg = c_forest/harden_params[14],\
               burgers = update_params[0], mfp = update_params[1], forest_stress = update_params[3]; 
               //reduced multiplication rate
        disloc_density = rho_H;
        double equi_strain_rate = calc_equivalent_value(strain_rate);
        gg = gg/factor(Current_intensity, ref_current_intensity_2); //reduced normalized energy 
        D = D ;//* factor(Current_intensity, ref_current_intensity_0); // give energy to drag stress;
        double rho_sat_new = c_forest * burgers / gg * (1-k_boltzmann * temperature/D/pow(burgers,3) * log(abs(equi_strain_rate)/ref_srate));
        rho_sat_new = max(pow(1/rho_sat_new,2), 0.5*disloc_density);
        if (rho_sat_new < 10 * rho_sat) rho_sat = rho_sat_new;
        if (rho_sat == 0.0) rho_sat = rho_sat_new;
        /* custom_vars[1] = max(custom_vars[1], rho_sat); */
        double tau_eff = max(abs(rss) - forest_stress, 0.);
        /* custom_vars[2] = max(custom_vars[2], tau_eff); */
        /* custom_vars[3] = max(custom_vars[3], abs(rss)); */
        /* custom_vars[4] = max(custom_vars[4], forest_stress); */
        double term_nuc = c_nuc * max(abs(rss)-tau_nuc,0.) / (shear_modulus * burgers * burgers);
        /* custom_vars[5] = max(custom_vars[5], term_nuc * abs(shear_rate) * dtime); */
        double term_multi = c_multi / mfp; 
        c_annih = (term_multi + term_nuc) / rho_sat;
        double disloc_incre = (term_multi + term_nuc - c_annih * disloc_density) * abs(shear_rate) * dtime;
        double debris_incre = c_annih * disloc_density * abs(shear_rate) * dtime * c_debri;
        debris_incre = debris_incre > rho_sat ? rho_debri * 0.1 : debris_incre;
        rho_debri += debris_incre;
        if (disloc_incre > 2 * rho_sat) {
            disloc_incre = 2 * rho_sat - disloc_density; 
            rho_debri -= debris_incre;
        }
        else if (disloc_incre + disloc_density < 0) disloc_incre = -0.1 * disloc_density; 
        disloc_density += disloc_incre;
        rho_mov = disloc_density;
        if(disloc_density < rho_init) rho_init = disloc_density;
    }
}

void Slip::update_ssd_coplanar_reaction(int modes_num, PMode** mode_sys, double time_incr){
    double d_term_coplanar = 0, coeff_coplanar = harden_params[15];
    vector<int> coplanar_indices;
    for (int g_num = 0; g_num < modes_num; g_num++) {
        if (mode_sys[g_num]->type != mode_type::slip) continue;
        if (mode_sys[g_num]->num == num) continue;
        int inter_mode = get_interaction_mode(burgers_vec, plane_norm, mode_sys[g_num]->burgers_vec, mode_sys[g_num]->plane_norm);
        if (inter_mode != 2) continue;
        coplanar_indices.push_back(g_num);
    }
    for (auto &index_1 : coplanar_indices) {
        for (auto &index_2 : coplanar_indices) {
            if (index_1 >= index_2) continue;
            double plus_term = mode_sys[index_1]->disloc_density * mode_sys[index_1]->velocity * sqrt(mode_sys[index_2]->disloc_density) + \
                               mode_sys[index_2]->disloc_density * mode_sys[index_2]->velocity * sqrt(mode_sys[index_1]->disloc_density);
            d_term_coplanar += plus_term;
        }
        double minus_term = mode_sys[index_1]->disloc_density * mode_sys[index_1]->velocity * sqrt(disloc_density) + \
                            disloc_density * velocity * sqrt(mode_sys[index_1]->disloc_density);
        d_term_coplanar -= minus_term;
    }
    d_term_coplanar = d_term_coplanar * coeff_coplanar * time_incr;
    if (d_term_coplanar + disloc_density < 0) d_term_coplanar = -disloc_density * 0.5;
    else if (d_term_coplanar > 2 * rho_sat) d_term_coplanar = disloc_density * 0.5;
    rho_H = disloc_density + d_term_coplanar;
    rho_init = min(rho_H, rho_init);
}

void Slip::update_lhparams(Matrix3d strain_rate){
    if (flag_harden == 1){ 
        double ref_srate = 1e-3, exp_lh = -0.1;
        lh_coeff = pow(calc_equivalent_value(strain_rate)/ref_srate, exp_lh);
        if (lh_coeff > 2) lh_coeff = 2;
    }
    else{}
}

void Slip::print(){
    logger.info("Slip system: " + to_string(num));
    logger.info("Burgers vector: ");
    logger.info(burgers_vec.transpose());
    logger.info("Plane normal: ");
    logger.info(plane_norm.transpose());
    logger.info("Dislocation density: " + to_string(disloc_density));
    logger.info("Shear modulus: " + to_string(shear_modulus));
    logger.info("Strain rate: " + to_string(shear_rate));
    logger.info("Accumulated strain: " + to_string(acc_strain));
    logger.info("Critical resolved shear stress: " + to_string(crss));
    logger.info("Harden parameters: ");
    string str = "";
    for (int i = 0; i < harden_params.size(); i++) str += to_string(harden_params[i]) + "  "; 
    logger.info(str);
    logger.info("Latent harden parameters: ");
    string str2 = "";
    for (int i = 0; i < latent_params.size(); i++) str2 += to_string(latent_params[i]) + "  "; 
    logger.info(str2);
}
