#include "EVPSC.h"
#include "global.h"
#include <string>

TwinG::TwinG() {};

TwinG::TwinG(json &j_twin)
{
    /* logger.info("Initializing twin system..."); */
    num = j_twin["id"];  shear_modulus = j_twin["G"];
    type = mode_type::twin; grain_link = -1;
    /* logger.info("Twin system: " + to_string(j_slip["id"])); */
    /* logger.info("Shear modulus: " + to_string(shear_modulus)); */
    shear_rate = 0; drate_dtau = 0; acc_strain = 0; disloc_velocity = 0, child_frac = 0.0; 

    VectorXd info_vec = to_vector(j_twin, "sn", 6);
    plane_norm = info_vec(seq(0,2)); //normal 
    burgers_vec = info_vec(seq(3,5)); //Burgers
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());

    rate_sen = j_twin["nrsx"];
    int len = j_twin["CRSS_p"].size();
    for (int i = 0; i != len; ++i) harden_params.push_back(j_twin["CRSS_p"][i]); 
    for (int i = 0; i != j_twin["hst"].size(); ++i) latent_params.push_back(j_twin["hst"][i]);
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    crss = harden_params[8];
    if (harden_params[7] != 0 || ~isnan(harden_params[7])) {
        ref_strain_rate = harden_params[7];
    }
    status = twin_status::inactive;
}

TwinG::TwinG(TwinG* t_twin, bool a)
{
    /* logger.info("Initializing twin system..."); */
    num = t_twin->num;  shear_modulus = t_twin->shear_modulus;
    type = mode_type::twin; grain_link = -1; 
    shear_rate = 0; drate_dtau = 0; acc_strain = 0; disloc_velocity = 0, child_frac = 0.0; 

    plane_norm = t_twin->plane_norm;
    burgers_vec = t_twin->burgers_vec;
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());

    rate_sen = t_twin->rate_sen;
    harden_params = t_twin->harden_params;
    latent_params = t_twin->latent_params;
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    crss = harden_params[8];
    if (harden_params[7] != 0 || ~isnan(harden_params[7])) {
        ref_strain_rate = harden_params[7];
    }
    status = twin_status::inactive;
    
}

void TwinG::check_sn_mode()
{
    logger.info("Twin system: " + to_string(num));
    logger.info("Normal: ");
    logger.info(plane_norm.transpose());
    logger.info("Burgers: ");
    logger.info(burgers_vec.transpose());
    logger.info("Pij: ");
    logger.info(Pij);
}

void TwinG::check_hardening_mode()
{
    logger.info("Twin system: " + to_string(num));
    string str = "";
    for (int i = 0; i != harden_params.size(); ++i) str += std::to_string(harden_params[i]) + "\t";
    logger.info("Hardening parameters: " + str);
}

void TwinG::cal_strain_rate(Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);
    rss_matrix = max(abs(rss_matrix) - harden_params[9], 0.0) * sign(rss_matrix);
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen)* sign(rss_matrix) * equivalent_frac; 
    switch (status) {
    inactive:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0; // rate_TN;
        break;
    growth:
        shear_rate = rate_;
        break;
    saturated:
        shear_rate = (rss_matrix < 0) ? rate_ : 0.0; // rate_MP;
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        break;
    }
    rss = rss_matrix;
}

void TwinG::cal_drate_dtau(Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);       
    rss_matrix = max(abs(rss_matrix) - harden_params[9], 1e-5) * sign(rss_matrix);
    double gradient_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen-1) * sign(rss_matrix) / rate_sen / crss * sign(rss_matrix) * equivalent_frac; 
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen) * sign(rss_matrix) * equivalent_frac; 
    switch (status){
    inactive:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0; // rate_TN;
        drate_dtau = (rss_matrix > 0) ? gradient_ : 0.0; // rate_TN;
        break;
    growth:
        shear_rate = rate_;
        drate_dtau = gradient_;
        break;
    saturated:
        shear_rate = (rss_matrix < 0) ? rate_ : 0.0; // rate_MP;
        drate_dtau = (rss_matrix < 0) ? gradient_ : 0.0; // rate_MP;
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        drate_dtau = (rss_matrix > 0) ? gradient_ : 0.0;
        break;
    }
}

void TwinG::update_status(grain &gr, double dtime){
    temperature = gr.temperature;
    double Gamma = 0;
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3], \
           twin_shear = harden_params[4], upper_bound = global_polycrys.twin_threshold;

    for(int i = 0; i < gr.modes_num; i++) Gamma += abs(gr.gmode[i]->acc_strain);
    if (gr.child_frac > upper_bound) status = twin_status::saturated;
    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(int i = 0; i < gr.modes_num; i++)
        crss += abs(gr.gmode[i]->shear_rate) * dtime * gr.lat_hard_mat[num][i] * dtau_by_dGamma;
    equivalent_frac = (1-gr.child_frac) + child_frac;
    euler_twin = get_twin_euler_vec(gr.get_Euler_M_g(), gr.weight_ref*child_frac, plane_norm);
}

void TwinG::update_ssd(Matrix3d strain_rate, double dtime){
    acc_strain += abs(shear_rate) * dtime;
    double lower_bound = 1e-3, upper_bound = global_polycrys.twin_threshold;
    double fraction_incre = shear_rate / harden_params[4] * dtime; 
    if (child_frac + fraction_incre > 1){
        fraction_incre = min(1 - child_frac, 0.5 * child_frac);
    } 
    else if (child_frac + fraction_incre < 0) {
        fraction_incre = sign(fraction_incre) * min(abs(child_frac), 0.5 * child_frac);
    }
    child_frac += fraction_incre;
    disloc_density = child_frac;
    if (status == twin_status::inactive && child_frac >= lower_bound) {
        status = twin_status::growth;
        crss = harden_params[0];
    } 
    if (child_frac < lower_bound) status = twin_status::inactive;
    else if (child_frac > upper_bound) status = twin_status::saturated;
    else status = twin_status::growth;
}

void TwinG::print(){
    logger.info("Twin system: " + to_string(num));
    string status_str = "";
    switch (status){
        case twin_status::growth:
            status_str = "growth";
            break;
        case twin_status::inactive:
            status_str = "inactive";
            break;
        case twin_status::saturated:
            status_str = "saturated";
            break;
        case twin_status::governed:
            status_str = "governed";
            break;
        default:
            status_str = "unknown";
            break;
    }
    logger.info("Twin status: " + status_str);
    logger.info("Burgers vector: ");
    logger.info(burgers_vec.transpose());
    logger.info("Plane normal: ");
    logger.info(plane_norm.transpose());
    logger.info("Dislocation density: " + to_string(disloc_density));
    logger.info("Shear modulus: " + to_string(shear_modulus));
    logger.info("Strain rate: " + to_string(shear_rate));
    logger.info("Accumulated strain: " + to_string(acc_strain));
    logger.info("Critical resolved shear stress: " + to_string(crss));
    logger.info("Child fraction: " + to_string(child_frac));
    logger.info("Harden parameters: ");
    string str = "";
    for (int i = 0; i < harden_params.size(); i++) str += to_string(harden_params[i]) + "  "; 
    logger.info(str);
    logger.info("Latent harden parameters: ");
    string str2 = "";
    for (int i = 0; i < latent_params.size(); i++) str2 += to_string(latent_params[i]) + "  "; 
    logger.info(str2);
}

void TwinG::set_status(twin_status s){
    status = s;
}
