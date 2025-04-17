#include "mechanism/PMode_common.h"

Twin::Twin() {};

Twin::Twin(json &j_twin)
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
    crss = harden_params[0];
    status = twin_status::inactive;
}

Twin::Twin(Twin* t_twin, bool a)
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
    crss = harden_params[0];
    status = twin_status::inactive;
}

void Twin::set_parent(int parent_id){
    grain_link = parent_id;
    status = twin_status::governed;
    link_variant = (Twin*) global_polycrys.g[parent_id].gmode[num];
    child_frac = 1;
}

void Twin::check_sn_mode()
{
    logger.info("Twin system: " + to_string(num));
    logger.info("Normal: ");
    logger.info(plane_norm.transpose());
    logger.info("Burgers: ");
    logger.info(burgers_vec.transpose());
    logger.info("Pij: ");
    logger.info(Pij);
}

void Twin::check_hardening_mode()
{
    logger.info("Twin system: " + to_string(num));
    string str = "";
    for (int i = 0; i != harden_params.size(); ++i) str += std::to_string(harden_params[i]) + "\t";
    logger.info("Hardening parameters: " + str);
}

void Twin::cal_strain_rate(Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen)* sign(rss_matrix); 
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
    governed:
        if (link_variant->status == twin_status::growth) shear_rate = rate_;
        else if (link_variant->status == twin_status::saturated) shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        else if (link_variant->status == twin_status::inactive) shear_rate = (rss_matrix < 0) ? rate_ : 0.0;
        else shear_rate = 0.0;
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        break;
    }
    rss = rss_matrix;
}

void Twin::cal_drate_dtau(Matrix3d stress_tensor){
    double rss_matrix = cal_rss(stress_tensor);       
    double gradient_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen-1) * sign(rss_matrix) / rate_sen / crss * sign(rss_matrix); 
    double rate_ = ref_strain_rate * pow(abs(rss_matrix / crss), 1/rate_sen) * sign(rss_matrix); 
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
    governed:
        if (link_variant->status == twin_status::growth) {
            shear_rate = rate_; drate_dtau = gradient_;
        }
        else if (link_variant->status == twin_status::saturated) {
            shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
            drate_dtau = (rss_matrix > 0) ? gradient_ : 0.0;
        }
        else if (link_variant->status == twin_status::inactive) {
            shear_rate = (rss_matrix < 0) ? rate_ : 0.0;
            drate_dtau = (rss_matrix < 0) ? gradient_ : 0.0;
        }
        else {
            shear_rate = 0.0;
            drate_dtau = 0.0;
        }
        break;
    default:
        shear_rate = (rss_matrix > 0) ? rate_ : 0.0;
        drate_dtau = (rss_matrix > 0) ? gradient_ : 0.0;
        break;
    }
}

void Twin::update_status(grain &gr, double dtime){
    temperature = gr.temperature;
    double Gamma = 0;
    double tau_0 = harden_params[0], tau_1 = harden_params[1], h_0 = harden_params[2], h_1 = harden_params[3], \
           twin_shear = harden_params[4], lower_bound = 1e-3, upper_bound = global_polycrys.twin_threshold;
    for(int i = 0; i < gr.modes_num; i++){
        Gamma += abs(gr.gmode[i]->acc_strain);
    }

    double eff_fraction = gr.get_weight_g();
    switch(status){
        case twin_status::growth:{
            if (grain_link != -1) break;
            Vector4d euler_vec = get_twin_euler_vec(gr.get_Euler_M_g(), gr.weight_ref, plane_norm);
            grain_link = global_polycrys.add_grain(euler_vec, gr.grain_i, num);
            (global_polycrys.g[grain_link]).set_weight_g(0.0);
        }
        case twin_status::inactive:
        case twin_status::saturated:{
            if (grain_link == -1) break;
            grain* g_link = &(global_polycrys.g[grain_link]);
            eff_fraction += g_link->get_weight_g();
            eff_fraction = eff_fraction / gr.weight_ref;
            break;
        }
        case twin_status::governed:{
            grain* g_link = &(global_polycrys.g[grain_link]);
            shear_rate = -1 * g_link->gmode[num]->shear_rate;
            eff_fraction += g_link->get_weight_g();
            eff_fraction = eff_fraction / gr.weight_ref;
            break;
        }
        default:
            break;
    }

    double dtau_by_dGamma = h_1 + (abs(h_0/tau_1)*tau_1 - h_1) * exp(-Gamma*abs(h_0/tau_1)) + abs(h_0/tau_1)*h_1*Gamma*exp(-Gamma*abs(h_0/tau_1));
    for(int i = 0; i < gr.modes_num; i++){
        double irate = gr.gmode[i]->shear_rate;
        double lt_hard = gr.lat_hard_mat[num][i];
        crss += abs(irate) * dtime * lt_hard * dtau_by_dGamma;
    }

    if (status == twin_status::governed) crss = link_variant->crss;

    child_frac += eff_fraction * shear_rate / twin_shear * dtime;
    disloc_density = child_frac;
    logger.debug("Twin " + to_string(num) + " child fraction: " + to_string(child_frac));
    
    if (status == twin_status::governed) return;
    if (child_frac < lower_bound) status = twin_status::inactive;
    else if (child_frac > upper_bound) status = twin_status::saturated;
    else status = twin_status::growth;
}

void Twin::update_ssd(Matrix3d strain_rate, double dtime){
    acc_strain += shear_rate * dtime;
}

void Twin::print(){
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

void Twin::set_status(twin_status s){
    status = s;
}
