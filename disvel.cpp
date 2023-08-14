#include "Slip.h"

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref, double temperature);
double running_time(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress, double v_c, double temperature);
vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                                 double Peierls_stress, double expo_kinkeng, double temperature_ref, double temperature);
vector<double> running_time_grad(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress, double v_c, double temperature);

double Slip::disl_velocity(double rss){
    /*
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11: v_c, 12: D, 13: ref_srate, 14: gg;
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_resist, 3: back_stress,
     * 4: barrier_distance, 5:c_annih;
     */
    double freq_Debye = harden_params[1], c_length = harden_params[2], kink_energy_ref = harden_params[3],\
           temperature_ref = harden_params[4], Peierls_stress = harden_params[5], expo_kinkeng = harden_params[6],\
           wave_speed = harden_params[7], c_drag = harden_params[8], v_c = harden_params[11];
    double burgers = update_params[0], disl_density_for = update_params[1],\
           back_stress = update_params[3], barrier_distance = update_params[4];
    if(abs(rss)-back_stress > 0.0){
        t_wait = waiting_time(rss, freq_Debye, c_length, burgers, disl_density_for, kink_energy_ref, back_stress,\
                    Peierls_stress, expo_kinkeng, temperature_ref, temperature); 
        t_run = running_time(rss, c_drag, wave_speed, barrier_distance, burgers, back_stress, v_c, temperature);
        return barrier_distance / (t_wait + t_run);
    }
    else{return 0.0;}
}

double waiting_time(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress,\
                    double Peierls_stress, double expo_kinkeng, double temperature_ref, double temperature){
    rss = rss*MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double freq_const = freq_Debye * pow(burgers,2) * c_length / sqrt(disl_density_for) * Peierls_stress / kink_energy_ref;
    double kink_energy = kink_energy_ref * (1-pow((abs(rss)-back_stress)/Peierls_stress,expo_kinkeng));
    kink_energy = min(kink_energy,500 * k_boltzmann * temperature);
    kink_energy = max(kink_energy,-500 * k_boltzmann * temperature);
    double arrh_term = exp(kink_energy/(k_boltzmann*temperature));
    return 1 / freq_const * arrh_term;
}

double running_time(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress, double v_c, double temperature){
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (wave_speed * burgers * burgers);
    double v_norm = 2 * burgers * (abs(rss)-back_stress) / (coeff_B * wave_speed) + v_c/wave_speed;
    double velocity = wave_speed * (sqrt(1 + 1/pow(v_norm,2))-1/v_norm);
    return barrier_distance / velocity;
}

vector<double> Slip::disl_velocity_grad(double rss, double crss, vector<double> harden_params, vector<double> update_params){
    /*
     * Also calculate the velocity.
     * harden parameters: 0: SSD_density,
     * 1: freq_Debye, 2: c_length, 3: kink_energy_ref, 4: temperature_ref,
     * 5: Peierls_stress, 6: expo_kinkeng, 7: wave_speed, 8: c_drag, 9: c_backstress,
     * 10: c_multi, 11:c_annih, 12:HP_stress.
     * 
     * update parameters:
     * 0: burgers, 1: disl_density_for, 2: disl_density_perp, 3: back_stress,
     * 4: barrier_distance
     */
    double freq_Debye = harden_params[1], c_length = harden_params[2], kink_energy_ref = harden_params[3],\
           temperature_ref = harden_params[4], Peierls_stress = harden_params[5], expo_kinkeng = harden_params[6],\
           wave_speed = harden_params[7], c_drag = harden_params[8], v_c = harden_params[11];
    double burgers = update_params[0], disl_density_for = update_params[1],\
           back_stress = update_params[3], barrier_distance = update_params[4];
    if(abs(rss)-back_stress > 1e-20){
        vector<double> dtwait_drss = waiting_time_grad(rss, freq_Debye, c_length, burgers, disl_density_for, kink_energy_ref, crss,\
                    Peierls_stress, expo_kinkeng, temperature_ref, temperature);
        vector<double> dtrun_drss = running_time_grad(rss, c_drag, wave_speed, barrier_distance, burgers, back_stress, v_c, temperature);
	double dvel_dtau = -barrier_distance / pow((dtwait_drss[1] + dtrun_drss[1]),2) * (dtwait_drss[0] + dtrun_drss[0]);
        double velocity = barrier_distance / (dtwait_drss[1] + dtrun_drss[1]);
        vector<double> result = {dvel_dtau,velocity};
        return result;
    }
    else{
        vector<double> result = {0.00,0.00};
        return result;
        }
}

vector<double> waiting_time_grad(double rss, double freq_Debye, double c_length, double burgers, double disl_density_for, double kink_energy_ref, double back_stress, double Peierls_stress, double expo_kinkeng, double temperature_ref, double temperature){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa, Peierls_stress = Peierls_stress*MPa_to_Pa, kink_energy_ref = kink_energy_ref*eV_to_J;
    double freq_const = freq_Debye * pow(burgers,2) * c_length / sqrt(disl_density_for) * Peierls_stress / kink_energy_ref;
    //double freq_const = 1e13;
    double kink_energy = kink_energy_ref * (1-pow((abs(rss)-back_stress)/Peierls_stress,expo_kinkeng));
    kink_energy = min(kink_energy,500 * k_boltzmann * temperature);
    double arrh_term = exp(kink_energy/(k_boltzmann*temperature));
    double waiting_time = 1 / freq_const * arrh_term;
    double grad_const = -1 * sign(rss) * expo_kinkeng * kink_energy_ref /(Peierls_stress*k_boltzmann*temperature)/freq_const;
    double exp_term = 0;
    if (kink_energy != 500*k_boltzmann*temperature && kink_energy != -500*k_boltzmann*temperature) exp_term = arrh_term * pow((abs(rss)-back_stress)/Peierls_stress,expo_kinkeng-1);
    vector<double> result ={ grad_const*exp_term*MPa_to_Pa, waiting_time};
    return result;
}

vector<double> running_time_grad(double rss, double c_drag, double wave_speed, double barrier_distance, double burgers, double back_stress, double v_c, double temperature){
    /* Return a vector: 0. dtr/dtau, 1. tr. */
    rss = rss* MPa_to_Pa, back_stress = back_stress*MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (wave_speed * burgers * burgers);
    double v_norm = 2 * burgers * (abs(rss)-back_stress) / (coeff_B * wave_speed) + v_c/wave_speed;
    double velocity = wave_speed * (sqrt(1 + 1/pow(v_norm,2))-1/v_norm); 
    double gradient = -1 * sign(rss) * (2*burgers*barrier_distance)/(pow(velocity,2)*coeff_B) / pow(v_norm,2) * (1-1/(v_norm * sqrt(1+1/pow(v_norm,2))));
    vector<double> result = {gradient*MPa_to_Pa, barrier_distance/max(velocity,1e-40)};
    return result;
}
