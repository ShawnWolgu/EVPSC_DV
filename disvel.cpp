#include "EVPSC.h"
#include "global.h"
#include <string>

double waiting_time(double rss, double freq_r, double act_energy_r, double forest_stress, double resistance_slip, double energy_expo, double temperature);
double running_time(double rss, double c_drag, double speed_sat, double mean_free_path, double burgers, double forest_stress, double temperature);
vector<double> waiting_time_grad(double rss, double freq_r, double act_energy_r, double forest_stress, double resistance_slip, double energy_expo, double temperature);
vector<double> running_time_grad(double rss, double c_drag, double speed_sat, double mean_free_path, double burgers, double forest_stress, double temperature);
const double vel_threshold = 1e-40;

double Slip::disl_velocity(double rss){
    // * [velocity parameters] 
    //  *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
    //  *  6. saturated speed, 7. drag coefficient
    //  * [hardening parameters] 
    //  *  8. forest hardening coefficient
    //  * [DD evolution parameters] 
    //  *  0. SSD_density, 9. nucleation coefficient, 10. nucleation threshold stress, 11. multiplication coefficient
    //  *  12. drag stress D, 13. reference strain rate, 14. c/g, 15. coplanar reaction coefficient
    //  *
    //  * update parameters:
    //  * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress
    //  */
    double freq_r = harden_params[2], act_energy_r = harden_params[3]/factor(Current_intensity), resistance_slip = harden_params[4],\
           energy_expo = harden_params[5], speed_sat = harden_params[6], c_drag = harden_params[7];
    resistance_slip = resistance_slip/factor(Current_intensity); //the renewed resistence by the current pulsing
    //freq_r = freq_r*factor(Current_intensity); //the renewed freq by the current pulsing
    double burgers = update_params[0], mean_free_path = update_params[1], forest_stress = update_params[3];
    t_wait = waiting_time(rss, freq_r, act_energy_r, forest_stress, resistance_slip, energy_expo, temperature); 
    t_run = running_time(rss, c_drag, speed_sat, mean_free_path, burgers, forest_stress, temperature);
    /* if(abs(rss) - forest_stress < 0.0) return mean_free_path / (t_wait); */
    /* else return mean_free_path / (t_wait + t_run); */
    return mean_free_path / (t_wait+t_run);
}

double waiting_time(double rss, double freq_r, double act_energy_r, double forest_stress, double resistance_slip,\
                    double energy_expo, double temperature){
    double stress_eff = (abs(rss) - forest_stress) * MPa_to_Pa; resistance_slip = resistance_slip*MPa_to_Pa, act_energy_r = act_energy_r*eV_to_J;
    double act_energy = 0.; freq_r = 1e14;
    if (stress_eff >= 0.0) act_energy = act_energy_r * (1-pow((stress_eff/resistance_slip), energy_expo)); // activation energy
    else act_energy = act_energy_r * (1+pow((abs(stress_eff)/resistance_slip), energy_expo)); // activation energy
    /* act_energy = min(act_energy, 500 * k_boltzmann * temperature); // avoid too large activation energy; */
    /* act_energy = max(act_energy, -500 * k_boltzmann * temperature); // avoid too small activation energy; */
    return 1 / freq_r * exp(act_energy / (k_boltzmann * temperature));
}


double running_time(double rss, double c_drag, double speed_sat, double mean_free_path, double burgers, double forest_stress, double temperature){
    double stress_eff = (abs(rss) - forest_stress) * MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * pow(burgers,2));
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + pow(v_norm,-2)) - 1/v_norm);
    velocity = max(velocity, vel_threshold); // avoid velocity = 0 (in case of v_norm = 0)
    return mean_free_path / velocity;
}

vector<double> Slip::disl_velocity_grad(double rss){
    /*
     * [velocity parameters] 
     *  1. MFP control coeffient, 2. reference frequency, 3. activation energy, 4. slip resistance, 5. energy exponent
     *  6. saturated speed, 7. drag coefficient
     * [hardening parameters] 
     *  8. forest hardening coefficient
     * [DD evolution parameters] 
     *  0. SSD_density, 9. multiplication coefficient, 10. drag stressD, 11. reference strain rate, 12. c/g 
     *
     * update parameters:
     * 0: burgers, 1: mean_free_path, 2: disl_density_resist, 3: forest_stress,
     */
    double freq_r = harden_params[2], act_energy_r = harden_params[3]/factor(Current_intensity), resistance_slip = harden_params[4], \
           energy_expo = harden_params[5], speed_sat = harden_params[6], c_drag = harden_params[7];
    double burgers = update_params[0], mean_free_path = update_params[1], forest_stress = update_params[3];
    resistance_slip = resistance_slip/factor(Current_intensity);//the renewed resistence by the current pulsing
    //freq_r = freq_r*(1+pow((Current_intensity/ref_current_intensity),2));//the renewed freq by the current pulsing
    vector<double> dtwait_drss = waiting_time_grad(rss, freq_r, act_energy_r, forest_stress, resistance_slip, energy_expo, temperature);
    vector<double> dtrun_drss = running_time_grad(rss, c_drag, speed_sat, mean_free_path, burgers, forest_stress, temperature);
    double dvel_dtau = 0.0, velocity = 0.0;
    /* if (abs(rss) - forest_stress < 0.0){ */
    /* if (true){ */
        /* dvel_dtau = -mean_free_path / pow((dtwait_drss[1]),2) * (dtwait_drss[0]); */
        /* velocity = mean_free_path / (dtwait_drss[1]); */
    /* } */
    /* else{ */
        dvel_dtau = -mean_free_path / pow((dtwait_drss[1] + dtrun_drss[1]),2) * (dtwait_drss[0] + dtrun_drss[0]);
        velocity = mean_free_path / (dtwait_drss[1] + dtrun_drss[1]);
    /* if (isnan(dvel_dtau)){ */
    /*     logger.debug(to_string(dtwait_drss[0])); */
    /* } */
    /* } */
    vector<double> result = {dvel_dtau,velocity};
    /* logger.debug("rss: " + to_string(rss) + "for_stress: " + to_string(forest_stress) + "tau_eff: " + to_string(abs(rss) - forest_stress) + "dvel_dtau: " + to_string(dvel_dtau) + ", velocity: " + to_string(velocity)); */
    /* logger.debug("reference_vel: " + to_string(mean_free_path * freq_r)); */
    return result;
}

vector<double> waiting_time_grad(double rss, double freq_r, double act_energy_r, double forest_stress, double resistance_slip, double energy_expo, double temperature){
    /* Return a vector: 0. dtw/dtau, 1. tw. */
    double stress_eff = (abs(rss) - forest_stress) * MPa_to_Pa; resistance_slip = resistance_slip*MPa_to_Pa, act_energy_r = act_energy_r*eV_to_J;
    double act_energy = 0.; freq_r = 1e14;
    if (stress_eff >= 0.0) act_energy = act_energy_r * (1-pow((stress_eff/resistance_slip), energy_expo)); // activation energy
    else act_energy = act_energy_r * (1+pow((abs(stress_eff)/resistance_slip), energy_expo)); // activation energy
    double waiting_time = 1 / freq_r * exp(act_energy / (k_boltzmann * temperature));
    if (isinf(waiting_time)){
        waiting_time = 1e100;
    }
    double grad_term = 0;
    if (stress_eff >= 0.0) 
        grad_term = - sign(rss) * act_energy_r * energy_expo * pow((stress_eff/resistance_slip), energy_expo-1) / (k_boltzmann * temperature * resistance_slip);
    else
        grad_term = sign(stress_eff) * sign(rss) * act_energy_r * energy_expo * pow((abs(stress_eff)/resistance_slip), energy_expo-1) / (k_boltzmann * temperature * resistance_slip);
    vector<double> result ={grad_term*waiting_time*MPa_to_Pa, waiting_time };
    return result;
}

vector<double> running_time_grad(double rss, double c_drag, double speed_sat, double mfp, double burgers, double forest_stress, double temperature){
    /* Return a vector: 0. dtr/dtau, 1. tr. */
    double stress_eff = (abs(rss) - forest_stress) * MPa_to_Pa;
    double coeff_B = (c_drag * k_boltzmann * temperature) / (speed_sat * burgers * burgers);
    double v_norm = 2 * burgers * stress_eff / (coeff_B * speed_sat);
    double velocity = speed_sat * (sqrt(1 + 1/pow(v_norm,2))-1/v_norm); 
    velocity = max(velocity, vel_threshold); // avoid velocity = 0 (in case of v_norm = 0)
    double gradient = 0.0;
    if (velocity != vel_threshold) gradient = -1 * sign(rss) * (2*burgers*mfp)/(pow(velocity,2)*coeff_B) * pow(v_norm,-2) * (1-1/(v_norm * sqrt(1+ pow(v_norm,-2))));
    vector<double> result = {gradient*MPa_to_Pa, mfp/velocity};
    return result;
}
