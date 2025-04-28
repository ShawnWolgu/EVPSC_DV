#include "mechanism/PMode_common.h"

PMode::PMode() {};

PMode::PMode(const ieMode &j_mode){
    /* logger.info("Initializing a PMode object"); */
    type = mode_type::undefined;  
    num = j_mode.id;
    shear_modulus = j_mode.shearModulus;
    shear_rate = 0; drate_dtau = 0; acc_strain = 0, disloc_density = 0;
    plane_norm = j_mode.plane_norm; //normal
    burgers_vec = j_mode.character_vec; //Burgers
    //
    Vector3d sense_vector = plane_norm.cross(burgers_vec/burgers_vec.norm());
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    JPij = 0.5 * (sense_vector * plane_norm.transpose() + plane_norm * sense_vector.transpose());

    rate_sen = j_mode.strainRateSens;
    int len = j_mode.control_params.size();
    for (int i = 0; i != len; ++i) harden_params.push_back(j_mode.control_params[i]); 
    int latent_len = j_mode.latent_params.size();
    for (int i = 0; i != latent_len; ++i) latent_params.push_back(j_mode.latent_params[i]);
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    crss = harden_params[0];
}

PMode::PMode(const PMode& other)
    : temperature(other.temperature),
      ref_strain_rate(other.ref_strain_rate),
      Rij(other.Rij),
      shear_rate_old(other.shear_rate_old),
      drate_dtau_old(other.drate_dtau_old),
      disloc_density_old(other.disloc_density_old),
      crss_old(other.crss_old),
      acc_strain_old(other.acc_strain_old),
      rss_old(other.rss_old),
      velocity_old(other.velocity_old),
      rho_init_old(other.rho_init_old),
      rho_H_old(other.rho_H_old),
      rho_debri_old(other.rho_debri_old),
      Pij(other.Pij),
      JPij(other.JPij),
      type(other.type),
      num(other.num),
      burgers_vec(other.burgers_vec),
      plane_norm(other.plane_norm),
      harden_params(other.harden_params),
      update_params(other.update_params),
      latent_params(other.latent_params),
      rate_sen(other.rate_sen),
      shear_rate(other.shear_rate),
      drate_dtau(other.drate_dtau),
      shear_modulus(other.shear_modulus),
      disloc_density(other.disloc_density),
      crss(other.crss),
      acc_strain(other.acc_strain),
      rss(other.rss),
      velocity(other.velocity),
      J_slipsystem(other.J_slipsystem),
      rho_init(other.rho_init),
      rho_H(other.rho_H),
      rho_debri(other.rho_debri){}

PMode& PMode::operator=(const PMode& other) {
    if (this != &other) {
        temperature = other.temperature;
        ref_strain_rate = other.ref_strain_rate;
        Rij = other.Rij;
        shear_rate_old = other.shear_rate_old;
        drate_dtau_old = other.drate_dtau_old;
        disloc_density_old = other.disloc_density_old;
        crss_old = other.crss_old;
        acc_strain_old = other.acc_strain_old;
        rss_old = other.rss_old;
        velocity_old = other.velocity_old;
        rho_init_old = other.rho_init_old;
        rho_H_old = other.rho_H_old;
        rho_debri_old = other.rho_debri_old;
        
        Pij = other.Pij;
        JPij = other.JPij;
        type = other.type;
        num = other.num;
        burgers_vec = other.burgers_vec;
        plane_norm = other.plane_norm;
        
        harden_params = other.harden_params;
        update_params = other.update_params;
        latent_params = other.latent_params;
        
        rate_sen = other.rate_sen;
        shear_rate = other.shear_rate;
        drate_dtau = other.drate_dtau;
        shear_modulus = other.shear_modulus;
        disloc_density = other.disloc_density;
        crss = other.crss;
        acc_strain = other.acc_strain;
        rss = other.rss;
        velocity = other.velocity;
        J_slipsystem = other.J_slipsystem;
        rho_init = other.rho_init;
        rho_H = other.rho_H;
        rho_debri = other.rho_debri;
    }
    return *this;
}

PMode::PMode(PMode&& other) noexcept
    : temperature(other.temperature),
      ref_strain_rate(other.ref_strain_rate),
      Rij(std::move(other.Rij)),
      shear_rate_old(other.shear_rate_old),
      drate_dtau_old(other.drate_dtau_old),
      disloc_density_old(other.disloc_density_old),
      crss_old(other.crss_old),
      acc_strain_old(other.acc_strain_old),
      rss_old(other.rss_old),
      velocity_old(other.velocity_old),
      rho_init_old(other.rho_init_old),
      rho_H_old(other.rho_H_old),
      rho_debri_old(other.rho_debri_old),
      Pij(std::move(other.Pij)),
      JPij(std::move(other.JPij)),
      type(other.type),
      num(other.num),
      burgers_vec(std::move(other.burgers_vec)),
      plane_norm(std::move(other.plane_norm)),
      harden_params(std::move(other.harden_params)),
      update_params(std::move(other.update_params)),
      latent_params(std::move(other.latent_params)),
      rate_sen(other.rate_sen),
      shear_rate(other.shear_rate),
      drate_dtau(other.drate_dtau),
      shear_modulus(other.shear_modulus),
      disloc_density(other.disloc_density),
      crss(other.crss),
      acc_strain(other.acc_strain),
      rss(other.rss),
      velocity(other.velocity),
      J_slipsystem(other.J_slipsystem),
      rho_init(other.rho_init),
      rho_H(other.rho_H),
      rho_debri(other.rho_debri)
{
    other.num = -1;
    other.type = mode_type::undefined;
    other.temperature = 0.0;
    other.ref_strain_rate = 0.001;
}

PMode& PMode::operator=(PMode&& other) noexcept {
    if (this != &other) {
        temperature = other.temperature;
        ref_strain_rate = other.ref_strain_rate;
        Rij = std::move(other.Rij);
        shear_rate_old = other.shear_rate_old;
        drate_dtau_old = other.drate_dtau_old;
        disloc_density_old = other.disloc_density_old;
        crss_old = other.crss_old;
        acc_strain_old = other.acc_strain_old;
        rss_old = other.rss_old;
        velocity_old = other.velocity_old;
        rho_init_old = other.rho_init_old;
        rho_H_old = other.rho_H_old;
        rho_debri_old = other.rho_debri_old;
        
        Pij = std::move(other.Pij);
        JPij = std::move(other.JPij);
        type = other.type;
        num = other.num;
        burgers_vec = std::move(other.burgers_vec);
        plane_norm = std::move(other.plane_norm);
        
        harden_params = std::move(other.harden_params);
        update_params = std::move(other.update_params);
        latent_params = std::move(other.latent_params);
        
        rate_sen = other.rate_sen;
        shear_rate = other.shear_rate;
        drate_dtau = other.drate_dtau;
        shear_modulus = other.shear_modulus;
        disloc_density = other.disloc_density;
        crss = other.crss;
        acc_strain = other.acc_strain;
        rss = other.rss;
        velocity = other.velocity;
        J_slipsystem = other.J_slipsystem;
        rho_init = other.rho_init;
        rho_H = other.rho_H;
        rho_debri = other.rho_debri;
        
        other.num = -1;
        other.type = mode_type::undefined;
        other.temperature = 0.0;
        other.ref_strain_rate = 0.001;
    }
    return *this;
}

PMode* PMode::clone() const {
    return new PMode(*this);
}

PMode::PMode(PMode* t_mode, bool a){
    type = mode_type::undefined;  num = t_mode->num;  shear_modulus = t_mode->shear_modulus;
    shear_rate = 0; drate_dtau = 0; acc_strain = 0, disloc_density = 0;

    plane_norm = t_mode->plane_norm; //normal
    burgers_vec = t_mode->burgers_vec; //Burgers
    Vector3d sense_vector = plane_norm.cross(burgers_vec/burgers_vec.norm());
    Pij = 0.5 * (burgers_vec / burgers_vec.norm() * plane_norm.transpose() + plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    Rij = 0.5 * (burgers_vec * plane_norm.transpose() / burgers_vec.norm() - plane_norm * burgers_vec.transpose() / burgers_vec.norm());
    JPij = 0.5 * (sense_vector * plane_norm.transpose() + plane_norm * sense_vector.transpose());

    rate_sen = t_mode->rate_sen;
    harden_params = t_mode->harden_params;
    latent_params = t_mode->latent_params;
    for (int i = 0; i != 5; ++i) update_params.push_back(0);
    crss = harden_params[0];
}

double PMode::get_gamma0(){ return ref_strain_rate; }

double PMode::get_nrsx(){ return rate_sen; }

double PMode::cal_rss(Matrix3d stress_tensor){ return (stress_tensor.cwiseProduct(Pij)).sum();}

double PMode::cal_relative_rss(Matrix3d stress_tensor){ 
    return (stress_tensor.cwiseProduct(Pij)).sum() / crss;
}

Matrix3d PMode::cal_dijpmode(Matrix3d sig){ 
    cal_strain_rate(sig);
    return shear_rate * Pij;
}

Matrix3d PMode::cal_rot_mode(){ return shear_rate * Rij;}

Matrix6d PMode::get_Fgradm(Matrix3d stress_grain)
{
    cal_drate_dtau(stress_grain);
    Vector6d Pijv = voigt(Pij);
    return Pijv * Pijv.transpose() * drate_dtau;
}

Matrix6d PMode::get_M_secant(Matrix3d stress_grain)
{
    cal_strain_rate(stress_grain);
    double secant = shear_rate / rss;
    Vector6d Pijv = voigt(Pij);
    return Pijv * Pijv.transpose() * secant;
}

double PMode::update_shear_strain_m() { return abs(shear_rate);}

void PMode::cal_shear_modulus(Matrix6d elastic_modulus){
    Matrix3d rotation_mat;
    Vector3d trav_direc = burgers_vec.cross(plane_norm);
    rotation_mat << (burgers_vec/burgers_vec.norm()), plane_norm/plane_norm.norm(), trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elastic_modulus, rotation_mat.transpose())(3,3);
}

void PMode::update_temperature(double temp_in){
    temperature = temp_in;
}

void PMode::print(){
    logger.info("Undefined mode number: " + to_string(num));
}

void PMode::save_status(){
    shear_rate_old = shear_rate;
    drate_dtau_old = drate_dtau;
    disloc_density_old = disloc_density;
    crss_old = crss;
    acc_strain_old = acc_strain;
    rss_old = rss;
    velocity_old = velocity;
    rho_init_old = rho_init;
    rho_H_old = rho_H;
    rho_debri_old = rho_debri;
}

void PMode::restore_status(){
    shear_rate = shear_rate_old;
    drate_dtau = drate_dtau_old;
    disloc_density = disloc_density_old;
    crss = crss_old;
    acc_strain = acc_strain_old;
    rss = rss_old;
    velocity = velocity_old;
    rho_init = rho_init_old;
    rho_H = rho_H_old;
    rho_debri = rho_debri_old;
}
