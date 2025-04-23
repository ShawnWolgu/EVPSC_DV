#pragma once

#include "common/base.h"

struct modeSys {
    MatrixXd sn_info; // 滑移面法向和滑移方向信息
    int mode_n;           // 模式中的系统数量
    int type;             // 类型：0-滑移，1-孪晶，2-其他
    double strainRateSens;          // 应变率敏感性
    vector<double> control_params; // 临界分切应力参数
    vector<double> latent_params;    // 潜在硬化参数
};

struct ieMode {
    int id;               // 模式ID
    int type;             // 类型：0-滑移，1-孪晶，2-其他
    double strainRateSens;          // 应变率敏感性
    double shearModulus;
    Vector3d plane_norm;
    Vector3d character_vec; 
    vector<double> sn; // 滑移面法向和滑移方向信息
    vector<double> control_params; // 临界分切应力参数
    vector<double> latent_params;    // 潜在硬化参数
};

struct materialPhase {
    // 晶体结构相关
    string phaseName; // 材料名称
    //
    string crystalSymmetry;
    vector<double> latticeConstants;
    int Miller_n;
    MatrixXd Mabc;
    MatrixXd Trans_Miller;
    
    MatrixXd elasticConstants; // Elastic constants
    VectorXd thermalExpansionCoeffs; // Thermal expansion coefficients
    MatrixXd piezoelectricConstants; // Piezoelectric constants
    VectorXd dielectricConstants; // Dielectric constants
    
    // 材料热物理特性
    double rho_material = 0;  // 密度
    double Cp_material = 0;   // 比热容
    double sigma_e_mat = 0;   // 电导率
    double h_ext = 0;         // 外部热传导系数
    double Surface = 0;       // 表面积
    double V_sample = 0;      // 样品体积
    double sigma_k = 0;       // 热导率
    
    // 晶粒尺寸
    double grainSize;
    
    // 变形模式相关
    int family_num;           // 考虑的模式族数量
    vector<int> modes_count_by_family; // 每个族的模式数量
    int modes_num;            // 总模式数量
    
    // mode system
    vector<modeSys> m_systems;
    vector<ieMode> info_modes; // 处理后的每个模式信息
    MatrixXd cal_sn_info(MatrixXd &Min, int system_n);
    void sx_info_postprocess();
    double cal_shear_modulus(ieMode &mode);
    
    void add_trans_miller(string crysym);
    void add_lattice_const(const VectorXd& consts);
};

