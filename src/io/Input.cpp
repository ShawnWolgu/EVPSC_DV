#include "common/common.h"
#include "io/Input.h"
#include "io/material.h"

void skipLine(fstream& file, int n) {
    if(n <= 0) return;
    string tp;
    for(int i = 0; i < n && file.good(); ++i) {
        getline(file, tp);
    }
}

int EVPSCinput(string &ftex,string &fsx,string &fload, Procs::Process &Proc)
{
    fstream ininp;
    logger.info("Loading input file EVPSC_CPP.in ...");
    ininp.open("EVPSC_CPP.in",ios::in);
    if (!ininp) {
        logger.error("Error code 0: loading file cannot be opened.");
        return EXIT_FAILURE;
    }
    //read the file path
    string tp;
    skipLine(ininp);
    getline(ininp, ftex);
    skipLine(ininp);
    getline(ininp, fsx); 
    skipLine(ininp);
    getline(ininp, fload); 

    //read the update control
    skipLine(ininp, 3);
    getline(ininp, tp); 
    VectorXd temp1 = getnum(tp, 4);
    Vector4i temp2;
    for(int i=0; i<4; i++) temp2(i) = int(temp1(i));
    set_control_flags(temp2);

    //read output control
    skipLine(ininp, 3);
    getline(ininp, tp); 
    VectorXd temp = getnum(tp, 1);
    Proc.Out_texset(int(temp(0)));

    ininp.close(); 
    return 0;
}

int loadinput(string fname, Procs::Process &Proc)
{
    fstream loadinp(fname, ios::in);
    logger.info("Loading process file " + fname + " ...");
    if (!loadinp) {
        logger.error("Error code 0: process file cannot be opened.");
        return EXIT_FAILURE;
    }

    string tp;
    getline(loadinp, tp); // 1st line: Loading control option
    Proc.load_ctrl(getnum(tp, 4));

    skipLine(loadinp); //Read Velocity Gradient Tensor control params
    Matrix3i IUdot;
    for (int i = 0; i < 3; ++i) {
        getline(loadinp, tp);
        Vector3d temp = getnum(tp, 3);
        IUdot.row(i) = temp.cast<int>();
    }
    Proc.get_IUdot(IUdot);

    skipLine(loadinp); //Read Velocity Gradient Tensor
    Matrix3d Udot;
    for (int i = 0; i < 3; ++i) {
        getline(loadinp, tp);
        Udot.row(i) = getnum(tp, 3);
    }
    Proc.get_Udot(Udot);

    skipLine(loadinp); //Read Stress Tensor control params
    Vector6i ISdot;
    auto readISdotComponents = [&](const vector<int>& indices, int count) {
        getline(loadinp, tp);
        VectorXd temp = getnum(tp, count);
        for (int i = 0; i < count; ++i) {
            ISdot(indices[i]) = int(temp(i));
        }
    };
    readISdotComponents({0,5,4}, 3);  // 第一行分量
    readISdotComponents({1,3}, 2);    // 第二行分量
    readISdotComponents({2}, 1);      // 第三行分量
    Proc.get_ISdot(ISdot);

    skipLine(loadinp); // Read Stress Tensor
    Vector6d Sig_rate;
    auto readSigComponents = [&](const vector<int>& indices, int count) {
        getline(loadinp, tp);
        VectorXd temp = getnum(tp, count);
        for (int i = 0; i < count; ++i) {
            Sig_rate(indices[i]) = temp[i];
        }
    };
    readSigComponents({0,5,4}, 3);
    readSigComponents({1,3}, 2);
    readSigComponents({2}, 1);
    Proc.get_Sdot(voigt(Sig_rate));
    getline(loadinp ,tp);//skip one line
    if (!loadinp.eof()) //if the file ends, return
    {
        if (tp.find("duty") != tp.npos){
            getline(loadinp, tp); VectorXd electric_coeff = getnum(tp, 11);
            duty_ratio_J = electric_coeff(0);
            Amplitude_J = electric_coeff(1);
            Frequency = electric_coeff(2);
            ref_current_intensity_0 = electric_coeff(3);
            ref_current_intensity_1 = electric_coeff(4);
            ref_current_intensity_2 = electric_coeff(5);
            bvalue = electric_coeff(6);
            shock_int = electric_coeff(7);
            shock_fin = electric_coeff(8);
            flag_emode = electric_coeff(9);
            K_ew = electric_coeff(10);
            logger.debug("duty_ratio_J = " + to_string(duty_ratio_J));
            logger.debug("Amplitude_J = " + to_string(Amplitude_J));
            logger.debug("Frequency = " + to_string(Frequency));
            getline(loadinp, tp); //skip one line
            for(int i = 0; i < 3; i++){
                getline(loadinp, tp);//获得电流张量
                J_tensor.row(i) = getnum(tp, 3);
            }
        }
        if (tp.find("temperature") != tp.npos){
            getline(loadinp, tp); 
            VectorXd tempK_coeff = getnum(tp, 2);
            double tempK_rate = tempK_coeff(0);
            double tempK_end = tempK_coeff(1);
            Proc.set_tempK_control(tempK_rate, tempK_end);
        }
    }
    else{
        flag_emode = 0; //Depress current control.
    }
    loadinp.close();
    return EXIT_SUCCESS;
}

int sxinput(string fname, Polycs::polycrystal &pcrys)
{
    fstream sxinp;
    materialPhase mat; 
    
    sxinp.open(fname,ios::in); //open .sx
    if (!sxinp) //checking whether the file is open
    {
        logger.error("Error code 0: .sx file cannot be opened");
        return 1;
    }
    logger.info("Loading sx file " + fname + " ...");
    string tp;     getline(sxinp, tp); //skip first line
    mat.phaseName = tp.substr(0, tp.find(" ")); //phase name
    logger.info("Phase name: " + mat.phaseName);
    string crysym; getline(sxinp, crysym); //crystal symmetry string
    mat.add_trans_miller(crysym.substr(0,5));
    logger.info("Crystal symmetry: " + crysym.substr(0,5));
    
    getline(sxinp, tp);
    VectorXd lattice_const = getnum(tp, 6);
    mat.add_lattice_const(lattice_const); //crystal constants
    logger.info("Lattice constants: ");
    logger.info(mat.latticeConstants);
    
    getline(sxinp, tp);  //skip a line;
    MatrixXd Cij6(6,6);  //Elastic constants;
    for (int i=0; i<6; i++) { 
        getline(sxinp, tp); 
        Cij6.row(i) = getnum(tp, 6);
    }
    mat.elasticConstants = Cij6;
    logger.info("Elastic constants: ");
    logger.info(mat.elasticConstants);

    getline(sxinp, tp);  //skip a line;
    getline(sxinp, tp);  
    VectorXd therm = getnum(tp, 6); //Thermal coefficients
    mat.thermalExpansionCoeffs = therm;
    logger.info("Thermal coefficients: ");
    logger.info(mat.thermalExpansionCoeffs);
    
    getline(sxinp, tp); // this line is for thermal coeffs check or read plasitcity modes
    if (tp.find("rho_material") != tp.npos){
        getline(sxinp, tp); 
        VectorXd thermal_coeff = getnum(tp, 7);
        mat.rho_material = thermal_coeff(0);
        mat.Cp_material = thermal_coeff(1);
        mat.sigma_e_mat = thermal_coeff(2);
        mat.h_ext = thermal_coeff(3);
        mat.Surface = thermal_coeff(4);
        mat.V_sample = thermal_coeff(5);
        mat.sigma_k = thermal_coeff(6);
        getline(sxinp, tp);  //skip a line;        //Start reading slip and twinning modes
        logger.info("rho_material = " + to_string(mat.rho_material));
        logger.info("Cp_material = " + to_string(mat.Cp_material));
        logger.info("sigma_e_mat = " + to_string(mat.sigma_e_mat));
        logger.info("h_ext = " + to_string(mat.h_ext));
        logger.info("Surface = " + to_string(mat.Surface));
        logger.info("V_sample = " + to_string(mat.V_sample));
        logger.info("sigma_k = " + to_string(mat.sigma_k));
    }
    
    getline(sxinp, tp);  int nmodesx = int(getnum(tp, 1)(0)); //total mode number in file
    getline(sxinp, tp);  int nmodes = int(getnum(tp, 1)(0));  //considered in current run
    getline(sxinp, tp);  VectorXd mode_i = getnum(tp, nmodes);  //the index of modes(mode_i)

    //Start reading slip and twinning modes
    int modes_num = 0; 
    mat.modes_count_by_family.clear();
    mat.m_systems.clear();
    
    bool f = 1;  
    for(int imode = 0; imode < nmodesx; ++imode){
        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);  VectorXd mode_info = getnum(tp, 4);
        /* mode_info 0: the serial number 
         * 1: number of deformation systems
         * 2: flag of slip (0 for twin; 1 for slip) 
         * 3: flag of twin (1 for twin; 0 for slip)
        */ 
        MatrixXd nor_dir(int(mode_info(1)), 2*mat.Miller_n); //normal and direction of slip plane
        for (int i = 0; i < int(mode_info(1)); i++) {
            getline(sxinp, tp); 
            nor_dir(i,all) = getnum(tp, 2*mat.Miller_n);
        }
        f = (mode_i.array() == imode+1).any();
        if(f){
            modeSys this_sys;
            MatrixXd sn_matrix = mat.cal_sn_info(nor_dir, int(mode_info(1)));
            this_sys.sn_info.resize(sn_matrix.rows(), sn_matrix.cols());
            for(int i = 0; i < sn_matrix.rows(); i++) {
                this_sys.sn_info.row(i) = sn_matrix.row(i);
            }
            this_sys.mode_n = int(mode_info(1));
            if (mode_info(2) == 1) this_sys.type = 0; //slip
            else if (mode_info(3) == 1) this_sys.type = 1; //twin
            else this_sys.type = 2; //other
            modes_num += int(mode_info(1));
            mat.modes_count_by_family.push_back(int(mode_info(1)));
            mat.m_systems.push_back(this_sys);
        }                
    }
    mat.family_num = nmodes;
    mat.modes_num = modes_num;

    getline(sxinp, tp);  //skip a line;
    getline(sxinp, tp);  int iharden = int(getnum(tp, 1)(0)); //hardening law(Voce=0, DV=1)
    getline(sxinp, tp);  bool irate = bool(getnum(tp, 1)(0)); //"rate sensitive" flag(1: Y; 0: N)
    getline(sxinp, tp);  mat.grainSize = getnum(tp, 1)(0); //grain size: um
    int harden_size;
    if(iharden == 0) harden_size = 4; else harden_size = 18;

    //Read hardening parameters of modes
    double nrsx; vector<double> CRSS_p, hst;
    for(int imode = 0; imode < nmodes; ++imode)
    {
        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);  nrsx = getnum(tp, 1)(0); //rate sensitivity
        getline(sxinp, tp);  //CRSS parameters
        if (mat.m_systems[imode].type == 0) CRSS_p = getnum_vec(tp, harden_size);
        else CRSS_p = getnum_vec(tp, 10);
        getline(sxinp, tp);  //latent hardening parameters
        if (mat.m_systems[imode].type == 0) hst = getnum_vec(tp, 6); //6 types of hardening
        else hst = getnum_vec(tp, 2);
        mat.m_systems[imode].strainRateSens = nrsx;
        mat.m_systems[imode].control_params = CRSS_p;
        mat.m_systems[imode].latent_params = hst;
    }
    mat.sx_info_postprocess();
    sxinp.close(); //close the file object.
    logger.debug("Finish reading the sx file.");
    pcrys.add_phase_from_material(mat);
    return 0;
}

int texinput(string fname, Polycs::polycrystal &pcrys)
{
    fstream texinp;
    texinp.open(fname,ios::in); //open .tex
    if (texinp.is_open())
        {   //checking whether the file is open
            logger.info("Reading texture file");
            string tp;
            //skip 3 lines;
            for(int i = 0; i < 3; i++)
            {
                getline(texinp, tp);
            }
            //number of grains 
            getline(texinp, tp);
            VectorXd Gn = getnum(tp, 1);
            //
            pcrys.grains_n(int(Gn(0)));
            // 
            //Euler angle and weighs
            Vector4d Euler_w(0,0,0,0);
            for(int i = 0; i < int(Gn(0)); i++)
            {
                getline(texinp, tp);
                Euler_w = getnum(tp, 4);
                pcrys.ini_euler(getnum(tp, 4),i);
            }
            texinp.close(); //close the file object.
            pcrys.Norm_weight();
            return 0;        
        }
    else
    {
        logger.error("Error code 0: texture file cannot be opened");
        return 1;
    }
}

VectorXd getnum(string strin, int num)
{
    int i = 0;
    VectorXd Vtemp(num);
    string pattern("[+-]?[\\d]+([\\.][\\d]*)?([Ee][+-]?[\\d]+)?");
    regex r(pattern);
    smatch results;

    string::const_iterator iter_begin = strin.cbegin();
    string::const_iterator iter_end = strin.cend();
    while (regex_search(iter_begin, iter_end,  results,  r))
    {
        if (i >= num) break;
        Vtemp(i)=stod(results[0].str());
        iter_begin = results[0].second;	
        i++;
    }
    if (i < num) {
        logger.error("Error code 1: getnum() cannot find enough numbers, expected " + to_string(num) + " numbers, but only found " + to_string(i) + " numbers.");
        logger.error("The involved line is " + strin);
        exit(1);
    }
    return Vtemp;
}

vector<double> getnum_vec(string strin, int num){
    int i = 0;
    vector<double> Vtemp;
    //string pattern("\\d+(\\.\\d+)?");
    string pattern("[+-]?[\\d]+([\\.][\\d]*)?([Ee][+-]?[\\d]+)?");
    regex r(pattern);
    smatch results;
    string::const_iterator iter_begin = strin.cbegin();
    string::const_iterator iter_end = strin.cend();
    while (regex_search(iter_begin, iter_end,  results,  r)){
        if (i >= num) break;
        Vtemp.push_back(stof(results[0].str()));
        iter_begin = results[0].second;	
        i++;
    }
    if (i < num) {
        logger.error("Error code 1: getnum() cannot find enough numbers, expected " + to_string(num) + " numbers, but only found " + to_string(i) + " numbers.");
        logger.error("The involved line is " + strin);
        exit(1);
    }
    return Vtemp;
}

vector<double> get_vector(MatrixXd &matrix){
    vector<double> v;
    for(int i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
            v.push_back(matrix(i,j));
        }
    }
    return v;
}

vector<double> get_vector(VectorXd &matrix){
    vector<double> v;
    for(int i = 0; i < matrix.size(); i++){
        v.push_back(matrix(i));
    }
    return v;
}

MatrixXd cal_sn_info(MatrixXd &Min, vector<double> m_abc, vector<double> transMl, int Miller_n, int system_n){
    MatrixXd Min_s, Min_n, Mabc, Trans_Miller;
    Mabc = Eigen::Map<Matrix3d>(m_abc.data());
    Trans_Miller = Eigen::Map<MatrixXd>(transMl.data(),3,Miller_n);

    Min_n = Min(all,seq(0,Miller_n-1)) * Trans_Miller.transpose();
    Min_s = Min(all,seq(Miller_n,2*Miller_n-1)) * Trans_Miller.transpose(); 
    //calculate the coordinate in Cartesian system
    /* logger.debug("Initial Min_n = "); logger.debug(Min_n); */
    MatrixXd Mtemp = Min_n.cwiseAbs().array().max(0.001).matrix();
    Min_n = (Min_n.array().sign() * Mtemp.array()).matrix();
    Min_n = Min_n.cwiseInverse();
    MatrixX3d Min_nv1(Min_n.rows(),3), Min_nv2(Min_n.rows(),3);
    Min_nv1 << -Min_n.col(0), Min_n.col(1), VectorXd::Zero(Min_n.rows());
    Min_nv2 << VectorXd::Zero(Min_n.rows()), -Min_n.col(1), Min_n.col(2);
    Min_nv1 = Min_nv1 * Mabc.transpose();
    Min_nv2 = Min_nv2 * Mabc.transpose();
    /* logger.debug("Min_nv1 = "); logger.debug(Min_nv1); */
    /* logger.debug("Min_nv2 = "); logger.debug(Min_nv2); */
    for (int i = 0; i < Min_n.rows(); ++i) {
        Min_n.row(i) = Min_nv1.row(i).cross(Min_nv2.row(i));
    }
    Min_s = Min_s*Mabc.transpose();
    //normalization
    for(int i = 0; i < system_n; i++)  Min_n.row(i) = Min_n.row(i).normalized();
    /* logger.debug("Min_n = "); logger.debug(Min_n); */

    MatrixXd Min_ns(system_n,6);
    Min_ns.block(0,0,system_n,3) = Min_n;
    Min_ns.block(0,3,system_n,3) = Min_s;

    for(int i = 0; i < system_n; i++)
        for(int j = 0; j < 6; j++)
            if(abs(Min_ns(i,j)) <= 1e-3 ) Min_ns(i,j) = 0.0;

    return Min_ns;
}
