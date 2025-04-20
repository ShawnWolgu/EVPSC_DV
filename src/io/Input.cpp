#include "common/common.h"
#include "io/Input.h"

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
    ininp.open("p-SCOPE.in",ios::in);
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
    
    skipLine(loadinp); //Read Electric Field Vector control params
    getline(loadinp, tp);
    Vector3i IEfdot = getnum(tp,3).cast<int>();
    Proc.get_IEfdot(IEfdot);

    skipLine(loadinp); //Read Electric Field Vector
    getline(loadinp, tp);
    Vector3d Efdot = getnum(tp, 3);
    Proc.get_Efdot(Efdot);

    skipLine(loadinp); //Read Electric Disp Vector control params
    getline(loadinp, tp);
    Vector3i IEddot = getnum(tp,3).cast<int>();
    Proc.get_IEddot(IEddot);

    skipLine(loadinp); //Read Electric Disp Vector
    getline(loadinp, tp);
    Vector3d Eddot = getnum(tp, 3);
    Proc.get_Eddot(Eddot);

    if (!loadinp.eof()) {
        getline(loadinp, tp);
        if (tp.find("duty") != tp.npos) {
            // read current control params
            VectorXd coeff = getnum(tp, 11);
            tie(duty_ratio_J, Amplitude_J, Frequency, 
                ref_current_intensity_0, ref_current_intensity_1,
                ref_current_intensity_2, bvalue, shock_int,
                shock_fin, flag_emode, K_ew) = 
                tie(coeff[0], coeff[1], coeff[2], coeff[3], coeff[4],
                   coeff[5], coeff[6], coeff[7], coeff[8], coeff[9], coeff[10]);
            if (flag_emode != 0|| flag_emode != 1 || flag_emode != 2) {
                logger.warn("Error code 1: emode is not supported.");
            }
            logger.debug("duty_ratio_J = " + to_string(duty_ratio_J));
            logger.debug("Amplitude_J = " + to_string(Amplitude_J));
            logger.debug("Frequency = " + to_string(Frequency));
        }
        getline(loadinp, tp);
        if (tp.find("Current") != tp.npos) {
            // read current tensor
            Matrix3d J_tensor;
            for (int i = 0; i < 3; ++i) {
                getline(loadinp, tp);
                J_tensor.row(i) = getnum(tp, 3);
            }
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
    fstream sxinp; json sx_json;
    sxinp.open(fname,ios::in); //open .sx

    if (sxinp.is_open()) //checking whether the file is open
    {  
        logger.info("Loading sx file " + fname + " ...");
        string tp;     getline(sxinp, tp); //skip first line
        string crysym; getline(sxinp, crysym); //crystal symmetry string
        add_trans_miller(crysym.substr(0,5), sx_json);
        getline(sxinp, tp);  add_lattice_const(getnum(tp, 6), sx_json); //cystal constants
        int Millern = sx_json["Miller_n"];

        getline(sxinp, tp);  //skip a line;
        MatrixXd Cij6(6,6);  //Elastic constants;
        for (int i=0; i<6; i++) { getline(sxinp, tp); Cij6.row(i) = getnum(tp, 6);}
        add_elastic_constant(Cij6, sx_json);

        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);  VectorXd therm = getnum(tp, 6); //Thermal coefficients
        add_thermal_coefficient(therm, sx_json);
        getline(sxinp ,tp); // this line is for thermal coeffs check or read plasitcity modes
        if (tp.find("rho_material") != tp.npos){
            getline(sxinp, tp); VectorXd thermal_coeff = getnum(tp, 7);
            rho_material = thermal_coeff(0);
            Cp_material = thermal_coeff(1);
            sigma_e_mat = thermal_coeff(2);
            h_ext = thermal_coeff(3);
            Surface = thermal_coeff(4);
            V_sample = thermal_coeff(5);
            sigma_k = thermal_coeff(6);
            getline(sxinp, tp);  //skip a line;        //Start reading slip and twinning modes
        }
        else{
            rho_material = 0;
            Cp_material = 0;
            sigma_e_mat = 0;
            h_ext = 0;
            Surface = 0;
            V_sample = 0;
            sigma_k = 0;
        }
        logger.debug("rho_material = " + to_string(rho_material));
        logger.debug("Cp_material = " + to_string(Cp_material));
        logger.debug("sigma_e_mat = " + to_string(sigma_e_mat));
        logger.debug("h_ext = " + to_string(h_ext));
        logger.debug("Surface = " + to_string(Surface));
        logger.debug("V_sample = " + to_string(V_sample));
        logger.debug("sigma_k = " + to_string(sigma_k));
        //关于传热的直接在这边赋值
        getline(sxinp, tp);  int nmodesx = int(getnum(tp, 1)(0)); //total mode number in file
        getline(sxinp, tp);  int nmodes = int(getnum(tp, 1)(0));  //considered in current run
        getline(sxinp, tp);  VectorXd mode_i = getnum(tp, nmodes);  //the index of modes(mode_i)

        //Start reading slip and twinning modes
        int modes_num = 0; vector<int> mode_count; vector<json> sx_modes;
        bool f = 1;  
        for(int imode = 0; imode < nmodesx; ++imode){
            getline(sxinp, tp);  //skip a line;
            getline(sxinp, tp);  VectorXd mode_info = getnum(tp, 4);
            /* mode_info 0: the serial number 
             * 1: number of deformation systems
             * 2: flag of slip (0 for twin; 1 for slip) 
             * 3: flag of twin (1 for twin; 0 for slip)
            */ 
            MatrixXd nor_dir(int(mode_info(1)),2*Millern); //normal and direction of slip plane
            /* if(int(mode_info(3))) getline(sxinp, tp); //special for twin */
            for (int i = 0; i < int(mode_info(1)); i++) {
                getline(sxinp, tp); nor_dir(i,all) = getnum(tp, 2*Millern);
            }
            f = (mode_i.array() == imode+1).any();
            if(f){
                json this_mode;
                MatrixXd sn_matrix = cal_sn_info(nor_dir, sx_json["Mabc"], sx_json["Trans_Miller"], Millern, int(mode_info(1)));
                this_mode["sn_info"] = get_vector(sn_matrix);
                this_mode["mode_n"] = int(mode_info(1));
                if (mode_info(2) == 1) this_mode["type"] = 0; //slip
                else if (mode_info(3) == 1) this_mode["type"] = 1; //twin
                else this_mode["type"] = 2; //other
                modes_num += int(mode_info(1));
                mode_count.push_back(int(mode_info(1)));
                sx_modes.push_back(this_mode);
            }                
        }
        sx_json["family_num"] = nmodes;
        sx_json["modes_count_by_family"] = mode_count;
        sx_json["modes_num"] = modes_num;

        getline(sxinp, tp);  //skip a line;
        getline(sxinp, tp);  int iharden = int(getnum(tp, 1)(0)); //hardening law(Voce=0, DV=1)
        getline(sxinp, tp);  bool irate = bool(getnum(tp, 1)(0)); //"rate sensitive" flag(1: Y; 0: N)
        getline(sxinp, tp);  sx_json["GZ"] = getnum(tp, 1)(0); //grain size: um
        int harden_size;
        if(iharden == 0) harden_size = 4; else harden_size = 18;

        //Read hardening parameters of modes
        double nrsx; vector<double> CRSS_p, hst;
        for(int imode = 0; imode < nmodes; ++imode)
        {
            getline(sxinp, tp);  //skip a line;
            getline(sxinp, tp);  nrsx = getnum(tp, 1)(0); //rate sensitivity
            getline(sxinp, tp);  //CRSS parameters
            if (sx_modes[imode]["type"] == 0) CRSS_p = getnum_vec(tp, harden_size);
            else CRSS_p = getnum_vec(tp, 10);
            getline(sxinp, tp);  //latent hardening parameters
            if (sx_modes[imode]["type"] == 0) hst = getnum_vec(tp, 6); //6 types of hardening
            else hst = getnum_vec(tp, 2);
            sx_modes[imode]["nrsx"] = nrsx;
            sx_modes[imode]["CRSS_p"] = CRSS_p;
            sx_modes[imode]["hst"] = hst;
        }
        sx_json["modes"] = sx_modes;
        json sx_per_mode = sx_info_postprocess(sx_json);
        /* conjugate_mode_config(sx_per_mode); */
        sx_json["sx_per_mode"] = sx_per_mode;
        sxinp.close(); //close the file object.
        pcrys.ini_from_json(sx_json);
        /* pcrys.g[0].gmode[0]->print(); */
        return 0;
    }
    else
    {
    logger.error("Error code 0: .sx file cannot be opened");
    return 1;
}
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

void add_trans_miller(string crysym, json &sx_json){
    //calculate conversion matrix of Miller indices according to the crysym
    //transform str to lower case
    for (int i = 0; i < crysym.size(); i++) crysym[i] = tolower(crysym[i]);
    vector<double> Mtemp; int Miller_n; 
    if(!crysym.compare("hexag"))
        {
            Miller_n = 4;
            Mtemp = {1, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 1};
        }
    else if(!crysym.compare("cubic")) 
        {
        Miller_n = 3;
        Mtemp = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    }
    else
    {
        logger.error("Error code 1: crystal symmetry is not supported.");
        exit(1);
    }
    sx_json["crysym"] = crysym;
    sx_json["Miller_n"] = Miller_n;
    sx_json["Trans_Miller"] = Mtemp;
}

void add_lattice_const(VectorXd ccon, json &sx_json){
    //add lattice constants to sx_json
    vector<double> Cdim, Cang, Mtemp; Matrix3d Mabc; 
    for(int i = 0; i < ccon.size(); i++){
        if (i < 3) Cdim.push_back(ccon(i));
        else Cang.push_back(ccon(i)/180*M_PI);
    }
    sx_json["Cdim"] = Cdim;
    sx_json["Cang"] = Cang;
    //calculate Mabc
    Mabc(0,0)=sin(Cang[1]);
    Mabc(1,0)=0.;
    Mabc(2,0)=cos(Cang[1]);
    Mabc(0,1)=(cos(Cang[2])-cos(Cang[0])*cos(Cang[1]))/sin(Cang[1]);
    Mabc(2,1)=cos(Cang[0]);
    Mabc(1,1)=sqrt(1.0-pow(Mabc(0,1),2)-pow(Mabc(2,1),2));
    Mabc(0,2)=0.;
    Mabc(1,2)=0.;
    Mabc(2,2)=1.;
    for(int i = 0; i < 3; i++){ for(int j = 0; j < 3; j++) Mtemp.push_back(Cdim[j] * Mabc(i,j));};
    sx_json["Mabc"] = Mtemp;
}

void add_elastic_constant(MatrixXd Cij6, json &sx_json){
    vector<double> ela_consts;
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            ela_consts.push_back(Cij6(i,j));
        }
    }
    sx_json["Cij6"] = ela_consts;
}

void add_thermal_coefficient(VectorXd ther, json &sx_json){
    vector<double> ther_consts;
    for(int i = 0; i < ther.size(); i++) ther_consts.push_back(ther(i));
    sx_json["therm"] = ther_consts;
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

json sx_info_postprocess(json &sx_json){
    vector<json> info_per_mode; int mode_id = 0;
    for (auto &json_mode : sx_json["modes"]){
        int system_n = json_mode["mode_n"];
        for (int crt_mode=0; crt_mode < system_n; crt_mode++){
            json this_mode;
            this_mode["id"] = mode_id++;
            this_mode["type"] = json_mode["type"]; //other
            vector<double> sn;
            for (int i = 6*crt_mode; i < 6*crt_mode+6; i++) sn.push_back(json_mode["sn_info"][i]);
            this_mode["nrsx"] = json_mode["nrsx"];
            this_mode["CRSS_p"] = json_mode["CRSS_p"];
            this_mode["hst"] = json_mode["hst"];
            this_mode["sn"] = sn;
            this_mode["G"] = cal_shear_modulus(sx_json["Cij6"], sn);
            info_per_mode.push_back(this_mode);
        }
    }
    return info_per_mode;
}

double cal_shear_modulus(vector<double> Cij6, vector<double> sn){
    Matrix3d slip_rotation; Matrix6d elastic_modulus;
    Vector3d plane_norm, burgers_vec, trav_direc;
    double shear_modulus;
    plane_norm << sn[0], sn[1], sn[2]; burgers_vec << sn[3], sn[4], sn[5];
    elastic_modulus = Eigen::Map<Matrix6d>(Cij6.data());
    trav_direc = burgers_vec.cross(plane_norm);
    slip_rotation << (burgers_vec/burgers_vec.norm()), plane_norm, trav_direc / trav_direc.norm();
    shear_modulus = rotate_6d_stiff_modu(elastic_modulus, slip_rotation.transpose())(3,3);
    return shear_modulus;
}
