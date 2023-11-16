#include "Grains.h"
#include "Eigen/src/Core/Matrix.h"
#include "Toolbox.h"
#include "func.h"
#include "global.h"
#include <string>

grain::grain()
{
    //initial the grain stress&strain
    eps_g = Matrix3d::Zero();   
    sig_g = Matrix3d::Zero();
    Dije_g = Matrix3d::Zero();
    Dijp_g = Matrix3d::Zero();

    //initial the shape of ellipsoid
    ell_axis_g = Vector3d::Ones();
    ell_axisb_g = Matrix3d::Identity();
    Fij_g = Matrix3d::Identity();
    d0_g = Vector5d::Zero();

    //initial the VP consistent
    Mpij6_g = 1e-10 * Matrix5d::Identity();

}

Vector3d grain::get_ell_axis_g(){return ell_axis_g;}

Matrix3d grain::get_ell_axisb_g(){return ell_axisb_g;}

void grain::ini_euler_g(Vector4d vin)
{
    Vector3d euler = vin(seq(0,2));
    Euler_M = Euler_trans(euler);
    weight = vin(3);
}

Vector3d grain::get_euler_g(){return Euler_trans(Euler_M);}

Matrix3d grain::get_Euler_M_g(){return Euler_M;}

double grain::get_weight_g(){return weight;}

void grain::set_weight_g(double w){weight = w;}

int grain::ini_gmode_g(int n)
{
    modes_num = n;
    gmode = new Slip[n];
    gamma_delta_gmode = new double[n];
    return 1;
}

int grain::ini_gmode_g(json &j)
{
    modes_num = j["modes_num"];
    gmode = new Slip[modes_num];
    int i = 0;
    for (json &j_mode: j["sx_per_mode"]){
        gmode[i] = Slip(j_mode); ++i;
    }
    gamma_delta_gmode = new double[modes_num];
    return 1;
}

int grain::check_gmode_g(){return modes_num;}

int grain::ini_sn_g(MatrixXd Min, int flag, int system_n, int modei, Matrix6d elastic_modulus)
{
    cout << "Min" << Min << "\n";
    for(int i = modei; i < modei+system_n; i++)
    {
        VectorXd Min_temp = Min(i-modei,all);
        cout << "line number: " << i-modei << "\n";
        cout << "Min_temp: " << Min_temp.transpose() << "\n";
        gmode[i].ini_sn_mode(Min_temp, flag, i);
        gmode[i].cal_shear_modulus(elastic_modulus);
    }
    return 0;
}

int grain::check_sn_g()
{
    for(int i = 0; i < modes_num; i++)
    {
        cout << "Mode " << i << ":\n";
        gmode[i].check_sn_mode();
    }
    return 0;
}

int grain::ini_hardening_g(double nrsx_in, VectorXd CRSS_p_in, VectorXd hst_in, int modei, int modes_num)
{
    for(int i = modei; i < modei+modes_num; i++)
    {
        gmode[i].ini_hardening_mode(nrsx_in, CRSS_p_in, hst_in);
    }
    return 0;
}

int grain::check_hardening_g()
{
    for(int i = 0; i < modes_num; i++)
    {
        logger.debug("Mode "+std::to_string(i)+":");
        gmode[i].check_hardening_mode();
    }
    return 0;    
}

void grain::Eshelby_E(double ESIM[3][3][3][3], double ESCR[3][3][3][3], \
                      Vector3d axis_t, Matrix6d C66,Integralpoint6 aa6, Integralpoint6 aaww6, Integralpoint3 alpha)
{
    int npoints = Intn * Intn; 

    double aixabc = axis_t(0)*axis_t(1)*axis_t(2);
    //double ratio1 = axis_t(1)/axis_t(2);
    //double ratio2 = axis_t(0)/axis_t(2);


    Matrix6d T66; //Tijkl
    T66 = Matrix6d::Zero();

    Vector6d aa2, aaww2, a1;
    Vector6d a1_inv;
    Matrix3d a1m;

    double Ro3;
    double abcoro3;
    for(int ny = 0; ny < npoints; ny++)
    {
        aa2 = aa6.col(ny); //Take out the Alpha*Alpha
        a1 = Mult_voigt(C66,aa2);
        //C66(i,j,k,l)*aa2(j,l)
        //If solving an elastic inclusion invert the system
        //A(3,3) x X(3,3) = C(3,3)
        a1m = voigt(a1).inverse();
        a1_inv = voigt(a1m);

        Ro3 =0;
        for(int i = 0; i < 3; i++)
            Ro3 += pow(alpha(i,ny)*axis_t(i), 2);
        Ro3 = pow(Ro3, 1.5);
        abcoro3 = aixabc/Ro3;

        //T(M,N,I,J)=A(M)*A(J)*G(N,I)
        for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
                T66(i,j) += aaww6(i,ny)*a1_inv(j)*abcoro3;
    }
    double Tijkl[3][3][3][3];
    voigt(T66, Tijkl);

    //   COMPUTE SYMMETRIC (DISTORTION) ESHELBY TENSOR FROM EQ.B9.
    //       ESIM(N,M,K,L)=0.5*(T(M,J,N,I)+T(N,J,M,I))*C4(I,J,K,L)
    //   COMPUTE ANTI-SYMMETRIC (ROTATION) ESHELBY TENSOR FROM EQ.B9.
    //       ESCR(N,M,K,L)=0.5*(T(M,J,N,I)-T(N,J,M,I))*C4(I,J,K,L)
    double C4[3][3][3][3];
    voigt(C66, C4);

    double dumsim, dumscr;
    for(int l = 0; l < 3; l++)
        for(int k = 0; k < 3; k++)
            for(int m = 0; m < 3; m++)
                for(int n = 0; n < 3; n++)
                {
                    dumsim = 0;
                    dumscr = 0;
                    for(int j = 0; j < 3; j++)
                        for(int i = 0; i < 3; i++)
                        {
                            dumsim=dumsim+(Tijkl[m][j][n][i]+Tijkl[n][j][m][i])*C4[i][j][k][l];
                            dumscr=dumscr+(Tijkl[m][j][n][i]-Tijkl[n][j][m][i])*C4[i][j][k][l];
                        }
                    ESIM[n][m][k][l]=0.5*dumsim;
                    ESCR[n][m][k][l]=0.5*dumscr;
                }
}

void grain::Eshelby_P(double ESIM[3][3][3][3],double ESCR[3][3][3][3],\
                      Vector3d axis_t, Matrix6d C66,Integralpoint6 aa6, Integralpoint6 aaww6, Integralpoint3 alpha, Integralpoint3 aww, Integralpoint1 ww)
{
    int npoints = Intn * Intn; 

    double aixabc = axis_t(0)*axis_t(1)*axis_t(2);
    //double ratio1 = axis_t(1)/axis_t(2);
    //double ratio2 = axis_t(0)/axis_t(2);

    Matrix6d T66; //Tijkl
    T66 = Matrix6d::Zero();

    Matrix3d P;
    P = Matrix3d::Zero();
    double PDIL = 0;

    Vector6d aa2, aaww2, a1;
    Vector6d a1_inv;
    Vector10d a1_10;
    Vector10d a1_10_inv;
    Matrix4d a1m;

    double Ro3;
    double abcoro3;
    for(int ny = 0; ny < npoints; ny++)
    {
        aa2 = aa6.col(ny); //Take out the Alpha*Alpha
        a1 = Mult_voigt(C66,aa2); //C66(i,j,k,l)*aa2(j,l)
        //refer to Equ[5-14b] in manual 7d)
        //solving a visco-plastic inclusion defines componets A1(7) TO A1(10).
        //solve the system given by A(4,4) X X(4,4) = C(4,4)
        a1_10(seq(0,5)) = a1;
        a1_10(6) = alpha(0,ny);
        a1_10(7) = alpha(1,ny);
        a1_10(8) = alpha(2,ny);
        a1_10(9) = 0;
        a1m = voigt(a1_10).inverse();
        a1_10_inv = voigt(a1m);

        Ro3 =0;
        for(int i = 0; i < 3; i++)
            Ro3 += pow(alpha(i,ny)*axis_t(i), 2);
        Ro3 = pow(Ro3, 1.5);
        abcoro3 = aixabc/Ro3;

        //T(M,N,I,J)=A(M)*A(J)*G(N,I)
        for(int i = 0; i < 6; i++)
            for(int j = 0; j < 6; j++)
                T66(i,j) += aaww6(i,ny)*a1_10_inv(j)*abcoro3;

        for(int j = 0; j < 3; j++)
            for(int i = 0; i < 3; i++)
                P(i,j) += aww(j,ny)*a1_10_inv(i+6)*abcoro3;
        PDIL += ww(ny)*a1_10_inv(9)*abcoro3;

    }
    double Tijkl[3][3][3][3];
    voigt(T66, Tijkl);

    //   COMPUTE SYMMETRIC (DISTORTION) ESHELBY TENSOR FROM EQ.B9.
    //       ESIM(N,M,K,L)=0.5*(T(M,J,N,I)+T(N,J,M,I))*C4(I,J,K,L)
    //   COMPUTE ANTI-SYMMETRIC (ROTATION) ESHELBY TENSOR FROM EQ.B9.
    //       ESCR(N,M,K,L)=0.5*(T(M,J,N,I)-T(N,J,M,I))*C4(I,J,K,L)
    double C4[3][3][3][3];
    voigt(C66, C4);

    double dumsim, dumscr;
    for(int l = 0; l < 3; l++)
        for(int k = 0; k < 3; k++)
            for(int m = 0; m < 3; m++)
                for(int n = 0; n < 3; n++)
                {
                    dumsim = 0;
                    dumscr = 0;
                    for(int j = 0; j < 3; j++)
                        for(int i = 0; i < 3; i++)
                        {
                            dumsim=dumsim+(Tijkl[m][j][n][i]+Tijkl[n][j][m][i])*C4[i][j][k][l];
                            dumscr=dumscr+(Tijkl[m][j][n][i]-Tijkl[n][j][m][i])*C4[i][j][k][l];
                        }
                    ESIM[n][m][k][l]=0.5*dumsim;
                    ESCR[n][m][k][l]=0.5*dumscr;
                }
}

Matrix3d grain::get_stress_g(){return sig_g;}
Matrix3d grain::get_strain_g(){return eps_g;}

Matrix3d grain::get_Dije_g(){return Dije_g;}
Matrix3d grain::get_Dijp_g(){return Dijp_g;}
Matrix3d grain::get_Udot_g(){return Udot_g;}

void grain::save_status_g(){
    sig_g_old = sig_g;
    Mij6_J_g_old = Mij6_J_g;  Metilde_g_old = Metilde_g;
    Mpij6_g_old = Mpij6_g;    Mptilde_g_old = Mptilde_g;
    Cij6_SA_g_old = Cij6_SA_g;
    Dij_g_old = Dij_g;        
    Dije_g_old = Dije_g;      Dijp_g_old = Dijp_g;
    save_RSinv_g();
}

void grain::save_RSinv_g(){
    for(int i = 0; i < 3; i++)  for(int j = 0; j < 3; j++)
        for(int m = 0; m < 3; m++)  for(int n = 0; n < 3; n++){
            RSinv_C_old[i][j][m][n] = RSinv_C[i][j][m][n];
            RSinv_VP_old[i][j][m][n] = RSinv_VP[i][j][m][n];
        }
}

void grain::restore_status_g(){
    sig_g = sig_g_old;
    Mij6_J_g = Mij6_J_g_old;  Metilde_g = Metilde_g_old;
    Mpij6_g = Mpij6_g_old;    Mptilde_g = Mptilde_g_old;
    Cij6_SA_g = Cij6_SA_g_old;
    Dij_g = Dij_g_old;
    Dije_g = Dije_g_old;      Dijp_g = Dijp_g_old;
    Update_RSinv_C_g(RSinv_C_old);
    Update_RSinv_VP_g(RSinv_VP_old);
}
//elastic consistent
void grain::Update_Cij6_SA_g(Matrix6d Min){Cij6_SA_g = Min;}
void grain::Update_Mij6_J_g(Matrix6d Min){Mij6_J_g = Min;}
void grain::Update_Metilde_g(Matrix6d Min){Metilde_g = Min;}
void grain::Update_RSinv_C_g(double A[3][3][3][3])
{
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int m = 0; m < 3; m++)
                for(int n = 0; n < 3; n++)
                    RSinv_C[i][j][m][n] = A[i][j][m][n];
}
Matrix6d grain::get_Mij6_J_g(){return Mij6_J_g;}


//visco-plastic consistent
void grain::Update_Mptilde_g(Matrix5d Min){Mptilde_g = Min;}
Matrix5d grain::get_Mpij6_g(){return Mpij6_g;}
Vector5d grain::get_d0_g(){return d0_g;}
void grain::Update_RSinv_VP_g(double A[3][3][3][3])
{
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int m = 0; m < 3; m++)
                for(int n = 0; n < 3; n++)
                    RSinv_VP[i][j][m][n] = A[i][j][m][n];
}

void grain::Update_Mpij6_g()
{
    Mpij6_g = 1e-10 * Matrix5d::Identity() + cal_Fgrad(sig_g) ;
    d0_g = Chg_basis5(Dijp_g) - Mpij6_g * Chg_basis5(sig_g);
}

void grain::Update_shape_g()
{
    //-1 F.transpose()*F
    //BX = Matrix3d::Zero();
    //for(int i = 0; i < 3; i++)
    //    for(int j = 0; j < 3; j++)
    //        BX(i,j) = Fij_g.row(i)*Fij_g.col(j);
    Matrix3d BX;
    Matrix3d FT = Fij_g.transpose();
    BX =  Fij_g*FT;
    //-1 F.transpose()*F

    //-2 solve the eigen vector of BX
    //and sort the value from largest to smallest
    EigenSolver<MatrixXd> es(BX);
    Matrix3d BX_vectors;
    Vector3d BX_value;
    BX_value = (es.eigenvalues()).real();
    BX_vectors = (es.eigenvectors()).real();
    Eigsrt(BX_vectors, BX_value);
    //-2 solve the eigen vector of BX
    //and sort the value from largest to smallest

    Matrix3d B = BX_vectors;
    Vector3d W = BX_value;
    //-3 
    //redefine Axis(1) (the second) to be the largest  
    //to improve the accuracy in calculation of Eshelby tensor
    //IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
    //HANDED BY EXCHANGING 1 AND 2.
    double Sign = -1;
    double temp;
    if(B.determinant() <= 0) Sign = 1;
    for(int i = 0; i < 3; i++)
    {
        temp = B(i,0);
        B(i,0) = B(i,1);
        B(i,1) = temp * Sign;
    }
    temp = W(0); W(0) = W(1); W(1) = temp;
    //-3 

    //-4 update the stretching of ellipsoid
    double Ratmax=sqrt(W(1)/W(2));
    double Ratmin=sqrt(W(0)/W(2));

    ell_axisb_g = B;
    if(!Iflat_g) //if Iflat_g = 0
        for(int i = 0; i < 3; i++)
            ell_axis_g(i) = sqrt(W(i));
    //if Iflat_g = 1, the axis of ellipsoid keeps unchange
    //-4 update the stretching of ellipsoid

    //-5 update the ellipsoid orientation
    Matrix3d BT = B.transpose();
    ellip_ang_g = Euler_trans(B);
    //-5

    //-6 Update the Iflat_g according to the Max axes ratio of ellipsoid
    if((!Iflat_g)&&(Ratmax >= ell_crit_shape_g))
        {
            Iflat_g = 1;
        }
        //-6

        //-7 Iflat_g = 1; recaculates the Fij_g in grain
    else if((Iflat_g)&&(Ratmax >= ell_crit_shape_g))
    {
    W(1) = W(1)/4;
    if(Ratmin >= 0.5 * ell_crit_shape_g)
        W(0) = W(0)/4;
    for(int i = 0; i < 3; i++)
        ell_axis_g(i) = sqrt(W(i));

    Matrix3d Fijx;
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            Fijx(i,j) = Iij(i,j) * ell_axis_g(i);

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            Fij_g(i,j) = 0;
            for(int k = 0; k < 3; k++)
                for(int l = 0; l < 3; l++)
                    Fij_g(i,j) = Fij_g(i,j) + B(i,k)*B(j,l)*Fijx(k,l);
        }
}
}

void grain::Update_Fij_g(double Tincr)
{
    Matrix3d Fnew;
    Fnew = Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
                Fnew(i,j) += (Tincr*Udot_g(i,k)+Iij(i,k))*Fij_g(k,j);
    Fij_g = Fnew;
}

Matrix3d grain::cal_Dijp(Matrix3d Min)
{
    //transform into the deviatoric tensor
    Matrix3d X = devia(Min);
    Matrix3d E = Euler_M;
    Matrix3d ET = Euler_M.transpose();
    X = E * X * ET;

    Matrix3d Dijp = Matrix3d::Zero();
    for(int i = 0; i < modes_num; i++)
    {
        Dijp += gmode[i].cal_dijpmode(X);
    }
    return ET*Dijp*E;
}

Matrix3d grain::cal_rotslip()
{
    Matrix3d E = Euler_M;
    Matrix3d ET = E.transpose();

    Matrix3d Wij = Matrix3d::Zero();
    for(int i = 0; i < modes_num; i++)
    {
        Wij += gmode[i].cal_rotslip_m();
    }
    return ET*Wij*E;
}

Matrix5d grain::cal_Fgrad(Matrix3d Min)
{
    Matrix3d X = devia(Min);
    Matrix3d E = Euler_M;
    Matrix3d ET = E.transpose();
    X = E * X * ET;

    Matrix6d Fgrad = Matrix6d::Zero();
    for(int i = 0; i < modes_num; i++)
        Fgrad += gmode[i].get_Fgradm(X);
    //return Chg_basis5(rotate_C66(Fgrad, E));
    return Chg_basis5(rotate_C66(Fgrad, ET));
}

double grain::cal_RSSxmax(Matrix3d Min)
{
    Matrix3d X = devia(Min);
    Matrix3d E = Euler_M;
    Matrix3d ET = Euler_M.transpose();
    X = E * X * ET;
    double RSSxmax = 0;
    for(int i = 0; i < modes_num; i++)
    {
        double temp = gmode[i].cal_relative_rss(X);
        if(RSSxmax < temp) RSSxmax = temp;
    }
    return RSSxmax;
}

double grain::cal_RSSxlim(Matrix3d D)
{   
    double nrsxmin = 0;
    for(int i = 0; i < modes_num; i++)
        if(nrsxmin > gmode[i].get_nrsx())
            nrsxmin = gmode[i].get_nrsx();

    double lim = 2* pow(D.norm()/gmode[0].get_gamma0(),1/nrsxmin);
    if(lim < 2) lim = 2;
    return  lim;           
}

void grain::grain_stress(double Tincr, Matrix3d Wij_m, Matrix3d Dij_m,\
                         Matrix3d Dije_AV, Matrix3d Dijp_AV, Matrix3d Sig_m, Matrix3d Sig_m_old)
{
    Matrix3d X = sig_g; 
    //according to Equ[5-30] in manual, calculate the wij~ and the Jaumann rate
    Matrix3d wg = Wij_m + mult_4th(RSinv_C,Dije_g-Dije_AV) + mult_4th(RSinv_VP,Dijp_g-Dijp_AV);
    Matrix3d Xjau = wg * X - X * wg;
    Matrix3d Sjau = Wij_m * Sig_m - Sig_m * Wij_m; 

    //solve the interaction equation of grain
    Vector6d DB = Chg_basis6(Dij_m);
    Matrix3d Mt = Sig_m - Sig_m_old + sig_g_old;
    DB += Metilde_g * Chg_basis6(Mt)/Tincr \
        + Mij6_J_g * Chg_basis6(sig_g_old)/Tincr; 
    Mt = Xjau-Sjau;
    DB += Mij6_J_g *Chg_basis6(Xjau) + Metilde_g * Chg_basis6(Mt);
    DB = Bbasisadd(DB, Mptilde_g * Chg_basis5(Sig_m));

    //Newton-Rapthon iteration
    Vector6d Xv = Chg_basis6(X), Xold = Xv, F;
    /* logger.debug("initial Xv:"); logger.debug(Xv); */
    Vector5d dijpgv;
    Matrix6d Fgrad, Mtemp = (Metilde_g + Mij6_J_g)/Tincr;
    double Errm = 1e-3, F_err = 1e-4, err_iter = Errm*10, coeff = 1;
    int NR_max_iter = 7, DH_max_iter = 25;
    for(int it = 0; it < NR_max_iter; it ++ )
    {
        F = -DB +  Mtemp * Xv;
        F =  Bbasisadd(F, Mptilde_g * Chg_basis5(Chg_basis(Xv)));
        dijpgv = Chg_basis5(cal_Dijp(Chg_basis(Xv)));
        F = Bbasisadd(F,dijpgv);
        Fgrad = -Bbasisadd(Mtemp,Mptilde_g+cal_Fgrad(Chg_basis(Xv)));
        Xold = Xv;
        Xv = Xold + Fgrad.inverse() * F;
        err_iter = Errorcal(Xv, Xold);
        /* if (isnan(err_iter)){ */
        /*     logger.debug("F_grad: "); logger.debug(Fgrad);  */
        /*     logger.debug("Mtemp: "); logger.debug(Mtemp); */
        /*     logger.debug("Mptilde_g: "); logger.debug(Mptilde_g); */
        /* } */
        /* if(grain_i == 0) logger.debug("Grain stress iteration: " + to_string(it) + ", F_norm: " + to_string(F.norm()) + ", err_iter: " + to_string(err_iter)); */
        if(err_iter < Errm) break;
    }
    if (isnan(err_iter)){
        logger.warn("Grain stress iteration failed (nan)!");
        Xv = Chg_basis6(X) + Chg_basis6(Sig_m) - Chg_basis6(Sig_m_old);
    }
    /* logger.debug("Xv after NR iteration: "); logger.debug(Xv); */
    /* logger.debug("F norm after NR iteration: " + to_string(F.norm())); */
    if(F.norm() >= F_err){
        /* logger.warn("Entering downhill ..."); */
        double F_norm = F.norm();
        for (int it = 0; F_norm >= F_err; it++) {
            F = -DB + Mtemp * Xv;
            F = Bbasisadd(F, Mptilde_g * Chg_basis5(Chg_basis(Xv)));
            dijpgv = Chg_basis5(cal_Dijp(Chg_basis(Xv)));
            F = Bbasisadd(F, dijpgv);
            if (F.norm() < F_norm) {
                F_norm = F.norm();
                coeff = min(1., coeff * 2);
                Xold = Xv;
                Fgrad = -Bbasisadd(Mtemp, Mptilde_g + cal_Fgrad(Chg_basis(Xv)));
                Xv = Xold + coeff * Fgrad.inverse() * F;
            }
            else{
                coeff *= 0.5;
                Xv = Xold + (Xv - Xold) * 0.5;
            }
            /* logger.debug("Grain stress (downhill) iteration: " + to_string(it)); */
            /* logger.debug(Xv); */
            /* if(grain_i == 0) logger.debug("Grain stress iteration: " + to_string(it) + ", F_norm: " + to_string(F_norm) + ", coeff: " + to_string(coeff)); */
            if ((it == DH_max_iter) || (coeff < 1e-4)) {
                logger.warn("Grain stress iteration failed !");
                Xv = Chg_basis6(X) + Chg_basis6(Sig_m) - Chg_basis6(Sig_m_old);
                break;
            }
        }
    }
    /**/
    sig_g = Chg_basis(Xv);
    /* logger.debug("Grain number = " + to_string(grain_i) + ", sig_g = "); logger.debug(sig_g); */
    Vector6d dijegv =  Mij6_J_g *((Xv - Chg_basis6(sig_g_old))/Tincr - Chg_basis6(Xjau)); 
    dijpgv = Chg_basis5(cal_Dijp(sig_g));
    Dije_g = Chg_basis(dijegv); Dijp_g = Chg_basis(dijpgv);
    Dij_g = Dije_g + Dijp_g;

    /* logger.debug("Grain number = %d" + to_string(grain_i)); */
    /* logger.debug(Mpij6_g); */
    Update_Mpij6_g();
}

void grain::Update_shear_strain(double Tincr)
{
    double temp = 0;
    gamma_delta = 0;
    for(int i = 0; i < modes_num; i++)
    {         
        temp = Tincr * gmode[i].update_shear_strain_m();
        gamma_delta_gmode[i] = temp;
        gamma_delta += temp;
    }    
}

void grain::update_strain(double Tincr)
{
    eps_g += Dij_g * Tincr;
}

Matrix3d grain::get_Wij_g()
{
    return Wij_g;
}

void grain::update_orientation(double Tincr, Matrix3d Wij_m, Matrix3d Dije_AV, Matrix3d Dijp_AV)
{
    Wij_g = Wij_m + mult_4th(RSinv_C,Dije_g-Dije_AV) + mult_4th(RSinv_VP,Dijp_g-Dijp_AV);

    Udot_g = Wij_g + Dij_g; //update the velocity gradient in grains

    Matrix3d Rot = (Wij_g - cal_rotslip())*Tincr;
    Matrix3d Euler_M_new = Euler_M*Rodrigues(Rot).transpose();
    if (Errorcal(Euler_M, Euler_M_new) < 0.33) Euler_M = Euler_M_new;
    /* Euler_M = Euler_M_new; */
}

void grain::update_modes(double Tincr)
{
    for(int i = 0; i < modes_num; i++)	gmode[i].update_ssd(Dij_g, Tincr);
    /* for(int i = 0; i < modes_num; i++)	gmode[i].update_lhparams(Dij_g); */
    for(int i = 0; i < modes_num; i++)	gmode[i].update_status(*this, Tincr);
    gamma_total += gamma_delta;
}

void grain::set_lat_hard_mat(){
    lat_hard_mat.resize(modes_num,modes_num);
    for (int i = 0; i < modes_num; i++){
        for (int j = 0; j < modes_num; j++){
            if (gmode[i].num == gmode[j].num) lat_hard_mat(gmode[i].num,gmode[j].num) = 1;
            else{
                int mode = get_interaction_mode(gmode[i].burgers_vec, gmode[i].plane_norm, gmode[j].burgers_vec, gmode[j].plane_norm);
                lat_hard_mat(gmode[i].num,gmode[j].num) = gmode[i].latent_params[mode];
            }
        }
    }
}

int grain::get_interaction_mode(Vector3d burgers_i, Vector3d plane_i, Vector3d burgers_j, Vector3d plane_j){
    /*
     * Return the dislocation interaction mode code between two slip system.
     * 0: No Junction, 1: Hirth Lock, 2: Coplanar Junction, 3: Glissile Junction, 4: Sessile Junction
     */
    double perp = 0.02, prll = 0.98;
    double cos_b_angle = cal_cosine(burgers_i, burgers_j);
    if(abs(cos_b_angle) < perp) return 1;
    else {
        if(abs(cos_b_angle) > prll) return 0;
        else{
            if (abs(cal_cosine(plane_i, plane_j)) > prll) return 2;
            else{
                bool if_glide_i = (abs(cal_cosine(plane_i, burgers_i+burgers_j)) < perp);
                bool if_glide_j = (abs(cal_cosine(plane_j, burgers_i+burgers_j)) < perp);
                if (if_glide_i || if_glide_j) return 3;
                else return 4;
            }
        }
    }
}
