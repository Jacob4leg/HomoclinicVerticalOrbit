#include "prove_fixed_point.h"
#include<fstream>

using namespace std;
using namespace capd;

CR3BP<long double> vf;
ICR3BP Ivf;

long double z = 0.05316795030019990465013731L;
LDVector v0 = vf.findVerticalLyapunovOrbit(LDVector({0.94690970780356629,0,z,0,-0.011166782657715382,0}));
long double E0 = vf.E(v0)[0];

LDMatrix T_E = energy_change_of_basis(v0,vf.E);

LDMatrix get_change_of_basis() {
    LDMatrix D = LDMatrix::Identity(6), T(6,6);
    LDVector u = v0;
    for(int i = 0; i < 2; i++) {
        u = vf.pm_y(u,T);
        D = vf.pm_y.computeDP(u,T) * D;
    }

    D = matrixAlgorithms::gaussInverseMatrix(T_E) * D * T_E;    
    D = matrix_add_cord(matrix_erase_cord(D,5),5);
    D[5][5] = 1.;
    
    LDMatrix D_new = matrix_erase_cord(D,1);
    

    LDMatrix L = us_change_of_basis(D_new);
    
    // cout << matrixAlgorithms::gaussInverseMatrix(L) * D_new * L << endl;
    return matrix_add_cord(L,1);
}


bool prove_fixed_point(long double eps_x, long double eps_y, long double E_eps) {
    
    LDMatrix T = get_change_of_basis();
    IMatrix T_total = IMatrix(T_E * T);

    Interval eps_interval_x(-eps_x,eps_x);
    Interval eps_interval_y(-eps_y,eps_y);
    Interval E(-E_eps,E_eps);
    IVector u{  eps_interval_x,
                0,
                eps_interval_y,
                eps_interval_y,
                eps_interval_y, E};

    IVector u0{0,0,0,0,0,E};

    IVector v = IVector(v0) + T_total * u;
    v[3] = 0.;
    v[5] = 0.;

    IMatrix D(6,6);
    C1HORect2Set S(v);
    IVector y = Ivf.pm_z(S,D);
    D = Ivf.pm_z.computeDP(y,D);
    
    IVector w{y[1],y[3]};
    IMatrix M({{D[1][0],D[1][4]},{D[3][0],D[3][4]}});
    
    IVector v1 = IVector(v0) + T_total * u0;

    C0HORect2Set S0(v1);
    y = Ivf.pm_z(S0);
    IVector w0{y[1],y[3]};

    IVector N = - matrixAlgorithms::gauss(M,w0);
    
    IVector X_ext = T_total * u;
    IVector X{X_ext[0], X_ext[4]};
    
    cout << X << endl;
    cout << N << endl;

    return subset(N,X);
}

Interval cone_coeff(long double eps_x, long double eps_y, long double E_eps) {
    LDMatrix T = get_change_of_basis();
    
    IMatrix T_total(T_E * T);
    IMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    Interval eps_interval_x(0,eps_x);
    Interval eps_interval_y(-eps_y,eps_y);
    Interval E(-E_eps,E_eps);

    int division_E = 10;
    int division_x = 1000;

    long double E_delta = 2 * E_eps / division_E;
    long double E_begin = -E_eps;
    long double E_end = -E_eps + E_delta;
    Interval E_i(E_begin,E_end);

    IMatrix D_total(5,5);

    for(int i = 0; i < division_E; i++) {
        long double x_delta = eps_x / division_x;
        Interval X_j(0,x_delta);

        for(int j = 0; j < division_x; j++) {
            IVector u{eps_interval_x,0,eps_interval_y,eps_interval_y,eps_interval_y, E_i};
            C1HORect2Set S(C1Rect2Set::C0BaseSet(IVector(v0),T_total,u), C1Rect2Set::C1BaseSet(T_total));
            IMatrix D(6,6);
            IVector y = Ivf.pm_y(S,D,2);
            D = Ivf.pm_y.computeDP(y,D);
            D = T_total_inv * D;
            IMatrix D_reduced = matrix_erase_cord(D,1);

            if(i == 0 && j == 0) D_total = D_reduced;
            else D_total = intervalHull(D_reduced,D_total);
            X_j += x_delta;
        }
        // if(i % 50 == 0) cout << i << endl;
        cout << i << endl;

        E_i += E_delta;
    }
    // [0.0002172074019, 0.0002172075871]

    // IVector u{  eps_interval_x,
    //             0,
    //             eps_interval_y,
    //             eps_interval_y,
    //             eps_interval_y, E};

    
    // C1HORect2Set S(C1Rect2Set::C0BaseSet(IVector(v0),T_total,u), C1Rect2Set::C1BaseSet(T_total));

    IMaxNorm max_norm;
    ISumNorm sum_norm;

    // IMatrix D(6,6);
    // IVector y = Ivf.pm_y(S,D,2);
    // D = Ivf.pm_y.computeDP(y,D);
    // D = T_total_inv * D;
    
    // IMatrix D_reduced = matrix_erase_cord(D,1);

    Interval A = abs(D_total[0][0]).left();
    Interval B = max_norm(IVector{D_total[0][1],D_total[0][2],D_total[0][3]}).right();
    Interval C = sum_norm(IVector{D_total[1][0],D_total[2][0],D_total[3][0]}).right();

    IMatrix H1 = matrix_erase_cord(D_total,4);
    Interval H = sum_norm(matrix_erase_cord(H1,0)).right();

    Interval delta = power(-A + H,2) - 4 * B * C;
    Interval alpha_1 = ( A - H + sqrt(delta) ) / (2 * B);
    Interval alpha_2 = ( A - H - sqrt(delta) ) / (2 * B);

    cout << alpha_1 << endl;
    cout << alpha_2 << endl;
    
    return 1;
}


IVector quick_eval(IVector W, IMatrix T) {
    C0HORect2Set S1(IVector(v0), T, W);
    IVector W1 = Ivf.pm_x(S1);
    C0HORect2Set S2(W1);
    IVector W2 = Ivf.pm_y(S2);
    return IVector{W2[3],W2[5]};
}

void test(Interval E) {
    long double alpha = 0.1103049783;
    // long double alpha = 0.;
    
    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    
    // cout << 9.199674025e-08 * 0.1103049783 << endl;

    LDMatrix T = get_change_of_basis();
    
    IMatrix T_total(T_E * T);
    LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    

    LDMatrix D = LDMatrix::Identity(6), Phi(6,6);
    LDVector u = v0;
    for(int i = 0; i < 2; i++) {
        u = vf.pm_y(u,Phi);
        D = vf.pm_y.computeDP(u,Phi) * D;
    }

    LDVector w_new = T_total_inv * (w0 - v0);
    // cout << w_new << endl;
    // cout << w_new << endl;
    
    // cout << W << endl;
    // cout << quick_eval(W,T_total) << endl;
    // IVector W1(W),W2(W),W3(W),W4(W);
    
    double x0 = w_new[0];
    cout << x0 << endl;
    double x_eps = 1e-8;
    double y_eps = (x0 - x_eps) * alpha; 
    y_eps = 1e-5;

    IVector W{x0,0,0,0,0,E};
    IVector W1{W[0] - x_eps, 0, 0, 0, 0, E};
    IVector W2{W[0] + x_eps, 0, 0, 0, 0, E};
    IVector W3{W[0], -y_eps, -y_eps, -y_eps, -y_eps, E};
    IVector W4{W[0], y_eps, y_eps, y_eps, y_eps, E};

    cout << quick_eval(W1,T_total) << endl;
    cout << quick_eval(W2,T_total) << endl;
    cout << quick_eval(W3,T_total) << endl;
    cout << quick_eval(W4,T_total) << endl;
    


    // cout << W1 << endl;
    // cout << W2 << endl;
    // cout << quick_eval(W1, T_total) << endl;
    // cout << quick_eval(W2, T_total) << endl;
    

    // return;
    
    // long double x_eps = 5e-9;
    // long double x0 = w_new[0];
    // long double x_minus = x0 - x_eps;
    // long double x_plus = x0 + x_eps;

    // // cout << x_minus << endl;
    // // cout << x_plus << endl;

    // long double y_lower_bound = -alpha * x_minus;
    // long double y_upper_bound = alpha * x_minus;
    
    // Interval x_bound(x_minus,x_plus);
    // Interval y_bound(y_lower_bound, y_upper_bound);
    // IVector left_cover{x_bound.left(),0,y_bound, y_bound, y_bound, E};
    // IVector right_cover{x_bound.right(),0,y_bound, y_bound, y_bound, E};
    // IVector lower_cover{x_bound,0,y_bound.left(),y_bound.left(),y_bound.left(),E};
    // IVector upper_cover{x_bound,0,y_bound.right(),y_bound.right(),y_bound.right(),E};

    // string file_name = "data.txt";
    // ofstream file(file_name);
    // file.precision(15);

    // int x_division = 2000;
    // long double x_delta = 2 * x_eps / x_division;
    // interval X_i(0, x_delta);
    // X_i += x_bound.left();

    // for(int i = 0; i < x_division; i++) {
        
    //     IVector lower_i{X_i,0,y_bound.left(),y_bound.left(),y_bound.left(),E};
    //     IVector upper_i{X_i,0,y_bound.right(),y_bound.right(),y_bound.right(),E};

    //     cout << lower_i << endl;
    //     cout << upper_i << endl;

    //     cout << quick_eval(lower_i,T_total) << endl;
    //     cout << quick_eval(upper_i,T_total) << endl;
    //     cout << "XD" << endl;

    //     break;

    //     C0HORect2Set S_lower(C1Rect2Set::C0BaseSet(IVector(v0),T_total,lower_i));
    //     C0HORect2Set S_upper(C1Rect2Set::C0BaseSet(IVector(v0),T_total,upper_i));

    //     IVector y_lower = Ivf.pm_x(S_lower);
    //     IVector y_upper = Ivf.pm_x(S_upper);

    //     C0HORect2Set S1_lower(y_lower);
    //     C0HORect2Set S1_upper(y_upper);

    //     IVector z_lower = Ivf.pm_y(S1_lower);
    //     IVector z_upper = Ivf.pm_y(S1_upper);

    //     file << z_lower[3] << " " << z_lower[5] << endl;
    //     file << z_upper[3] << " " << z_upper[5] << endl;

    //     X_i += x_delta;
    // }

    // file.close();
}