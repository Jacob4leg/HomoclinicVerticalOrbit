#include "prove_fixed_point.h"
#include<chrono>
#include<fstream>

using namespace std;
using namespace capd;

CR3BP<long double> vf;
ICR3BP Ivf;

long double z = 0.05316795030019990465013731L;
LDVector v0 = vf.findVerticalLyapunovOrbit(LDVector({0.94690970780356629,0,z,0,-0.011166782657715382,0}));
// LDMatrix B{ {1,0,0,0,0,0},
//             {0,1,0,0,0,0},
//             {0,0,1,0,0,0},
//             {0,0,0,1008.04179513,0,85.38967319},
//             {0,0,0,0,1,0},
//             {0,0,0,-188.97139923,0,-112.39859377}};
// LDMatrix B_inv = matrixAlgorithms::gaussInverseMatrix(B);
long double E0 = vf.E(v0)[0];

LDMatrix T0 = energy_change_of_basis(v0,vf.E);

tuple<LDMatrix,LDMatrix> get_change_of_basis(LDVector v_E) {
    LDMatrix T_E = energy_change_of_basis(v_E,vf.E);
    LDMatrix D = LDMatrix::Identity(6), T(6,6);
    LDVector u = v_E;
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
    return {T_E,matrix_add_cord(L,1)};
}


bool prove_fixed_point(LDVector v_E, Interval X, Interval Y, Interval E, const LDMatrix &T) {
    
    // LDMatrix T = get_change_of_basis();
    IMatrix T_total = IMatrix(T);

    IVector w{X,0,Y,Y,Y, E};
    IVector u0(v_E);
    

    IVector u = u0 + T_total * w;
    // v[3] = 0.;
    // v[5] = 0.;

    IMatrix D(6,6);
    C1HORect2Set S(u);
    IVector y = Ivf.pm_z(S,D);
    D = Ivf.pm_z.computeDP(y,D);
    
    // IVector w{y[1],y[3]};
    IMatrix M({{D[1][0],D[1][4]},{D[3][0],D[3][4]}});
    

    C0HORect2Set S0(u0);
    y = Ivf.pm_z(S0);
    IVector w0{y[1],y[3]};

    IVector N = - matrixAlgorithms::gauss(M,w0);
    
    IVector Z_ext = T_total * u;
    IVector Z{u[0] - u0[0], u[4] - u0[4]};

    // cout << "N=" << N << endl;
    // cout << "Z=" << Z << endl;
    // cout << endl;
    return subset(N,Z);
}

Interval cone_coeff(int division_x, int division_E, Interval E, const LDVector &v_E, const LDMatrix &T) {
    
    long double eps_x = 1e-7;
    long double eps_y = 1e-8;

    IMatrix T_total(T);
    IMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T,1)),1));

    Interval eps_interval_x(0,eps_x);
    Interval eps_interval_y(-eps_y,eps_y);
    // Interval E(-E_eps,E_eps);

    long double x_delta = eps_x / division_x;
    long double E_delta = (E.rightBound() - E.leftBound()) / division_E;
    long double E_begin = E.leftBound();
    long double E_end = E.leftBound() + E_delta;
    Interval E_i(E_begin,E_end);

    IMatrix D_total(5,5);

    for(int i = 0; i < division_E; i++) {
        
        Interval X_j(0,x_delta);
        
        for(int j = 0; j < division_x; j++) {
            
            IVector u{X_j,0,eps_interval_y,eps_interval_y,eps_interval_y, E_i};
            C1HORect2Set S(C1Rect2Set::C0BaseSet(IVector(v_E),T_total,u), C1Rect2Set::C1BaseSet(T_total));
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
        // cout << i << endl;

        E_i += E_delta;
    }
    // [0.0002172074019, 0.0002172075871]

    // IVector u{  eps_interval_x,
    //             0,
    //             eps_interval_y,
    //             eps_interval_y,
    //             eps_interval_y, E};

    
    IMaxNorm max_norm;
    ISumNorm sum_norm;

    Interval A = abs(D_total[0][0]).left();
    Interval B = sum_norm(IVector{D_total[0][1],D_total[0][2],D_total[0][3]}).right();
    Interval C = max_norm(IVector{D_total[1][0],D_total[2][0],D_total[3][0]}).right();

    IMatrix H1 = matrix_erase_cord(D_total,4);
    Interval H = max_norm(matrix_erase_cord(H1,0)).right();

    Interval delta = power(-A + H,2) - 4 * B * C;
    Interval alpha_1 = ( A - H + sqrt(delta) ) / (2 * B);
    Interval alpha_2 = ( A - H - sqrt(delta) ) / (2 * B);

    // cout << alpha_1 << endl;
    // cout << alpha_2 << endl;
    
    return alpha_2.right();
}


IVector quick_eval(IVector W, IMatrix T) {
    C0HORect2Set S1(IVector(v0), T, W);
    IVector W1 = Ivf.pm_x(S1);
    C0HORect2Set S2(W1);
    IVector W2 = Ivf.pm_y(S2);
    return IVector{W2[3],W2[5]};
}

LDVector quick_eval_non_rig(LDVector w, const LDMatrix &T, const LDMatrix &T_inv) {
    LDVector u = T * w + v0;
    LDVector u1 = vf.pm_x(u);
    LDVector u2 = vf.pm_y(u1);
    return LDVector{u2[3],u2[5]};
}


tuple<LDMatrix,LDMatrix,LDMatrix,LDMatrix> get_poincare_changes_of_basis(LDMatrix T_total_0, long double x0, long double x_delta, long double E_delta) {
    LDVector v_E = vf.findVerticalLyapunovOrbit(v0 + T0 * LDVector{0,0,0,0,0,E_delta});
    
    auto Q = get_change_of_basis(v_E);
    LDMatrix T_total_E(get<0>(Q) * get<1>(Q));

    LDVector w0 = v0 + T_total_0 * LDVector{x0,0,0,0,0,0};

    LDMatrix D(6,6);
    LDVector u = vf.pm_x(w0,D);
    D = vf.pm_x.computeDP(u,D);

    // cout << D << endl;

    // LDVector rV(6), iV(6);
    // LDMatrix rVec(6,6), iVec(6,6);
    // computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);
    D = us_change_of_basis(D);
    // cout << D << endl;

    LDVector w1 = v0 + T_total_0 * LDVector{x0 + x_delta,0,0,0,0,0};
    LDVector w2 = v_E + T_total_E * LDVector{x0,0,0,0,0,0};

    LDVector Pw1_x = vf.pm_x(w1);
    LDVector Pw2_x = vf.pm_x(w2);
    Pw1_x[0] = 0.; Pw2_x[0] = 0.;

    LDSumNorm norm;

    LDMatrix Phi_X = LDMatrix::Identity(6);
    Phi_X = D;
    // Phi_X.column(0) = LDVector{1,0,0,0,0,0};
    Phi_X.row(0) = LDVector{1,0,0,0,0,0};
    Phi_X.column(0) = u / norm(u);

    // Phi_X.column(1) = Pw1_x;
    // Phi_X.column(2) = Pw2_x;
    

    LDVector Pw1_y = vf.pm_y(Pw1_x);
    LDVector Pw2_y = vf.pm_y(Pw2_x);

    LDMatrix Phi_Y = LDMatrix::Identity(6);

    Phi_Y[3][3] = Pw1_y[3]; Phi_Y[3][5] = Pw2_y[3];
    Phi_Y[5][3] = Pw1_y[5]; Phi_Y[5][5] = Pw2_y[5];

    return {matrixAlgorithms::gaussInverseMatrix(Phi_X),Phi_X,matrixAlgorithms::gaussInverseMatrix(Phi_Y),Phi_Y};
}

void rectangle_non_rig() {
    LDMatrix O1(6,6), O2(6,6), O(6,6);
    LDVector v1 = vf.pm_y(v0,O1);
    LDVector v2 = vf.pm_y(v1,O2);

    O1 = vf.pm_y.computeDP(v1,O1);
    O2 = vf.pm_y.computeDP(v2,O2);

    O = O2 * O1;
    auto [lambda, v_eig] = get_dominant_eigenvalue_and_eigenvector(O);

    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    
    auto [T_E, T] = get_change_of_basis(v0);
    LDMatrix T_total(T_E * T);
    LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    LDVector w_new = T_total_inv * (w0 - v0);

    long double x0 = w_new[0];
    long double x_radius = 1e-8;
    long double E_radius = 1e-4;

    auto [B,B_inv,C,C_inv] = get_poincare_changes_of_basis(T_total,x0,x_radius,E_radius);

    Interval X(x0 - x_radius, x0 + x_radius);
    Interval E(-E_radius,E_radius);

    int x_div = 300;
    int E_div = 300;

    ofstream file("rectangle_non_rig.txt");

    long double x_delta = 2 * x_radius / x_div;
    long double E_delta = 2 * E_radius / E_div;
    
    long double x_i = X.leftBound();
    long double E_i = E.leftBound();

    LDVector w{0,0,0,0,0,-E_radius};
    LDVector v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    T_total = get<0>(get_change_of_basis(v_E)) * get<1>(get_change_of_basis(v_E));

    // A = LDMatrix::Identity(2);
    // LDMatrix A_inv = matrixAlgorithms::gaussInverseMatrix(A);
    
    for(int i = 0; i < x_div; i++) {
        LDVector v_i = v_E + T_total * LDVector{x_i,0,0,0,0,0};
        LDVector v1 = vf.pm_x(v_i);
        LDVector v2 = vf.pm_y(v1);

        auto res = v2;

        file << res[3] << " " << res[5] << endl;
        x_i += x_delta;
    }

    w = LDVector{0,0,0,0,0,E_radius};
    v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    T_total = get<0>(get_change_of_basis(v_E)) * get<1>(get_change_of_basis(v_E));
    x_i = X.leftBound();

    for(int i = 0; i < x_div; i++) {
        LDVector v_i = v_E + T_total * LDVector{x_i,0,0,0,0,0};
        LDVector v1 = vf.pm_x(v_i);
        LDVector v2 = vf.pm_y(v1);

        auto res = v2;

        file << res[3] << " " << res[5] << endl;
        x_i += x_delta;
    }

    for(int i = 0; i < E_div; i++) {
        w = LDVector{0,0,0,0,0,E_i};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        T_total = get<0>(get_change_of_basis(v_E)) * get<1>(get_change_of_basis(v_E));
        LDVector v_i = v_E + T_total * LDVector{X.leftBound(),0,0,0,0,0};
        LDVector v1 = vf.pm_x(v_i);
        LDVector v2 = vf.pm_y(v1);

        auto res = v2;

        file << res[3] << " " << res[5] << endl;
        E_i += E_delta;
    }
    
    E_i = E.leftBound();

    for(int i = 0; i < E_div; i++) {
        w = LDVector{0,0,0,0,0,E_i};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        T_total = get<0>(get_change_of_basis(v_E)) * get<1>(get_change_of_basis(v_E));
        LDVector v_i = v_E + T_total * LDVector{X.rightBound(),0,0,0,0,0};
        LDVector v1 = vf.pm_x(v_i);
        LDVector v2 = vf.pm_y(v1);
        
        auto res = v2;

        file << res[3] << " " << res[5] << endl;
        E_i += E_delta;
    }
}

void rectangle() {
    // rectangle_non_rig();
    // return;

    LDMatrix O1(6,6), O2(6,6), O(6,6);
    LDVector v1 = vf.pm_y(v0,O1);
    LDVector v2 = vf.pm_y(v1,O2);

    O1 = vf.pm_y.computeDP(v1,O1);
    O2 = vf.pm_y.computeDP(v2,O2);

    O = O2 * O1;

    LDVector rV(6), iV(6);
    LDMatrix rVec(6,6), iVec(6,6);
    computeEigenvaluesAndEigenvectors(O, rV, iV, rVec, iVec);
    LDVector v_eig = rVec.column(0);

    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    
    auto [T_E, T] = get_change_of_basis(v0);
    LDMatrix T_total(T_E * T);
    LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    cout << T_total_inv * O * T_total << endl;
    return;

    LDVector w_new = T_total_inv * (w0 - v0);

    long double x0 = w_new[0];
    long double x_to_change_of_basis = 1e-10;
    long double E_to_change_of_basis = 1e-5;

    long double x_radius = 1e-10;
    long double E_radius = 1e-7;

    auto [B,B_inv,C,C_inv] = get_poincare_changes_of_basis(T_total,x0,x_radius,E_radius);

    // auto [B,B_inv,C,C_inv] = get_poincare_changes_of_basis(T_total,x0,x_to_change_of_basis,E_to_change_of_basis);

    // return;
    
    
    Interval X(x0 - x_radius, x0 + x_radius);
    Interval E(-E_radius,E_radius);
        
    int x_div = 1000;
    int E_div = 50000;

    ofstream file("rectangle_lin_cone.txt");

    Interval x_delta(0,(X.rightBound() - X.leftBound()) / x_div);
    Interval E_delta(0,(E.rightBound() - E.leftBound()) / E_div);
    
    Interval x_i = X.left() + x_delta;
    Interval E_i = E.left() + E_delta;
    
    Interval E_left = E.left();
    Interval E_right = E.right();
    
    LDVector w{0,0,0,0,0,E_left.leftBound()};
    LDVector v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    IVector V_E(v_E);
    auto T_tuple = get_change_of_basis(v_E);
    T_E = get<0>(T_tuple);
    T = get<1>(T_tuple);
    T_total = T_E * T;

    Interval X_test(-1e-13,1e-13);
    Interval Y_test(-1e-13,1e-13);
    IVector W_test{X_test,0,Y_test,Y_test,Y_test,0};

    Interval unit_interval(-1,1);
    IVector unit_center_stable_cube{0,0,unit_interval,unit_interval,unit_interval,0};

    if(!prove_fixed_point(v_E,X_test,Y_test,Interval(-1e-13,1e-13),T_total)) {
        throw runtime_error("Fixed point not proved.");
    }

    Interval returnTime;
    IVector zero_vector{0,0,0,0,0,0};
    IVector PW0(vf.pm_x(w0));

    

    Interval alpha = cone_coeff(20,1,Interval(-1e-13,1e-13),v_E,T_total);
    cout << alpha << endl;
    
    chrono::steady_clock::time_point beg = chrono::steady_clock::now();
    chrono::steady_clock::time_point end;

    for(int i = 0; i < x_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        C0HOTripletonSet S1(V_E,T_total, IVector{x_i,0,0,0,0,0} + W_test + unit_center_stable_cube * x_i * alpha);
        Ivf.pm_x(S1);
        IVector u2 = Ivf.pm_y(S1, zero_vector, C, returnTime);
        

        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;
        x_i += x_delta.right();
    }
    
    end = chrono::steady_clock::now();
    cout << "First edge calculated after " << chrono::duration_cast<chrono::minutes>(end - beg).count() << " min " 
                                            << chrono::duration_cast<chrono::seconds>(end - beg).count() % 60 << " sec." << endl;

    x_i = X.left() + x_delta;
    w = LDVector{0,0,0,0,0,E_right.leftBound()};
    v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    V_E = IVector(v_E);
    T_tuple = get_change_of_basis(v_E);
    T_E = get<0>(T_tuple);
    T = get<1>(T_tuple);
    T_total = T_E * T;
    
    if(!prove_fixed_point(v_E,X_test,Y_test,E_right,T_total)) {
        throw runtime_error("Fixed point not proved.");
    }

    alpha = cone_coeff(20,1,0,v_E,T_total);
    cout << alpha << endl;

    for(int i = 0; i < x_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        C0HOTripletonSet S1(V_E,T_total, IVector{x_i,0,0,0,0,0} + W_test + unit_center_stable_cube * x_i * alpha);
        Ivf.pm_x(S1);
        IVector u2 = Ivf.pm_y(S1, zero_vector, C, returnTime);

        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;
        x_i += x_delta.right();
    }

    chrono::steady_clock::time_point beg1 = chrono::steady_clock::now();
    end = chrono::steady_clock::now();
    cout << "Two edges calculated after " << chrono::duration_cast<chrono::minutes>(end - beg).count() << " min " 
                                            << chrono::duration_cast<chrono::seconds>(end - beg).count() % 60 << " sec." << endl;

    Interval x_left = X.left();
    Interval x_right = X.right();

    for(int i = 0; i < E_div; i++) {
        if(i % 500 == 0 && i > 0) {
            end = chrono::steady_clock::now();
            cout << "Third edge calculated with " << (i * 100) / E_div << "% progress." << 
            "Time elapsed: " << chrono::duration_cast<chrono::minutes>(end - beg1).count() << " min " 
                                << chrono::duration_cast<chrono::seconds>(end - beg1).count() % 60 << " sec." <<  endl;
        } 
        
        w = LDVector{0,0,0,0,0,E_i.leftBound()};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        V_E = IVector(v_E);
        T_tuple = get_change_of_basis(v_E);
        T_E = get<0>(T_tuple);
        T = get<1>(T_tuple);
        T_total = T_E * T;
        
        if(!prove_fixed_point(v_E,X_test,Y_test,E_delta,T_total)) {
            throw runtime_error("Fixed point not proved.");
        }

        alpha = cone_coeff(1,1,E_delta,v_E,T_total);

        C0HOTripletonSet S1(V_E,T_total, IVector{x_left,0,0,0,0,E_delta} + W_test + unit_center_stable_cube * x_left * alpha);
        Ivf.pm_x(S1);
        IVector u2 = Ivf.pm_y(S1, zero_vector, C, returnTime);

        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;

        E_i += E_delta.right();
    }

    E_i = E.left() + E_delta;

    beg1 = chrono::steady_clock::now();
    end = chrono::steady_clock::now();
    cout << "Three edges calculated after " << chrono::duration_cast<chrono::minutes>(end - beg).count() << " min " 
                                            << chrono::duration_cast<chrono::seconds>(end - beg).count() % 60 << " sec." << endl;

    for(int i = 0; i < E_div; i++) {
        if(i % 500 == 0 && i > 0) {
            end = chrono::steady_clock::now();
            cout << "Fourth edge calculated with " << (i * 100) / E_div << "% progress." << 
            "Time elapsed: " << chrono::duration_cast<chrono::minutes>(end - beg1).count() << " min " 
                                << chrono::duration_cast<chrono::seconds>(end - beg1).count() % 60 << " sec." <<  endl;
        }

        w = LDVector{0,0,0,0,0,E_i.leftBound()};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        V_E = IVector(v_E);
        T_tuple = get_change_of_basis(v_E);
        T_E = get<0>(T_tuple);
        T = get<1>(T_tuple);
        T_total = T_E * T;

        if(!prove_fixed_point(v_E,X_test,Y_test,E_delta,T_total)) {
            throw runtime_error("Fixed point not proved.");
        }
        
        alpha = cone_coeff(1,1,E_delta,v_E,T_total);

        C0HOTripletonSet S1(V_E,T_total, IVector{x_right,0,0,0,0,E_delta} + W_test + unit_center_stable_cube * x_right * alpha);
        Ivf.pm_x(S1);
        IVector u2 = Ivf.pm_y(S1, zero_vector, C, returnTime);
        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;

        E_i += E_delta.right();
    }

    end = chrono::steady_clock::now();
    cout << "All edges calculated after " << chrono::duration_cast<chrono::minutes>(end - beg).count() << " min " 
                                            << chrono::duration_cast<chrono::seconds>(end - beg).count() % 60 << " sec." << endl;

    file.close();

    // return;

    // for(int i = 0; i < x_div; i++) {

    //     LDVector w1{x_i,0,0,0,0,E.leftBound()};
    //     LDVector w2{x_i,0,0,0,0,E.rightBound()};
    //     LDVector w3{X.leftBound(),0,0,0,0,E_i};
    //     LDVector w4{X.rightBound(),0,0,0,0,E_i};

    //     // auto res1 = DD_inv * LDVector{w1[0],w1[5]};
    //     // auto res2 = DD_inv * LDVector{w2[0],w2[5]};
    //     // auto res3 = DD_inv * LDVector{w3[0],w3[5]};
    //     // auto res4 = DD_inv * LDVector{w4[0],w4[5]};

    //     auto res1 = quick_eval_non_rig(w1,T_total,T_total_inv);
    //     auto res2 = quick_eval_non_rig(w2,T_total,T_total_inv);
    //     auto res3 = quick_eval_non_rig(w3,T_total,T_total_inv);
    //     auto res4 = quick_eval_non_rig(w4,T_total,T_total_inv);

    //     file << x_i << " " << E.leftBound() << " " << res1[0] << " " << res1[1] << endl;
    //     file << x_i << " " << E.rightBound() << " " << res2[0] << " " << res2[1] << endl;
    //     file << X.leftBound() << " " << E_i << " " << res3[0] << " " << res3[1] << endl;
    //     file << X.rightBound() << " " << E_i << " " << res4[0] << " " << res4[1] << endl;

    //     E_i += E_delta;
    //     x_i += x_delta;
        
    // }

    // file.close();
}

void eval_rectangle(long double x_eps, long double y_eps, Interval E) {
    rectangle();
    return;


    // LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
    //     -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    
    // LDMatrix T = get_change_of_basis();
    // IMatrix T_total(T_E * T);
    // LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    // LDVector w_new = T_total_inv * (w0 - v0);
    // long double x0 = w_new[0];

    // // Interval E_left = E.left();
    // // Interval E_right = E.right();
    // // Interval alpha_1 = cone_coeff(x_eps,y_eps, E_left);
    // // Interval alpha_2 = cone_coeff(x_eps,y_eps,E_right);
    

    // long double x_delta = 2e-8;
    // Interval X(x0 - x_delta, x0 + x_delta);
    
    // int x_div = 10000;
    // int E_div = 10000;

    // // ofstream file("data1.txt");
    // Interval alpha = cone_coeff(x_eps, y_eps, E.left());
    // // cout << alpha << endl;


    // Interval X_delta(0, (X.rightBound() - X.leftBound()) / x_div );
    // Interval X_i = X.left() + X_delta;
    // Interval E_delta(0, (E.rightBound() - E.leftBound()) / E_div);
    // Interval E_i = E.left() + E_delta;

    // for(int i = 0; i < x_div; i++) {
        
    //     Interval Y( -(X_i * alpha).leftBound(), (X_i * alpha).rightBound() );
    //     Y = 0;

    //     IVector W1{X_i,0,Y,Y,Y,E.left()};
    //     IVector W2{X_i,0,Y,Y,Y,E.right()};
    //     auto res1 = quick_eval(W1,T_total);
    //     auto res2 = quick_eval(W2,T_total);
        
    //     // file << res1[0].leftBound() << " " << res1[0].rightBound() << " " << res1[1].leftBound() << " " << res1[1].rightBound() << endl;
    //     cout << res2[0].leftBound() << " " << res2[0].rightBound() << " " << res2[1].leftBound() << " " << res2[1].rightBound() << endl;
    //     X_i += X_delta.right();
    // }
    
    // alpha = cone_coeff(x_eps, y_eps, E.left() + 100 * E_delta);
    // for(int i = 0; i < E_div; i++) {
    //     // if(i % 100 == 0) {
    //     //     alpha = cone_coeff(x_eps, y_eps, E_i + 100 * E_delta);
    //     //     cout << alpha << endl;
    //     // }
        
    //     Interval Y1( -(X.left() * alpha).leftBound(), (X.left() * alpha).rightBound());
    //     Interval Y2( -(X.right() * alpha).leftBound(), (X.right() * alpha).rightBound());
        
    //     Y1 = 0;
    //     Y2 = 0;

    //     IVector W1{X.left(),0,Y1,Y1,Y1,E_i};
    //     IVector W2{X.right(),0,Y2,Y2,Y2,E_i};
    //     auto res1 = quick_eval(W1,T_total);
    //     auto res2 = quick_eval(W2,T_total);

    //     // file << res1[0].leftBound() << " " << res1[0].rightBound() << " " << res1[1].leftBound() << " " << res1[1].rightBound() << endl;
    //     // file << res2[0].leftBound() << " " << res2[0].rightBound() << " " << res2[1].leftBound() << " " << res2[1].rightBound() << endl;
    //     E_i += E_delta.right();
    // }


    // // file.close();

    // // Interval Y1 = (x0 - x_delta) * alpha_1; // x_left x E_left
    // // Interval Y2 = (x0 - x_delta) * alpha_2; // x_left x E_right
    // // Interval Y3 = (x0 + x_delta) * alpha_1; // x_right x E_left
    // // Interval Y4 = (x0 + x_delta) * alpha_2; // x_right x E_right

    // // IVector W1{x0 - x_delta, 0, Y1, Y1, Y1, E_left};
    // // IVector W2{x0 - x_delta, 0, Y2, Y2, Y2, E_right};
    // // IVector W3{x0 + x_delta, 0, Y3, Y3, Y3, E_left};
    // // IVector W4{x0 + x_delta, 0, Y4, Y4, Y4, E_right};

    // // cout << quick_eval(W1,T_total) << endl;
    // // cout << quick_eval(W2,T_total) << endl;
    // // cout << quick_eval(W3,T_total) << endl;
    // // cout << quick_eval(W4,T_total) << endl;
    
}