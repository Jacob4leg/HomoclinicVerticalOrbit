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


bool prove_fixed_point(LDVector v_E, Interval X, Interval Y, Interval E) {
    
    LDMatrix T = get_change_of_basis();
    IMatrix T_total = IMatrix(T_E * T);

    IVector w{X,0,Y,Y,Y, E};
    IVector u0(v_E);
    

    IVector u = IVector(v_E) + T_total * w;
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
    
    // cout << Z << endl;
    // cout << N << endl;
    // cout << diam(N) << endl;
    
    return subset(N,Z);
}

Interval cone_coeff(long double eps_x, long double eps_y, Interval E) {
    LDMatrix T = get_change_of_basis();
    
    IMatrix T_total(T_E * T);
    IMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    Interval eps_interval_x(0,eps_x);
    Interval eps_interval_y(-eps_y,eps_y);
    // Interval E(-E_eps,E_eps);

    int division_E = 1;
    int division_x = 10;

    long double E_delta = (E.rightBound() - E.leftBound()) / division_E;
    long double E_begin = E.leftBound();
    long double E_end = E.leftBound() + E_delta;
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

void rectangle() {


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
    
    
    LDMatrix T = get_change_of_basis();
    LDMatrix T_total(T_E * T);
    LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    LDVector w_new = T_total_inv * (w0 - v0);

    long double x0 = w_new[0];
    long double x_radius = 1e-8;
    
    Interval X(x0 - x_radius, x0 + x_radius);
    Interval E(-1e-4,1e-4);

    int x_div = 3000;
    int E_div = 3000;

    ofstream file("rectangle_lin_fixed_point.txt");

    Interval x_delta(0,(X.rightBound() - X.leftBound()) / x_div);
    Interval E_delta(0,(E.rightBound() - E.leftBound()) / E_div);
    
    Interval x_i = X.left() + x_delta;
    Interval E_i = E.left() + E_delta;
    
    Interval E_left = E.left();
    Interval E_right = E.right();
    
    LDVector w{0,0,0,0,0,E_left.leftBound()};
    LDVector v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    IVector V_E(v_E);

    Interval X_test(-1e-13,1e-13);
    Interval Y_test(-1e-13,1e-13);
    IVector W_test{X_test,0,Y_test,Y_test,Y_test,0};

    if(!prove_fixed_point(v_E,X_test,Y_test,E_left)) {
        throw runtime_error("Fixed point not proved.");
    }

    for(int i = 0; i < x_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        C0HORect2Set S1(V_E,T_total, IVector{x_i,0,0,0,0,0} + W_test);

        IVector u1 = Ivf.pm_x(S1);

        C0HORect2Set S2(u1);
        IVector u2 = Ivf.pm_y(S2);

        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;
        x_i += x_delta.right();
    }

    x_i = X.left() + x_delta;
    w = LDVector{0,0,0,0,0,E_right.leftBound()};
    v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
    V_E = IVector(v_E);

    if(!prove_fixed_point(v_E,X_test,Y_test,E_right)) {
        throw runtime_error("Fixed point not proved.");
    }

    for(int i = 0; i < x_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        C0HORect2Set S1(V_E,T_total, IVector{x_i,0,0,0,0,0} + W_test);
        IVector u1 = Ivf.pm_x(S1);

        C0HORect2Set S2(u1);
        IVector u2 = Ivf.pm_y(S2);

        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;
        x_i += x_delta.right();
    }

    Interval x_left = X.left();
    Interval x_right = X.right();

    for(int i = 0; i < E_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        w = LDVector{0,0,0,0,0,E_i.leftBound()};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        V_E = IVector(v_E);

        if(!prove_fixed_point(v_E,X_test,Y_test,E_i)) {
            throw runtime_error("Fixed point not proved.");
        }

        C0HORect2Set S1(V_E,T_total, IVector{x_left,0,0,0,0,0} + W_test);
        IVector u1 = Ivf.pm_x(S1);

        C0HORect2Set S2(u1);
        IVector u2 = Ivf.pm_y(S2);
        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;

        E_i += E_delta.right();
    }

    E_i = E.left() + E_delta;

    for(int i = 0; i < E_div; i++) {
        if(i % 50 == 0) cout << i << endl;

        w = LDVector{0,0,0,0,0,E_i.leftBound()};
        v_E = vf.findVerticalLyapunovOrbit(T_total * w + v0);
        V_E = IVector(v_E);

        if(!prove_fixed_point(v_E,X_test,Y_test,E_i)) {
            throw runtime_error("Fixed point not proved.");
        }

        C0HORect2Set S1(V_E,T_total, IVector{x_right,0,0,0,0,0} + W_test);
        IVector u1 = Ivf.pm_x(S1);

        C0HORect2Set S2(u1);
        IVector u2 = Ivf.pm_y(S2);
        file << u2[3].leftBound() << " " << u2[3].rightBound() << " " << u2[5].leftBound() << " " << u2[5].rightBound() << endl;

        E_i += E_delta.right();
    }



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


    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    
    LDMatrix T = get_change_of_basis();
    IMatrix T_total(T_E * T);
    LDMatrix T_total_inv(matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(T_E * T,1)),1));

    LDVector w_new = T_total_inv * (w0 - v0);
    long double x0 = w_new[0];

    // Interval E_left = E.left();
    // Interval E_right = E.right();
    // Interval alpha_1 = cone_coeff(x_eps,y_eps, E_left);
    // Interval alpha_2 = cone_coeff(x_eps,y_eps,E_right);
    

    long double x_delta = 2e-8;
    Interval X(x0 - x_delta, x0 + x_delta);
    
    int x_div = 10000;
    int E_div = 10000;

    // ofstream file("data1.txt");
    Interval alpha = cone_coeff(x_eps, y_eps, E.left());
    // cout << alpha << endl;


    Interval X_delta(0, (X.rightBound() - X.leftBound()) / x_div );
    Interval X_i = X.left() + X_delta;
    Interval E_delta(0, (E.rightBound() - E.leftBound()) / E_div);
    Interval E_i = E.left() + E_delta;

    for(int i = 0; i < x_div; i++) {
        
        Interval Y( -(X_i * alpha).leftBound(), (X_i * alpha).rightBound() );
        Y = 0;

        IVector W1{X_i,0,Y,Y,Y,E.left()};
        IVector W2{X_i,0,Y,Y,Y,E.right()};
        auto res1 = quick_eval(W1,T_total);
        auto res2 = quick_eval(W2,T_total);
        
        // file << res1[0].leftBound() << " " << res1[0].rightBound() << " " << res1[1].leftBound() << " " << res1[1].rightBound() << endl;
        cout << res2[0].leftBound() << " " << res2[0].rightBound() << " " << res2[1].leftBound() << " " << res2[1].rightBound() << endl;
        X_i += X_delta.right();
    }
    
    alpha = cone_coeff(x_eps, y_eps, E.left() + 100 * E_delta);
    for(int i = 0; i < E_div; i++) {
        // if(i % 100 == 0) {
        //     alpha = cone_coeff(x_eps, y_eps, E_i + 100 * E_delta);
        //     cout << alpha << endl;
        // }
        
        Interval Y1( -(X.left() * alpha).leftBound(), (X.left() * alpha).rightBound());
        Interval Y2( -(X.right() * alpha).leftBound(), (X.right() * alpha).rightBound());
        
        Y1 = 0;
        Y2 = 0;

        IVector W1{X.left(),0,Y1,Y1,Y1,E_i};
        IVector W2{X.right(),0,Y2,Y2,Y2,E_i};
        auto res1 = quick_eval(W1,T_total);
        auto res2 = quick_eval(W2,T_total);

        // file << res1[0].leftBound() << " " << res1[0].rightBound() << " " << res1[1].leftBound() << " " << res1[1].rightBound() << endl;
        // file << res2[0].leftBound() << " " << res2[0].rightBound() << " " << res2[1].leftBound() << " " << res2[1].rightBound() << endl;
        E_i += E_delta.right();
    }


    // file.close();

    // Interval Y1 = (x0 - x_delta) * alpha_1; // x_left x E_left
    // Interval Y2 = (x0 - x_delta) * alpha_2; // x_left x E_right
    // Interval Y3 = (x0 + x_delta) * alpha_1; // x_right x E_left
    // Interval Y4 = (x0 + x_delta) * alpha_2; // x_right x E_right

    // IVector W1{x0 - x_delta, 0, Y1, Y1, Y1, E_left};
    // IVector W2{x0 - x_delta, 0, Y2, Y2, Y2, E_right};
    // IVector W3{x0 + x_delta, 0, Y3, Y3, Y3, E_left};
    // IVector W4{x0 + x_delta, 0, Y4, Y4, Y4, E_right};

    // cout << quick_eval(W1,T_total) << endl;
    // cout << quick_eval(W2,T_total) << endl;
    // cout << quick_eval(W3,T_total) << endl;
    // cout << quick_eval(W4,T_total) << endl;
    
}

void test(Interval E) {
    long double alpha = 0.1103049783L;
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
    
    double x0 = w_new[0];
    cout << x0 << endl;
    double x_eps = 5e-8;
    double y_eps = (x0 - x_eps) * alpha; 
    Interval Y(-y_eps,y_eps);
    // y_eps = 1e-5;

    cout << x0 << endl;
    IVector W{x0,0,0,0,0,0};
    IVector W1{W[0] - x_eps, 0, Y, Y, Y, E.left()};
    IVector W2{W[0] + x_eps, 0, Y, Y, Y, E.left()};
    IVector W3{W[0] - x_eps, 0, Y, Y, Y, E.right()};
    IVector W4{W[0] + x_eps, 0, Y, Y, Y, E.right()};

    cout << quick_eval(W,T_total) << endl;
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