#include<iostream>
#include<future>
#include<numeric>
#include<thread>
#include "capd/capdlib.h"
#include "prove_fixed_point.h"

using namespace std;
using namespace capd;

int main() {
    // CR3BP<long double> vf;
    // LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
    //     -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    
    // ofstream file("homoclinic.txt");
    // LDTimeMap tm(vf.solver);
    // LDTimeMap::SolutionCurve sol(0.);
    // long double T = 20;
    // auto res = tm(T,w0,sol);
    
    // long double h = 0.01;
    // int n = T / h;

    // for(int i = 0; i < n; i++) {
    //     long double x = h * i;
    //     LDVector u = sol(x);
    //     file << u[0] << " " << u[1] << " " << u[2] << endl;
    // }
    // file.close();

    cout.precision(15);

    // long double eps = 1e-6;
    // long double E_eps = 1e-8;

    // long double eps_x = 5e-10;
    // long double eps_y = 1e-8;
    // test();
    // return 0;

    // long double eps_x = 15e-8;
    // long double eps_y = 1.2e-8;
    // long double E_eps = 1e-8;

    long double eps_x = 3e-10;
    long double eps_y = 1e-12;
    long double E_eps = 1e-13;

    // cout << eps_x << endl;
    // cout << eps_y << endl;
    // cout << E_eps << endl;

    Interval X(-eps_x,eps_x);
    Interval Y(-eps_y,eps_y);
    Interval E(-E_eps,E_eps);
    // cout << Y << endl;
    // cout << prove_fixed_point(X,Y,E_eps) << endl;
    // cone_coeff(eps_x,eps_y,E_eps);

    eval_rectangle(eps_x,eps_y, E);
    // test(E);

    // LDVector v0{0.9468923401720671061132517L,
                    // -4.072102120831082499146823e-24L,
                    // 0.05316795353478707980175375L,
                    // -5.553112274845604899656097e-08L,
                    // -0.01115319054270743243833971L,
                    // -9.199674025000000205149813e-08L};

                    
    // long double z = 0.05316795030019990465013731L;
    // LDVector lap{0.94690970780356629,0,z,0,-0.011166782657715382,0};

    // CR3BP<long double> vf;
    // LDMatrix D(6,6), D1(6,6), D2(6,6);
    // LDVector v0 = vf.findVerticalLyapunovOrbit(lap);
    // LDVector v1 = vf.pm_y(v0,D1);
    // LDVector v2 = vf.pm_y(v1,D2);

    // D1 = vf.pm_y.computeDP(v1,D1);
    // D2 = vf.pm_y.computeDP(v2,D2);
    // D = D2 * D1;

    // LDVector w0{0, 0, 0, 0, 0, 0.0001};
    // LDMatrix T = energy_change_of_basis(v0,vf.E);
    // // cout << v0 + T * w0 << endl;
    // cout << matrixAlgorithms::gaussInverseMatrix(T) * D * T << endl;
}