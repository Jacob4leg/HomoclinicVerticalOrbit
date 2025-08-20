#include<iostream>
#include<future>
#include<numeric>
#include<thread>
#include "capd/capdlib.h"
#include "prove_fixed_point.h"

using namespace std;
using namespace capd;

int main() {
    
    cout.precision(15);

    // long double eps = 1e-6;
    // long double E_eps = 1e-8;

    // long double eps_x = 5e-10;
    // long double eps_y = 1e-8;
    // test();
    // return 0;

    long double eps_x = 1e-7;
    // long double eps_y = 2e-9;
    long double eps_y = 1.2e-8;
    // long double E_eps = 1e-12;
    long double E_eps = 1e-12;
    Interval E(-E_eps,E_eps);
    // cout << prove_fixed_point(eps_x,eps_y,E_eps) << endl;
    // cone_coeff(eps_x,eps_y,E_eps);

    test(E);

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