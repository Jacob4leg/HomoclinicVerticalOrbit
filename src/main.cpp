#include<iostream>
#include "capd/capdlib.h"
#include "cr3bp.h"

using namespace std;
using namespace capd;

int main() {
    cout.precision(20);

    // LDVector v0{0.9468923401720671061132517L,
                    // -4.072102120831082499146823e-24L,
                    // 0.05316795353478707980175375L,
                    // -5.553112274845604899656097e-08L,
                    // -0.01115319054270743243833971L,
                    // -9.199674025000000205149813e-08L};

    
    long double z = 0.05316795030019990465013731L;
    LDVector lap{0.94690970780356629,0,z,0,-0.011166782657715382,0};

    CR3BP<long double> vf;
    auto v0 = lap;
    auto v1 = vf.pm_y(v0);
    auto v2 = vf.pm_y(v1);

    cout << vf.E(v0) << endl;
    cout << vf.E(v1) << endl;
    cout << vf.E(v2) << endl;
    

}