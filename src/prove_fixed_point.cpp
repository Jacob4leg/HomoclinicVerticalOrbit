#include "prove_fixed_point.h"

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
    
    
    return L;
}


bool prove_fixed_point(long double eps, Interval E) {
    
    auto T = IMatrix(get_change_of_basis());
    
    Interval eps_interval(-eps,eps);
    IVector u{  eps_interval,
                eps_interval,
                eps_interval,
                eps_interval, E};
    
    u = T * u;
    IVector v{u[0],0,u[1],u[2],u[3],u[4]};
    v = IVector(v0) + IMatrix(T_E) * v;
    v[3] = 0.;
    v[5] = 0.;


    IMatrix D(6,6);
    C2Rect2Set S(v);
    IVector y = Ivf.pm_z(S,D);
    D = Ivf.pm_z.computeDP(y,D);
    
    
    IVector w{y[1],y[3]};
    IMatrix M({{D[1][0],D[1][4]},{D[3][0],D[3][4]}});
    
    IVector v1{v0[0],0,v[2],0,v0[4],0};
    C2Rect2Set S1(v1);
    y = Ivf.pm_z(S1,D);
    D = Ivf.pm_z.computeDP(y,D);
    IVector w0{y[1],y[3]};

    IVector N = - matrixAlgorithms::gauss(M,w0);
    cout << v - IVector(v0) << endl;
    cout << N << endl;

    return 1;
}