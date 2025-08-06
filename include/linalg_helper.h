#include<iostream>
#include "capd/capdlib.h"


#ifndef _LINALG_HELPER_H_
#define _LINALG_HELPER_H_

capd::LDMatrix energy_change_of_basis(capd::LDVector p, capd::LDMap &E) {

    capd::LDMatrix DE = E.derivative(p);
    capd::LDVector E_grad{DE[0][0],DE[0][1], DE[0][3], DE[0][4], DE[0][5],-1};
    capd::LDVector z_grad = -E_grad / DE[0][2];

    capd::LDMatrix T = {{1,0,0,0,0,0},
                        {0,1,0,0,0,0},
                        {z_grad[0],z_grad[1],z_grad[2],z_grad[3],z_grad[4],z_grad[5]},
                        {0,0,1,0,0,0},
                        {0,0,0,1,0,0},
                        {0,0,0,0,1,0}};
    return T;
}

capd::LDMatrix us_change_of_basis(capd::LDMatrix A) {
    std::pair<int,int> dim = A.dimension();
    int n = std::get<0>(dim);
    int m = std::get<1>(dim);
    assert(n == m);

    capd::LDVector rV(n), iV(n);
    capd::LDMatrix rVec(n,n), iVec(n,n);
    capd::alglib::computeEigenvaluesAndEigenvectors(A,rV,iV,rVec,iVec);

    for(int i = 0; i < n; i++) {
        if(iV[i] != 0) {
            rVec.column(i+1) = iVec.column(i);
            i++;
        }
    }
    return rVec;
}


#endif