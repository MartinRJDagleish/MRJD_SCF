#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio> 

// new for math
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    // -------------------------------------------------------------
    //! INITIAL VERSION OF CODE
    // ifstream input("geom.dat");

    // int num_atoms;
    // input >> num_atoms;

    // int *Z_val = new int[num_atoms];
    // double *x = new double[num_atoms];
    // double *y = new double[num_atoms];
    // double *z = new double[num_atoms];

    // for(int i = 0; i < num_atoms; i++) {
    //     input >> Z_val[i] >> x[i] >> y[i] >> z[i]; 
    // }
    // input.close();

    // // For testing and learning 
    // cout << "Number of atoms: " << num_atoms << endl;
    // cout << "Input Cartesian coordinates:\n" << endl;

    // for(int i = 0; i < num_atoms; i++) {
    //     printf("%d %20.12f %20.12f %20.12f\n", (int) Z_val[i], x[i], y[i], z[i]);
    // }
    // //! WICHTIG: of course outside of fstream block
    // delete [] Z_val; 
    // delete [] x; delete [] y; delete [] z;
    // std::cout << "Hello my, world!\n";

    //* VERSION 2 for local R-Matrix + bonds + angles 

    //* Initialize R-Matrix for distance between atoms 
    // double **R = new double*[mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     R[i] = new double[mol.num_atoms];
    // }

    // //* Calculate distance between atoms
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < i; j++) {
    //         R[i][j] = R[j][i] = sqrt( pow((mol.geom[i][0] - mol.geom[j][0]),2) 
    //                                 + pow((mol.geom[i][1] - mol.geom[j][1]), 2) 
    //                                 + pow((mol.geom[i][2] - mol.geom[j][2]), 2));
    //     }
    // }

    // //* Print R-Matrix
    // cout << "Distance matrix:\n" << endl;
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < mol.num_atoms; j++) {
    //         printf("%20.12f", R[i][j]);
    //     }
    //     cout << endl;
    // }

    // cout << "Distance matrix:\n" << endl;
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < i; j++) {
    //         printf("%d %d %20.12f \n", i, j, R[i][j]);
    //     }
    // }

    // //* Bond Angles 
    // //* Initialize e_vecs
    // double **e_x = new double* [mol.num_atoms];
    // double **e_y = new double* [mol.num_atoms];
    // double **e_z = new double* [mol.num_atoms];

    // for (int i = 0; i < mol.num_atoms; i++){
    //     e_x[i] = new double[mol.num_atoms];
    //     e_y[i] = new double[mol.num_atoms];
    //     e_z[i] = new double[mol.num_atoms];
    // }

    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < i; j++){
    //         e_x[i][j] = e_x[j][i] = - (mol.geom[i][0] - mol.geom[j][0]) / R[i][j];
    //         e_y[i][j] = e_y[j][i] = - (mol.geom[i][1] - mol.geom[j][1]) / R[i][j];
    //         e_z[i][j] = e_z[j][i] = - (mol.geom[i][2] - mol.geom[j][2]) / R[i][j];
    //     }
    // }


    // //* Print e_x, e_y, e_z Matrix
    // // for (int i = 0; i < mol.num_atoms; i++){
    // //     for (int j = 0; j < i; j++){
    // //         printf("%10.10f %10.10f %10.10f\n", e_x[i][j], e_y[i][j], e_z[i][j]);
    // //     }
    // // }

    // double ***phi = new double** [mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++){
    //     phi[i] = new double* [mol.num_atoms];
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         phi[i][j] = new double[mol.num_atoms];
    //     }
    // }


    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < i; j++){
    //         for (int k = 0; k < j; k++){
    //             if (R[i][j] < 4.0 && R[j][k] < 4.0){
    //                 phi[i][j][k] = acos(e_x[j][i] * e_x[j][k] 
    //                         + e_y[j][i] * e_y[j][k] 
    //                         + e_z[j][i] * e_z[j][k]);
    //             }
    //         }
    //     }
    // }


    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < i; j++){
    //         for (int k = 0; k < j; k++){
    //             printf("%d-%d-%d %20.10f\n", i, j, k, phi[i][j][k]);
    //         }
    //     }
    // }



    // //* Deallocation of ALL 
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++)
    //         delete[] phi[i][j];
    //     delete[] phi[i];
    // }
    // delete[] phi;
    // // phi = nullptr;

    // //* Delete allocated mem for R-Matrix 
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     delete R[i];
    // }
    // delete [] R;

    // for(int i = 0; i < mol.num_atoms; i++){
    //     delete e_x[i]; delete e_y[i]; delete e_z[i];
    // }
    // delete[] e_x; delete[] e_y; delete[] e_z; 

    
    // -------------------------------------------------------------

    Molecule mol("../Project1_Geometries/Acetaldehyd.dat", 0);

    cout << "Number of atoms: " << mol.num_atoms << endl;
    cout << "Input Cartesian coordinates:\n" << endl;
    mol.print_geom();

    cout << "\nInteratomic distances (bohr units):" << endl;
    for (int i = 0; i < mol.num_atoms; i++){
        for (int j = 0; j < i; j++){
            printf("%d %d %8.5f \n", i, j, mol.bond(i, j));
        }
    }

    cout << "\nBond angles (degree):" << endl;
    for (int i = 0; i < mol.num_atoms; i++){
        for (int j = 0; j < i; j++){
            for (int k = 0; k < j; k++){
                if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0){
                    printf("%d %d %d %8.5f \n", i, j, k, 180 * mol.angle(i, j, k) / M_PI);
                }
            }
        }
    }

    cout << "\nDihedral angles (degree):" << endl;



    return 0;
}
