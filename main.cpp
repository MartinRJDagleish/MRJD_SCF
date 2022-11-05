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
    
    // //* Initialize R-Matrix for distance between atoms 
    // double **R = new double*[mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     R[i] = new double[mol.num_atoms];
    // }

    // //* Calculate distance between atoms
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < mol.num_atoms; j++) {
    //         R[i][j] = sqrt( pow((mol.geom[i][0] - mol.geom[j][0]),2) 
    //                         + pow((mol.geom[i][1] - mol.geom[j][1]), 2) 
    //                         + pow((mol.geom[i][2] - mol.geom[j][2]), 2));
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

    // //* Delete allocated mem for R-Matrix 
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     delete R[i];
    // }
    // delete [] R;
    // -------------------------------------------------------------

    Molecule mol("../Project1_Geometries/Acetaldehyd.dat", 0);

    cout << "Number of atoms: " << mol.num_atoms << endl;
    cout << "Input Cartesian coordinates:\n" << endl;
    mol.print_geom();

    cout << "Distance matrix:\n" << endl;
    for (int i = 0; i < mol.num_atoms; i++) {
        for (int j = 0; j < i; j++) {
            printf("%d %d %20.12f \n", i, j, mol.bond(i,j));
        }
    }

    //* Bond Angles MRJD try
    //* Initialize e_vecs
    double *e_vec_1 = new double[3];
    double *e_vec_2 = new double[3];

    for (int i = 0; i < num_atoms; i++) {
        for (int j = 0; j < i; j++){
            // TODO: evtl. Matrix erstellen und da die Vektoren speichern? -> Komponenten berechen 
        }
    }


    delete[] e_vec_1;
    delete[] e_vec_2;


    return 0;
}
