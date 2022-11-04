#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio> 

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
    
    // -------------------------------------------------------------
    Molecule mol("geom.dat", 0);

    cout << "Number of atoms: " << mol.num_atoms << endl;
    cout << "Input Cartesian coordinates:\n" << endl;
    mol.print_geom();
    return 0;
}
