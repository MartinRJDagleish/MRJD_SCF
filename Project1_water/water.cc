#include "molecule.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){

    // Molecule h2o(3,0);
    //! Not needed with new version of code
    // h2o.num_atoms = 3;
    // h2o.charge = 0;
    // h2o.Z_vals = new int[h2o.num_atoms]; // defines array of atomic numbers
    // h2o.geom = new double*[h2o.num_atoms]; // defines array of atomic coordinates

    // for (int i = 0; i < h2o.num_atoms; i++){
    //     h2o.geom[i] = new double[3];
    // }

    // h2o.Z_vals[0] = 8;
    // h2o.geom[0][0] =  0.000000000000;
    // h2o.geom[0][1] =  0.000000000000;
    // h2o.geom[0][2] = -0.122368916506;
    // h2o.Z_vals[1] = 1;
    // h2o.geom[1][0] =  0.000000000000;
    // h2o.geom[1][1] =  1.414995841403;
    // h2o.geom[1][2] =  0.971041753535;
    // h2o.Z_vals[2] = 1;
    // h2o.geom[2][0] =  0.000000000000;
    // h2o.geom[2][1] = -1.414995841403;
    // h2o.geom[2][2] =  0.971041753535;

    // h2o.print_geom();
    // cout << "\n" << endl;
    // h2o.translate(1.0, 1.0, 1.0);
    // h2o.print_geom();

    //! Not needed with new version of code
    // delete [] h2o.Z_vals;
    // for (int i = 0; i < h2o.num_atoms; i++){
    //     delete [] h2o.geom[i];
    // }
    // delete [] h2o.geom;

    Molecule h2o("water.dat", 0);

    h2o.print_geom();
    cout << "\n" << endl;
    h2o.translate(5, 2, 3);
    h2o.print_geom();


    return 0;
}