#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>

void Molecule::print_geom(){
    for(int i = 0; i < num_atoms; i++) {
        printf("%d %8.5f %8.5f %8.5f\n", Z_vals[i], geom[i][0], geom[i][1], geom[i][2]);
    }
}

void Molecule::translate(double x, double y, double z){
    for(int i = 0; i < num_atoms; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

double Molecule::bond(int a, int b){
    return sqrt( pow((geom[a][0] - geom[b][0]),2) 
                + pow((geom[a][1] - geom[b][1]), 2) 
                + pow((geom[a][2] - geom[b][2]), 2));
}

//* Returns the val of unit vector between atoms a and b
//* in the cartesian dir (cart=0=x, cart=1=y, …)
double Molecule::unit(int cart, int a, int b){
    return - (geom[a][cart] - geom[b][cart]) / bond(a,b);
}

//* Returns the angle between atoms a, b and c in rad 
double Molecule::angle(int a, int b, int c){
    return acos(unit(0,b,a) * unit(0,b,c) + unit(1,b,a) * unit(1,b,c) + unit(2,b,a) * unit(2,b,c));
}

//* Constructor with fstream as input  
Molecule::Molecule(const char *filename, int q){
    charge = q;

    //* open filename
    std::ifstream ifs(filename);
    assert(ifs.good());

    //* read the number of atoms from filename
    ifs >> num_atoms;

    //* allocate space for Z_vals and geom
    Z_vals = new int[num_atoms];
    geom = new double*[num_atoms];
    for(int i = 0; i < num_atoms; i++){
        geom[i] = new double[3];
    }

    //* read the data from file
    for(int i = 0; i < num_atoms; i++){
        ifs >> Z_vals[i];
        ifs >> geom[i][0];
        ifs >> geom[i][1];
        ifs >> geom[i][2];
    }

    //* close filestream -> Python: f.close() 
    ifs.close();
}

// //* Constructor with 2 ints as arguments
// Molecule::Molecule(int n, int q){
//     num_atoms = n;
//     charge = q;
//     Z_vals = new int[num_atoms];
//     geom = new double*[num_atoms];
//     for(int i = 0; i < num_atoms; i++){
//         geom[i] = new double[3];
//     }
// }

Molecule::~Molecule(){
    delete[] Z_vals;
    for(int i = 0; i < num_atoms; i++){
        delete[] geom[i];
    }
    delete[] geom;
}