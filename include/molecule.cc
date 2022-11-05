#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>

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