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

//* Returns the out-of-plane angle between atoms a, b, c and d in rad
double Molecule::oop(int a, int b, int c, int d){

    double ebcd_x = unit(1,c,b) * unit(2,d,b) - unit(2,c,b) * unit(1,d,b);
    double ebcd_y = unit(2,c,b) * unit(0,d,b) - unit(0,c,b) * unit(2,d,b);
    double ebcd_z = unit(0,c,b) * unit(1,d,b) - unit(1,c,b) * unit(0,d,b);

    double exx = ebcd_x * unit(0,c,a);
    double eyy = ebcd_y * unit(1,c,a);
    double ezz = ebcd_z * unit(2,c,a);

    double theta;
    double sin_theta = (exx + eyy + ezz) / sin(angle(b,c,d));

    if (sin_theta < -1.0) theta = asin(-1.0);
    else if (sin_theta > 1.0) theta = asin(1.0);
    else theta = asin(sin_theta);

    return theta;
}

//* Returns the torsion angle between atoms a, b, c and d in rad
double Molecule::torsion(int a, int b, int c, int d){
    double eabc_x = unit(1, b, a) * unit(2, a, c) - unit(2, b, a) * unit(1, a, c);
    double eabc_y = unit(2, b, a) * unit(0, a, c) - unit(0, b, a) * unit(2, a, c);
    double eabc_z = unit(0, b, a) * unit(1, a, c) - unit(1, b, a) * unit(0, a, c);

    double ebcd_x = unit(1, a, c) * unit(2, c, d) - unit(2, a, c) * unit(1, c, d);
    double ebcd_y = unit(2, a, c) * unit(0, c, d) - unit(0, a, c) * unit(2, c, d);
    double ebcd_z = unit(0, a, c) * unit(1, c, d) - unit(1, a, c) * unit(0, c, d);

    double exx = eabc_x * ebcd_x;
    double eyy = eabc_y * ebcd_y;
    double ezz = eabc_z * ebcd_z;
    double numerator = exx + eyy + ezz;
    double denominator = sin(angle(a,b,c)) * sin(angle(b,c,d));

    double tau = 0.0;
    double cos_tau = numerator / denominator;
    if (cos_tau > 1.0) tau = acos(1.0);
    else if (cos_tau < -1.0) tau = acos(-1.0);
    else tau = acos(cos_tau);

    //* Sign of torsion angle
    double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
    double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
    double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
    double norm = sqrt(pow(cross_x,2) + pow(cross_y,2) + pow(cross_z,2));

    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;

    double sign = 1.0; //* positive by default
    double dot = cross_x * unit(0,b,c) + cross_y * unit(1,b,c) + cross_z * unit(2,b,c);
    if (dot < 0.0) sign = -1.0;

    return tau*sign;
}

//* Global array with the atomic masses
// double mass[] = {
// }

//* Returns the mass of atom a
double Molecule::mass(int a){
    return Z_vals[a];
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