#include "molecule.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>

#include <cmath>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//* Eigen library
#include "Eigen/Dense" 
#include "Eigen/Eigenvalues" 
#include "Eigen/Core" 
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

//* NISTConst
#define NISTCONST_COMMON_SYMBOLS_NAMES // Common symbols and names for constants.
#include "include/NISTConst.hpp"

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

    //* VERSION 2 for local R-Matrix + bonds + bond-angles 

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

    
    // //* VERSION 3 for local R-Matrix + bonds + bond-angles + out-of-plane angles

    // //* Initialize R-Matrix for distance between atoms 
    // double **R = new double*[mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     R[i] = new double[mol.num_atoms];
    // }

    // // //* Calculate distance between atoms
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     for (int j = 0; j < i; j++) {
    //         R[i][j] = R[j][i] = sqrt( pow((mol.geom[i][0] - mol.geom[j][0]),2) 
    //                                 + pow((mol.geom[i][1] - mol.geom[j][1]), 2) 
    //                                 + pow((mol.geom[i][2] - mol.geom[j][2]), 2));
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

    // //* Out-of-plane angles
    // double ****theta_oop = new double*** [mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++){
    //     theta_oop[i] = new double** [mol.num_atoms];
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         theta_oop[i][j] = new double*[mol.num_atoms];
    //         for(int k = 0; k < mol.num_atoms; k++){
    //             theta_oop[i][j][k] = new double[mol.num_atoms];
    //         }
    //     }
    // }


    // //* with this implementation the phi matrix is not needed 
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             for (int l = 0; l < mol.num_atoms; l++){
    //                 theta_oop[i][j][k][l] = (e_x[k][i] * (e_y[k][j] * e_z[k][l] - e_z[k][j] * e_y[k][l]) 
    //                                          + e_y[k][i] * (e_z[k][j] * e_x[k][l] - e_x[k][j] * e_z[k][l]) 
    //                                          + e_z[k][i] * (e_x[k][j] * e_y[k][l] - e_y[k][j] * e_x[k][l])) 
    //                                             / sin(mol.angle(j,k,l)); //* technically, this is sin(theta_oop), but we need the theta check 
    //                 if (theta_oop[i][j][k][l] < -1.0) theta_oop[i][j][k][l] = 180*asin(-1.0)/M_PI;
    //                 else if(theta_oop[i][j][k][l] > 1.0) theta_oop[i][j][k][l] = 180*asin(1.0)/M_PI;
    //                 else theta_oop[i][j][k][l] = 180*asin(theta_oop[i][j][k][l])/M_PI;
    //                 // printf("%20.10f \n", theta_oop[i][j][k][l]);
    //                 } 
    //             }
    //         }
    //     }

    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             for (int l = 0; l < mol.num_atoms; l++){
    //                 if(i!=j && i!=k && i!=l && j!=k && j!=l && k!=l){
    //                     if(R[i][k] < 4.0 && R[k][j] < 4.0 && R[k][l] < 4.0) {
    //                         printf("%d-%d-%d-%d %f\n", i, j, k, l, theta_oop[i][j][k][l]);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }



    // //* Deallocation of ALL 
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             delete[] theta_oop[i][j][k];
    //         }
    //         delete[] theta_oop[i][j];
    //     }
    //     delete[] theta_oop[i];
    // }
    // delete[] theta_oop;

    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++)
    //         delete[] phi[i][j];
    //     delete[] phi[i];
    // }
    // delete[] phi;

    // //* Delete allocated mem for R-Matrix 
    // for (int i = 0; i < mol.num_atoms; i++) {
    //     delete R[i];
    // }
    // delete [] R;

    // for(int i = 0; i < mol.num_atoms; i++){
    //     delete e_x[i]; delete e_y[i]; delete e_z[i];
    // }
    // delete[] e_x; delete[] e_y; delete[] e_z; 

    // //* VERSION 4 for local R-Matrix + bonds + bond-angles + out-of-plane angles + torsion angles
    // //* Allocation of tau Matrix
    // double ****tau = new double***[mol.num_atoms];
    // for (int i = 0; i < mol.num_atoms; i++){
    //     tau[i] = new double**[mol.num_atoms];
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         tau[i][j] = new double*[mol.num_atoms];
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             tau[i][j][k] = new double[mol.num_atoms];
    //             }
    //         }
    //     }

    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             for (int l = 0; l < mol.num_atoms; l++){
    //                 if (i != j && i != k && i != l && j != k && j != l && k != l){
    //                     if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0){

    //                         double eijk_x = mol.unit(1, i, j) * mol.unit(2, j, k) - mol.unit(2, i, j) * mol.unit(1, j, k);
    //                         double eijk_y = mol.unit(2, i, j) * mol.unit(0, j, k) - mol.unit(0, i, j) * mol.unit(2, j, k);
    //                         double eijk_z = mol.unit(0, i, j) * mol.unit(1, j, k) - mol.unit(1, i, j) * mol.unit(0, j, k);

    //                         double ejkl_x = mol.unit(1, j, k) * mol.unit(2, k, l) - mol.unit(2, j, k) * mol.unit(1, k, l);
    //                         double ejkl_y = mol.unit(2, j, k) * mol.unit(0, k, l) - mol.unit(0, j, k) * mol.unit(2, k, l);
    //                         double ejkl_z = mol.unit(0, j, k) * mol.unit(1, k, l) - mol.unit(1, j, k) * mol.unit(0, k, l);

    //                         double exx = eijk_x * ejkl_x;
    //                         double eyy = eijk_y * ejkl_y;
    //                         double ezz = eijk_z * ejkl_z;
    //                         double numerator = exx + eyy + ezz;
    //                         double denominator = sin(mol.angle(i,j,k)) * sin(mol.angle(j,k,l));

    //                         double cos_tau = numerator / denominator;
    //                         if (cos_tau > 1.0) tau[i][j][k][l] = 180*acos(1.0)/M_PI;
    //                         else if (cos_tau < -1.0) tau[i][j][k][l] = 180*acos(-1.0)/M_PI;
    //                         else tau[i][j][k][l] = 180*acos(cos_tau)/M_PI;
    //                         printf("%d-%d-%d-%d %8.5f \n", i, j, k, l, tau[i][j][k][l]);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }



    // //* Deallocation of ALL
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             delete[] tau[i][j][k];
    //         }
    //         delete[] tau[i][j];
    //     }
    //     delete[] tau[i];
    // }
    // delete[] tau;


    // -------------------------------------------------------------

    // Linux 
    Molecule mol("../Project1_Geometries/Allene.dat", 0);
    // Windows
    // Molecule mol("../../Project1_Geometries/Acetaldehyd.dat", 0);

//* ASCII SOURCE: 
// https://patorjk.com/software/taag/#p=testall&h=1&f=Blocks&t=%20SCF%20
    string header = R"(
.----------------.  .----------------.  .----------------.  .----------------.  
| .--------------. || .--------------. || .--------------. || .--------------. | 
| | ____    ____ | || |  _______     | || |     _____    | || |  ________    | | 
| ||_   \  /   _|| || | |_   __ \    | || |    |_   _|   | || | |_   ___ `.  | |
| |  |   \/   |  | || |   | |__) |   | || |      | |     | || |   | |   `. \ | | 
| |  | |\  /| |  | || |   |  __ /    | || |   _  | |     | || |   | |    | | | |  
| | _| |_\/_| |_ | || |  _| |  \ \_  | || |  | |_' |     | || |  _| |___.' / | | 
| ||_____||_____|| || | |____| |___| | || |  `.___.'     | || | |________.'  | | 
| |              | || |              | || |              | || |              | | 
| '--------------' || '--------------' || '--------------' || '--------------' | 
 '----------------'  '----------------'  '----------------'  '----------------'  )";

    string header2 = R"(

__/\\\\____________/\\\\_____/\\\\\\\\\___________/\\\\\\\\\\\___/\\\\\\\\\\\\____        
 _\/\\\\\\________/\\\\\\___/\\\///////\\\________\/////\\\///___\/\\\////////\\\__       
  _\/\\\//\\\____/\\\//\\\__\/\\\_____\/\\\____________\/\\\______\/\\\______\//\\\_      
   _\/\\\\///\\\/\\\/_\/\\\__\/\\\\\\\\\\\/_____________\/\\\______\/\\\_______\/\\\_     
    _\/\\\__\///\\\/___\/\\\__\/\\\//////\\\_____________\/\\\______\/\\\_______\/\\\_    
     _\/\\\____\///_____\/\\\__\/\\\____\//\\\____________\/\\\______\/\\\_______\/\\\_   
      _\/\\\_____________\/\\\__\/\\\_____\//\\\____/\\\___\/\\\______\/\\\_______/\\\__  
       _\/\\\_____________\/\\\__\/\\\______\//\\\__\//\\\\\\\\\_______\/\\\\\\\\\\\\/___ 
        _\///______________\///___\///________\///____\/////////________\////////////_____

 ________________/\\\\\\\\\\\___________/\\\\\\\\\___/\\\\\\\\\\\\\\\____________         
  ______________/\\\/////////\\\______/\\\////////___\/\\\///////////_____________        
   _____________\//\\\______\///_____/\\\/____________\/\\\________________________       
    ______________\////\\\___________/\\\______________\/\\\\\\\\\\\________________      
     _________________\////\\\_______\/\\\______________\/\\\///////_________________     
      ____________________\////\\\____\//\\\_____________\/\\\________________________    
       _____________/\\\______\//\\\____\///\\\___________\/\\\________________________   
        ____________\///\\\\\\\\\\\/_______\////\\\\\\\\\__\/\\\________________________  
         ______________\///////////____________\/////////___\///_________________________ 

    )";

    string header3 = R"(
███╗   ███╗██████╗      ██╗██████╗ 
████╗ ████║██╔══██╗     ██║██╔══██╗
██╔████╔██║██████╔╝     ██║██║  ██║
██║╚██╔╝██║██╔══██╗██   ██║██║  ██║
██║ ╚═╝ ██║██║  ██║╚█████╔╝██████╔╝
╚═╝     ╚═╝╚═╝  ╚═╝ ╚════╝ ╚═════╝ 
)";                                   


    cout << header2 << "\n" << endl;

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

    cout << "\nOut-of-plane angles (degree):" << endl;
    for (int i = 0; i < mol.num_atoms; i++){
        for (int j = 0; j < mol.num_atoms; j++){
            for (int k = 0; k < mol.num_atoms; k++){
                for (int l = 0; l < mol.num_atoms; l++){
                    if (i != j && i != k && i != l && j != k && j != l && k != l){
                        if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0){
                            printf("%d-%d-%d-%d %8.5f \n", i, j, k, l, 180 * mol.oop(i, j, k, l) / M_PI);
                        }
                    }
                }
            }
        }
    }
    
    cout << "\nDihedral angles (degree):" << endl;
    for (int i = 0; i < mol.num_atoms; i++){
        for (int j = 0; j < i; j++){
            for (int k = 0; k < j; k++){
                for (int l = 0; l < k; l++){
                    if (i != j && i != k && i != l && j != k && j != l && k != l){
                        if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0){
                            printf("%d-%d-%d-%d %8.5f \n", i, j, k, l, 180 * mol.torsion(i, j, k, l) / M_PI);
                        }
                    }
                }
            }
        }
    }

    cout << "\nCenter of mass:" << endl;
    double CM_x = 0.0;
    double CM_y = 0.0;
    double CM_z = 0.0;

    double tot_mass = 0.0;

    for (int i = 0; i < mol.num_atoms; i++){
        CM_x += mol.mass(mol.Z_vals[i]) * mol.geom[i][0];
        CM_y += mol.mass(mol.Z_vals[i]) * mol.geom[i][1];
        CM_z += mol.mass(mol.Z_vals[i]) * mol.geom[i][2];
        tot_mass += mol.mass(mol.Z_vals[i]);
    }
    CM_x /= tot_mass;
    CM_y /= tot_mass;
    CM_z /= tot_mass;

    printf("  Total mass: %10.5f\n", tot_mass);
    printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", CM_x, CM_y, CM_z);

    mol.translate(-CM_x, -CM_y, -CM_z);
    cout << "\nGeometry after translation:" << endl;
    mol.print_geom();

    cout << "\nMoment of inertia tensor (u * bohr**2):" << endl;

    //* ---------------------------------------------- 
    //* Standard C++ without Eigen 
    //* Allocation of tensor and stuff
    // double **inertia_tensor = new double*[3];
    // for (int i = 0; i < 3; i++){
    //     inertia_tensor[i] = new double[3];
    // }

    // for (int i = 0; i < 3; i++){
    //     for (int j = 0; j < 3; j++){
    //         if (i == j){
    //             for (int k = 0; k < mol.num_atoms; k++){
    //                 inertia_tensor[i][i] += mol.mass(mol.Z_vals[k]) * (pow(mol.geom[k][1], 2) + pow(mol.geom[k][2], 2));
    //             }
    //         }
    //         else{
    //             for (int k = 0; k < mol.num_atoms; k++){
    //                 inertia_tensor[i][j] -= mol.mass(mol.Z_vals[k]) * mol.geom[k][i] * mol.geom[k][j];
    //             }
    //         }
    //     }
    // }

    // for (int i = 0; i < 3; i++){
    //     for (int j = 0; j < 3; j++){
    //         printf("%12.8f ", inertia_tensor[i][j]);
    //     }
    //     printf("\n");
    // }

    // //* Deallocation of tensor and stuff
    // for (int i = 0; i < 3; i++){
    //     delete[] inertia_tensor[i];
    // }
    // delete[] inertia_tensor;
    //* ---------------------------------------------- 
    //* C++ without Eigen 

    Matrix Inertia_T(3,3);

    //! WRONG!!!!
    // for (int i = 0; i < 3; i++){
    //     for (int j = 0; j < 3; j++){
    //         if (i == j){
    //             if (i == j ==1){
    //                 for (int k = 0; k < mol.num_atoms; k++){
    //                     Inertia_T(i,j) = Inertia_T(j,i) += mol.mass(mol.Z_vals[k]) 
    //                             * (pow(mol.geom[k][1], 2) + pow(mol.geom[k][2], 2));
    //                     }
    //             }
    //         }
    //         else{
    //             for (int k = 0; k < mol.num_atoms; k++){
    //                 Inertia_T(i,j) = Inertia_T(j,i) -= mol.mass(mol.Z_vals[k]) 
    //                      * mol.geom[k][i] * mol.geom[k][j];
    //             }
    //         }
    //     }
    // }
    for(int i = 0; i < mol.num_atoms; i++){
        Inertia_T(0,0) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][1], 2) + pow(mol.geom[i][2], 2));
        Inertia_T(1,1) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][0], 2) + pow(mol.geom[i][2], 2));
        Inertia_T(2,2) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][0], 2) + pow(mol.geom[i][1], 2));
        Inertia_T(0,1) = Inertia_T(1,0) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][0] * mol.geom[i][1];
        Inertia_T(0,2) = Inertia_T(2,0) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][0] * mol.geom[i][2];
        Inertia_T(1,2) = Inertia_T(2,1) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][1] * mol.geom[i][2];
    }

    cout << "\nInertia tensor:" << endl;
    cout << Inertia_T << endl;

    //* Working eigenvals and eigenvectors
    Eigen::SelfAdjointEigenSolver<Matrix> solver(Inertia_T);
    Matrix eigenvecs = solver.eigenvectors();
    Matrix eigenvals = solver.eigenvalues();

    //* Debug printing of eigenvals and eigenvectors
    // cout << "\nEigenvalues:" << endl;
    // cout << eigenvals << endl;
    // cout << "\nEigenvectors:" << endl;
    // cout << eigenvecs << endl;

    cout << "\nPrincipal moments of inertia (u * bohr**2):" << endl;
    cout << eigenvals << endl;

    cout << "\nPrincipal moments of inertia (u * Angstrom**2):" << endl;
    cout << eigenvals 
        * 1/pow(NISTConst::AngstromStar * pow(10,10),2) 
        * pow(NISTConst::BohrRadius * pow(10,10),2) << endl;
    // printf("Test again: %12.8f\n ", NISTConst::AngstromStar * pow(10,10));

    cout << "\nPrincipal moments of inertia (g * cm**2):" << endl;
    cout << eigenvals * NISTConst::atomicMassConstant 
         * pow(10,3) * pow(NISTConst::BohrRadius, 2) 
         * pow(100,2) << endl; 
    
    //* Classification of the rotor type
    if(mol.num_atoms == 2) cout << "\nMolecule is diatomic.\n";
    else if(eigenvals(0) < 1e-4) cout << "\nMolecule is linear.\n";
    else if((fabs(eigenvals(0) - eigenvals(1)) < 1e-4) && 
            (fabs(eigenvals(1) - eigenvals(2)) < 1e-4)) {
        cout << "\nMolecule is a spherical top.\n";
        }
    else if((fabs(eigenvals(0) - eigenvals(1)) < 1e-4) && 
            (fabs(eigenvals(1) - eigenvals(2)) > 1e-4)){
        cout << "\nMolecule is an oblate symmetric top.\n";
    }
    else if((fabs(eigenvals(0) - eigenvals(1)) > 1e-4) && 
            (fabs(eigenvals(1) - eigenvals(2)) < 1e-4)){
        cout << "\nMolecule is a prolate symmetric top.\n";
    }
    else cout << "\nMolecule is an asymmetric top.\n";

    //* Calculation of the rotational constants
    cout << "\nRotational constants (cm**-1):" << endl;
    Vector rot_consts_cm(3);
    Vector eigenvals_cm(3);

    for (int i = 0; i < 3; i++){
        //* unit conversions first
        eigenvals_cm(i) = eigenvals(i) * NISTConst::atomicMassConstant 
            * pow(NISTConst::BohrRadius, 2) * pow(10,4); 
        rot_consts_cm(i) = NISTConst::h * pow(10,2) / 
            (8 * pow(M_PI,2) * NISTConst::c * eigenvals_cm(i));
    }
    cout << "A = " << rot_consts_cm(0) << "  B = " << rot_consts_cm(1) <<
        "  C = " << rot_consts_cm(2) << endl;

    cout << "\nRotational constants (MHz):" << endl;
    Vector rot_consts_MHz(3);
    Vector eigenvals_MHz(3);

    for (int i = 0; i < 3; i++){
        eigenvals_MHz(i) = eigenvals(i) * NISTConst::atomicMassConstant 
            * pow(NISTConst::BohrRadius, 2); 
        rot_consts_MHz(i) = NISTConst::h * pow(10,-6) / 
            (8 * pow(M_PI,2) * eigenvals_MHz(i));

    }
    cout << "A = " << rot_consts_MHz(0) << "  B = " << rot_consts_MHz(1) <<
        "  C = " << rot_consts_MHz(2) << endl;

    //* ---------------------------------------------- 
    //* Test of NISTConst
    // const double test_const = NISTConst::h;
    // cout << test_const << endl;

    return 0;
}
