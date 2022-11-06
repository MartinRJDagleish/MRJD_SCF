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

int main(int argc, char *argv[])
{
    // -------------------------------------------------------------
    // Linux
    Molecule mol("../Project1_Geometries/water.dat", 0);
    // Windows
    // Molecule mol("../../Project1_Geometries/Acetaldehyd.dat", 0);

    //* ASCII SOURCE:
    //* https://patorjk.com/software/taag/#p=testall&h=1&f=Blocks&t=%20SCF%20
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

    cout << header2 << "\n"
         << endl;

    cout << "Number of atoms: " << mol.num_atoms << endl;
    cout << "Input Cartesian coordinates:\n"
         << endl;
    mol.print_geom();

    cout << "\nInteratomic distances (bohr units):" << endl;
    for (int i = 0; i < mol.num_atoms; i++)
    {
        for (int j = 0; j < i; j++)
        {
            printf("%d %d %8.5f \n", i, j, mol.bond(i, j));
        }
    }

    cout << "\nBond angles (degree):" << endl;
    for (int i = 0; i < mol.num_atoms; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < j; k++)
            {
                if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
                {
                    printf("%d %d %d %8.5f \n", i, j, k, 180 * mol.angle(i, j, k) / M_PI);
                }
            }
        }
    }

    //* -------------------------------------------------------------
    //* Project 2 starts here
    //* Step 1 & 2

    const char *hess_filename = "../Project2_Hessian/water.hess";

    //* open filename
    std::ifstream hess_ifs(hess_filename);
    assert(hess_ifs.good());

    //* read the number of atoms from filename
    hess_ifs >> mol.num_atoms;

    //* allocate space for Z_vals and geom
    // int hess_rows = mol.num_atoms;
    Matrix hessian(mol.num_atoms * 3, mol.num_atoms * 3);

    //* read the data from file
    for (int i = 0; i < mol.num_atoms * 3; i++)
    {
        for (int j = 0; j < mol.num_atoms * 3; j++)
        {
            hess_ifs >> hessian(i, j);
        }
    }

    //* close filestream -> Python: f.close()
    hess_ifs.close();

    cout << "\nHessian matrix (hartree/bohr^2):" << endl;
    cout << hessian << endl;

    //* Step 3
    for (int i = 0; i < mol.num_atoms; i++)
    {
        for (int j = 0; j < mol.num_atoms; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    hessian(i + k, j + l) = hessian(i + k, j + l) * 1 / sqrt(mol.mass(mol.Z_vals[i]) * mol.mass(mol.Z_vals[j]));
                }
            }
        }
    }

    cout << "\n" << endl;
    cout << hessian << endl;

    //* Working, but long print out (oop angles)
    // cout << "\nOut-of-plane angles (degree):" << endl;
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < mol.num_atoms; j++){
    //         for (int k = 0; k < mol.num_atoms; k++){
    //             for (int l = 0; l < mol.num_atoms; l++){
    //                 if (i != j && i != k && i != l && j != k && j != l && k != l){
    //                     if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0){
    //                         printf("%d-%d-%d-%d %8.5f \n", i, j, k, l, 180 * mol.oop(i, j, k, l) / M_PI);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    //* Working, but long print out (Dihedral angles)
    // cout << "\nDihedral angles (degree):" << endl;
    // for (int i = 0; i < mol.num_atoms; i++){
    //     for (int j = 0; j < i; j++){
    //         for (int k = 0; k < j; k++){
    //             for (int l = 0; l < k; l++){
    //                 if (i != j && i != k && i != l && j != k && j != l && k != l){
    //                     if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0){
    //                         printf("%d-%d-%d-%d %8.5f \n", i, j, k, l, 180 * mol.torsion(i, j, k, l) / M_PI);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    //* Center of mass + Inertia + Rotor
    // cout << "\nCenter of mass:" << endl;
    // double CM_x = 0.0;
    // double CM_y = 0.0;
    // double CM_z = 0.0;
    // double tot_mass = 0.0;
    // for (int i = 0; i < mol.num_atoms; i++)
    // {
    //     CM_x += mol.mass(mol.Z_vals[i]) * mol.geom[i][0];
    //     CM_y += mol.mass(mol.Z_vals[i]) * mol.geom[i][1];
    //     CM_z += mol.mass(mol.Z_vals[i]) * mol.geom[i][2];
    //     tot_mass += mol.mass(mol.Z_vals[i]);
    // }
    // CM_x /= tot_mass;
    // CM_y /= tot_mass;
    // CM_z /= tot_mass;
    //
    // printf("  Total mass: %10.5f\n", tot_mass);
    // printf("\nMolecular center of mass: %12.8f %12.8f %12.8f\n", CM_x, CM_y, CM_z);
    //
    // mol.translate(-CM_x, -CM_y, -CM_z);
    // cout << "\nGeometry after translation:" << endl;
    // mol.print_geom();
    //
    // cout << "\nMoment of inertia tensor (u * bohr**2):" << endl;
    //
    // Matrix Inertia_T(3, 3);
    //
    // //! WRONG!!!!
    // // for (int i = 0; i < 3; i++){
    // //     for (int j = 0; j < 3; j++){
    // //         if (i == j){
    // //             if (i == j ==1){
    // //                 for (int k = 0; k < mol.num_atoms; k++){
    // //                     Inertia_T(i,j) = Inertia_T(j,i) += mol.mass(mol.Z_vals[k])
    // //                             * (pow(mol.geom[k][1], 2) + pow(mol.geom[k][2], 2));
    // //                     }
    // //             }
    // //         }
    // //         else{
    // //             for (int k = 0; k < mol.num_atoms; k++){
    // //                 Inertia_T(i,j) = Inertia_T(j,i) -= mol.mass(mol.Z_vals[k])
    // //                      * mol.geom[k][i] * mol.geom[k][j];
    // //             }
    // //         }
    // //     }
    // // }
    // for (int i = 0; i < mol.num_atoms; i++)
    // {
    //     Inertia_T(0, 0) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][1], 2) + pow(mol.geom[i][2], 2));
    //     Inertia_T(1, 1) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][0], 2) + pow(mol.geom[i][2], 2));
    //     Inertia_T(2, 2) += mol.mass(mol.Z_vals[i]) * (pow(mol.geom[i][0], 2) + pow(mol.geom[i][1], 2));
    //     Inertia_T(0, 1) = Inertia_T(1, 0) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][0] * mol.geom[i][1];
    //     Inertia_T(0, 2) = Inertia_T(2, 0) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][0] * mol.geom[i][2];
    //     Inertia_T(1, 2) = Inertia_T(2, 1) -= mol.mass(mol.Z_vals[i]) * mol.geom[i][1] * mol.geom[i][2];
    // }
    //
    // cout << "\nInertia tensor:" << endl;
    // cout << Inertia_T << endl;
    //
    // //* Working eigenvals and eigenvectors
    // Eigen::SelfAdjointEigenSolver<Matrix> solver(Inertia_T);
    // Matrix eigenvecs = solver.eigenvectors();
    // Matrix eigenvals = solver.eigenvalues();
    //
    // //* Debug printing of eigenvals and eigenvectors
    // // cout << "\nEigenvalues:" << endl;
    // // cout << eigenvals << endl;
    // // cout << "\nEigenvectors:" << endl;
    // // cout << eigenvecs << endl;
    //
    // cout << "\nPrincipal moments of inertia (u * bohr**2):" << endl;
    // cout << eigenvals << endl;
    //
    // cout << "\nPrincipal moments of inertia (u * Angstrom**2):" << endl;
    // cout << eigenvals * 1 / pow(NISTConst::AngstromStar * pow(10, 10), 2) * pow(NISTConst::BohrRadius * pow(10, 10), 2) << endl;
    // // printf("Test again: %12.8f\n ", NISTConst::AngstromStar * pow(10,10));
    //
    // cout << "\nPrincipal moments of inertia (g * cm**2):" << endl;
    // cout << eigenvals * NISTConst::atomicMassConstant * pow(10, 3) * pow(NISTConst::BohrRadius, 2) * pow(100, 2) << endl;
    //
    // //* Classification of the rotor type
    // if (mol.num_atoms == 2)
    //     cout << "\nMolecule is diatomic.\n";
    // else if (eigenvals(0) < 1e-4)
    //     cout << "\nMolecule is linear.\n";
    // else if ((fabs(eigenvals(0) - eigenvals(1)) < 1e-4) &&
    //          (fabs(eigenvals(1) - eigenvals(2)) < 1e-4))
    // {
    //     cout << "\nMolecule is a spherical top.\n";
    // }
    // else if ((fabs(eigenvals(0) - eigenvals(1)) < 1e-4) &&
    //          (fabs(eigenvals(1) - eigenvals(2)) > 1e-4))
    // {
    //     cout << "\nMolecule is an oblate symmetric top.\n";
    // }
    // else if ((fabs(eigenvals(0) - eigenvals(1)) > 1e-4) &&
    //          (fabs(eigenvals(1) - eigenvals(2)) < 1e-4))
    // {
    //     cout << "\nMolecule is a prolate symmetric top.\n";
    // }
    // else
    //     cout << "\nMolecule is an asymmetric top.\n";
    //
    // //* Calculation of the rotational constants
    // cout << "\nRotational constants (cm**-1):" << endl;
    // Vector rot_consts_cm(3);
    // Vector eigenvals_cm(3);
    //
    // for (int i = 0; i < 3; i++)
    // {
    //     //* unit conversions first
    //     eigenvals_cm(i) = eigenvals(i) * NISTConst::atomicMassConstant * pow(NISTConst::BohrRadius, 2) * pow(10, 4);
    //     rot_consts_cm(i) = NISTConst::h * pow(10, 2) /
    //                        (8 * pow(M_PI, 2) * NISTConst::c * eigenvals_cm(i));
    // }
    // cout << "A = " << rot_consts_cm(0) << "  B = " << rot_consts_cm(1) << "  C = " << rot_consts_cm(2) << endl;
    //
    // cout << "\nRotational constants (MHz):" << endl;
    // Vector rot_consts_MHz(3);
    // Vector eigenvals_MHz(3);
    //
    // for (int i = 0; i < 3; i++)
    // {
    //     eigenvals_MHz(i) = eigenvals(i) * NISTConst::atomicMassConstant * pow(NISTConst::BohrRadius, 2);
    //     rot_consts_MHz(i) = NISTConst::h * pow(10, -6) /
    //                         (8 * pow(M_PI, 2) * eigenvals_MHz(i));
    // }
    // cout << "A = " << rot_consts_MHz(0) << "  B = " << rot_consts_MHz(1) << "  C = " << rot_consts_MHz(2) << endl;

    return 0;
}
