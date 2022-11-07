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
    // Molecule mol("../Project1_Geometries/water.dat", 0);
    // Windows
    Molecule mol("../../Project1_Geometries/Acetaldehyd.dat", 0);

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

    // cout << "\nInteratomic distances (bohr units):" << endl;
    // for (int i = 0; i < mol.num_atoms; i++)
    // {
    //     for (int j = 0; j < i; j++)
    //     {
    //         printf("%d %d %8.5f \n", i, j, mol.bond(i, j));
    //     }
    // }

    // cout << "\nBond angles (degree):" << endl;
    // for (int i = 0; i < mol.num_atoms; i++)
    // {
    //     for (int j = 0; j < i; j++)
    //     {
    //         for (int k = 0; k < j; k++)
    //         {
    //             if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0)
    //             {
    //                 printf("%d %d %d %8.5f \n", i, j, k, 180 * mol.angle(i, j, k) / M_PI);
    //             }
    //         }
    //     }
    // }

    //* -------------------------------------------------------------
    //* Project 3 starts here
    //* Step 1 and Step 2

    //*****************************
    //* E_nuc_rep read from file
    //*****************************
    // Linux
    const char *E_nuc_rep_filename = "../Project3_files/h2o/STO-3G/enuc.dat";
    // Windows
    // const char *E_nuc_rep_filename = "../../Project3_files/h2o/STO-3G/enuc.dat";
    double E_nuc_rep = 0.0;
    std::ifstream E_nuc_rep_fs(E_nuc_rep_filename);
    assert(E_nuc_rep_fs.good());
    E_nuc_rep_fs >> E_nuc_rep;
    E_nuc_rep_fs.close();

    //*****************************
    //* Overlap matrix S read from file
    //*****************************
    // Linux
    const char *overlap_mat_filename = "../Project3_files/h2o/STO-3G/s.dat";
    // Windows
    // const char *overlap_mat_filename = "../../Project3_files/h2o/STO-3G/s.dat";
    int num_tot_orbitals = 7; // TODO: implementation of number of orbitals
    Matrix S_overlap_mat(num_tot_orbitals, num_tot_orbitals);

    std::ifstream overlap_mat_fs(overlap_mat_filename);
    assert(overlap_mat_fs.good());
    for (int i = 0; i < num_tot_orbitals; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int rubbish;
            double buffer; // * only works like this
            overlap_mat_fs >> rubbish;
            overlap_mat_fs >> rubbish;
            overlap_mat_fs >> buffer;
            S_overlap_mat(i, j) = S_overlap_mat(j, i) = buffer;
        }
    }
    overlap_mat_fs.close();

    cout << "Overlap matrix S: \n"
         << S_overlap_mat << "\n"
         << endl;

    //*****************************
    //* Kinetic energy matrix read from file
    //*****************************
    // Linux
    const char *E_kin_filename = "../Project3_files/h2o/STO-3G/t.dat";
    // Windows
    // const char *E_kin_filename = "../../Project3_files/h2o/STO-3G/t.dat";
    Matrix E_kin_mat(num_tot_orbitals, num_tot_orbitals);

    std::ifstream E_kin_fs(E_kin_filename);
    assert(E_kin_fs.good());
    for (int i = 0; i < num_tot_orbitals; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int rubbish;
            double buffer; // * only works like this
            E_kin_fs >> rubbish;
            E_kin_fs >> rubbish;
            E_kin_fs >> buffer;
            E_kin_mat(i, j) = E_kin_mat(j, i) = buffer;
        }
    }
    E_kin_fs.close();
    cout << "Kinetic energy matrix: \n"
         << E_kin_mat << "\n"
         << endl;

    //*****************************
    //* Nuclear attraction matrix read from file
    //*****************************
    // Linux
    const char *Nuc_Att_filename = "../Project3_files/h2o/STO-3G/v.dat";
    // Windows
    // const char *Nuc_Att_filename = "../../Project3_files/h2o/STO-3G/v.dat";
    Matrix Nuc_Att_mat(num_tot_orbitals, num_tot_orbitals);

    std::ifstream Nut_Att_fs(Nuc_Att_filename);
    assert(Nut_Att_fs.good());
    for (int i = 0; i < num_tot_orbitals; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int rubbish;
            double buffer; // * only works like this
            Nut_Att_fs >> rubbish;
            Nut_Att_fs >> rubbish;
            Nut_Att_fs >> buffer;
            Nuc_Att_mat(i, j) = Nuc_Att_mat(j, i) = buffer;
        }
    }
    Nut_Att_fs.close();
    cout << "Nuclear Attraction matrix: \n"
         << Nuc_Att_mat << "\n"
         << endl;

    //*****************************
    //* Core Hamilonian Matrix
    //*****************************
    Matrix Hamiliton_Mat(num_tot_orbitals, num_tot_orbitals);
    for (int mu = 0; mu < num_tot_orbitals; mu++)
    {
        for (int nu = 0; nu < num_tot_orbitals; nu++)
        {
            Hamiliton_Mat(mu, nu) = E_kin_mat(mu, nu) + Nuc_Att_mat(mu, nu);
        }
    }
    cout << "Core Hamilitonian matrix: \n"
         << Hamiliton_Mat << "\n"
         << endl;

    // * Step 3
    //* Store two-electron repulsion integrals in vec:
    Vector Elec_rep_ints(228); // TODO: implement way of figuring out how many ints there are
    int cmp_idx = 0;
    int mu, nu, lambda, sigma, munu, lambsig, munulambsig;

    if (mu > nu) munu = mu * (mu+1)/2 + nu;
    else munu = nu*(nu+1)/2 + mu;

    if (lambda > sigma) lambsig = lambda * (lambda+1)/2 + sigma;
    else lambsig = sigma*(sigma+1)/2 + lambda;

    if (munu > lambsig) munulambsig = munu * (munu+1)/2 + lambsig;
    else munulambsig = lambsig*(lambsig+1)/2 + munu;


    // for (int mu = 0; mu < num_tot_orbitals; mu++)
    // {
    //     for (int nu = 0; nu <= mu; nu++)
    //     {
    //         for (int lambda = 0; lambda <= num_tot_orbitals; lambda++)
    //         {
    //             for (int sigma = 0; sigma <= lambda; sigma++)
    //             {

    //                 // lambsig = lambda * (lambda + 1) / 2 + sigma;
    //                 // cmp_idx = munu * (munu + 1) / 2 + lambsig;
    //                 // cout << cmp_idx << endl;
    //             }
    //         }
    //     }
    // }

    //* -------------------------------------------------------------
    return 0;
}
