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
????????????   ?????????????????????????????????      ?????????????????????????????? 
??????????????? ???????????????????????????????????????     ?????????????????????????????????
?????????????????????????????????????????????????????????     ??????????????????  ?????????
???????????????????????????????????????????????????????????????   ??????????????????  ?????????
????????? ????????? ??????????????????  ?????????????????????????????????????????????????????????
?????????     ??????????????????  ????????? ?????????????????? ????????????????????? 
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
    //* Project 3 starts here
    //* Step 1 & 2


    //* -------------------------------------------------------------
    return 0;
}
