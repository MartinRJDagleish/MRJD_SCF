#include "molecule.h"

using namespace std;

int main(int argc, char *argv[]){
    Molecule h2o;

    h2o.num_atoms = 3;
    h2o.charge = 0;
    h2o.Z_vals = new int[h2o.num_atoms]; // defines array of atomic numbers
    
}