#include <string>

using namespace std;

class Molecule {
    public:
        int num_atoms;
        int charge;
        int *Z_vals; // array of atomic numbers
        double **geom; // array of atomic coordinates
        string point_group; // point group of molecule

        void print_geom();
        void rotate(double phi);
        void translate(double x, double y, double z);
        double bond(int atom1, int atom2);
        double angle(int atom1, int atom2, int atom3);
        double torsion(int atom1, int atom2, int atom3, int atom4);
        double unit(int cart, int atom1, int atom2);

        Molecule(const char *filename, int q); //* default constructor
        ~Molecule(); //* destructor of class Molecule 
};