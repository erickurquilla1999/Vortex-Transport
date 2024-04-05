#include <iostream>
#include <vector>

#include "Parameters.H"
#include "Utilities.H"
#include "Meshgeneration.H"

int main(int argc, char* argv[]) {
    
    // read simulation paramaters
    parameters parms = read_input_files(argc, argv);
    
    // generate mesh
    generate_mesh(parms);

    

    return 0;
}