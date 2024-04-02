#include <iostream>

#include "Parameters.H"

int main(int argc, char* argv[]) {
    
    parameters parms = read_input_files(argc, argv);
    
    std::cout << "p: " << parms.p << std::endl;

    return 0;
}









