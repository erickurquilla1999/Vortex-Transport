#include <iostream>
#include <vector>

#include "Parameters.H"
#include "Utilities.H"
#include "Meshgeneration.H"

int main(int argc, char* argv[]) {
    
    // read simulation paramaters
    parameters parms = read_input_files(argc, argv);
    
    // generate mesh
    mesh simulation_mesh = generate_mesh(parms);

    std::vector<std::string> lines(11);
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        lines[0]="\nelement_number=" + std::to_string(simulation_mesh.element_number[i]);
        lines[1]="coordinate_1_x=" + std::to_string(simulation_mesh.element_coordinates[i][0][0]);
        lines[2]="coordinate_1_y=" + std::to_string(simulation_mesh.element_coordinates[i][0][1]);
        lines[3]="coordinate_2_x=" + std::to_string(simulation_mesh.element_coordinates[i][1][0]);
        lines[4]="coordinate_2_y=" + std::to_string(simulation_mesh.element_coordinates[i][1][1]);
        lines[5]="coordinate_3_x=" + std::to_string(simulation_mesh.element_coordinates[i][2][0]);
        lines[6]="coordinate_3_y=" + std::to_string(simulation_mesh.element_coordinates[i][2][1]);
        lines[7]="type=" + std::to_string(simulation_mesh.element_type[i]);
        lines[8]="right_element=" + std::to_string(simulation_mesh.elements_at_boundary[i][0]);
        lines[9]="left_element=" + std::to_string(simulation_mesh.elements_at_boundary[i][1]);
        lines[10]="vertical_element=" + std::to_string(simulation_mesh.elements_at_boundary[i][2]);
        for (int j = 0; j < 11; ++j) {
            std::cout << lines[j] << std::endl;
        }
    }

    return 0;
}