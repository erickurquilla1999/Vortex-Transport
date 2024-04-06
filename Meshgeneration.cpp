#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>


#include "Meshgeneration.H"
#include "Utilities.H"

mesh generate_mesh(const parameters& parms){

    std::vector<double> grids_cord_x(parms.num_element_in_x+1);
    std::vector<double> grids_cord_y(parms.num_element_in_y+1);

    double element_size_x = parms.domain_x / parms.num_element_in_x;
    double element_size_y = parms.domain_y / parms.num_element_in_y;

    // generate x and y coordinates of the grid on the grid boundary
    for (int i = 0; i < parms.num_element_in_x + 1; ++i) {
        grids_cord_x[i] = -1 * parms.domain_x / 2 + element_size_x * i ;
    }
    for (int i = 0; i < parms.num_element_in_y + 1; ++i) {
        grids_cord_y[i] = -1 * parms.domain_y / 2 + element_size_y * i ;
    }

    std::vector<double> allgridpoints_x( ( parms.num_element_in_x + 1 ) * ( parms.num_element_in_y + 1 ) );
    std::vector<double> allgridpoints_y( ( parms.num_element_in_x + 1 ) * ( parms.num_element_in_y + 1 ) ); 

    // generate all the grid points cordinates
    int counter = 0;
    for (int i = 0; i < parms.num_element_in_y+1; ++i){
        for (int j = 0; j < parms.num_element_in_x+1; ++j){
            allgridpoints_x[counter] = grids_cord_x[j] + (i == 0 ? 0 : 1) * (i == parms.num_element_in_y ? 0 : 1) * (j == 0 ? 0 : 1) * (j == parms.num_element_in_x ? 0 : 1) * parms.perturb_grid * 0.2 * element_size_x * sin ( pow( grids_cord_x[j] + 6 , grids_cord_y[i] + 6 ) ) ;
            allgridpoints_y[counter] = grids_cord_y[i] + (i == 0 ? 0 : 1) * (i == parms.num_element_in_y ? 0 : 1) * (j == 0 ? 0 : 1) * (j == parms.num_element_in_x ? 0 : 1) * parms.perturb_grid * 0.2 * element_size_y * cos ( pow( grids_cord_y[i] + 6 , grids_cord_x[j] + 6 ) ) ;
            counter++;
        }
    }

    std::cout << "\nGrid coodinates information ( x , y )" << std::endl;
    for (int i = 0; i < ( parms.num_element_in_x+1 ) * ( parms.num_element_in_y+1 ) ; ++i) {
        std::cout << i <<" : ( " << allgridpoints_x[i] << " , " << allgridpoints_y[i] << " )"<< std::endl;
    }

    std::vector<int> el_to_nod_1( 2 * parms.num_element_in_x * parms.num_element_in_y ); // node one is the square angle
    std::vector<int> el_to_nod_2( 2 * parms.num_element_in_x * parms.num_element_in_y ); // node two is the next to the node one counter clockwise
    std::vector<int> el_to_nod_3( 2 * parms.num_element_in_x * parms.num_element_in_y ); // node two is the next to the node two counter clockwise

    // generate element to node information 
    counter = 0;
    // node one is the square angle
    for (int i = 1; i < parms.num_element_in_y+1; ++i){
        for (int j = 0; j < parms.num_element_in_x; ++j){
            el_to_nod_1[counter] =  j + ( parms.num_element_in_x + 1 ) * ( i - 1) ;
            counter++; 
            el_to_nod_1[counter] =  ( j + ( parms.num_element_in_x + 1 ) * i ) + 1 ;
            counter++; 
        }
    }
    // node two is the next to the node one counter clockwise
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        el_to_nod_2[i] = el_to_nod_1[i] + (i % 2 == 0 ? 1 : -1);
    }
    // node two is the next to the node two counter clockwise
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        el_to_nod_3[i] = el_to_nod_2[(i % 2 == 0 ? i+1 : i-1)];
    }

    std::cout << "\nElement to node information\nThe coordinates of the squere angle vertex is the first number (see below for the coodinates in the physical space)\nThe other indices represent the other vertices going counter clockwise" << std::endl;
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        std::cout << i << " : " << el_to_nod_1[i] << " , " << el_to_nod_2[i] << " , " << el_to_nod_3[i] << std::endl;
    }

    // element type 0: square angle down, 1: squeare angle up    
    std::vector<int> element_type(2 * parms.num_element_in_x * parms.num_element_in_y);
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        element_type[i] = (i % 2 == 0 ? 0 : 1);
    }

    std::cout << "\nElement type information\n0 represents an element with the square angle down\n1 represents an element with the squere angle up" << std::endl;
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        std::cout << i << " : " << element_type[i] << std::endl;
    }

    // element boundary information
    std::vector<int> element_right   ( 2 * parms.num_element_in_x * parms.num_element_in_y);
    std::vector<int> element_left    ( 2 * parms.num_element_in_x * parms.num_element_in_y);
    std::vector<int> element_vertical( 2 * parms.num_element_in_x * parms.num_element_in_y);
    
    // element to the right
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i){
        element_right[i] = (i % ( 2 * parms.num_element_in_x ) == 0 ? i - 1 + 2 * parms.num_element_in_x : i - 1 );
    }    

    // element to the left
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i){
        element_left[i] = ( ( i + 1 ) % ( 2 * parms.num_element_in_x ) == 0 ? i + 1 - 2 * parms.num_element_in_x : i + 1 );
    }

    // vertical element: up for even element number and down for odd element number
    for (int i = 0; i < 2 * parms.num_element_in_x; ++i){
        element_vertical[i] = ( i % 2 == 0 ?  i + 1 + 2 * parms.num_element_in_x * ( parms.num_element_in_y - 1 ) : i - 1 + 2 * parms.num_element_in_x );
    }
    for (int i = 2 * parms.num_element_in_x ; i < 2 * parms.num_element_in_x * ( parms.num_element_in_y - 1 ) ; ++i){
        element_vertical[i] = ( i % 2 == 0 ?  i + 1 - 2 * parms.num_element_in_x : i - 1 + 2 * parms.num_element_in_x );
    }
    for (int i = 2 * parms.num_element_in_x * ( parms.num_element_in_y - 1 ); i < 2 * parms.num_element_in_x * parms.num_element_in_y ; ++i){
        element_vertical[i] = ( i % 2 == 0 ?  i + 1 - 2 * parms.num_element_in_x : i - 1 - 2 * parms.num_element_in_x * ( parms.num_element_in_y - 1 ) );
    }

    std::cout << "\nElement boundaries\nFirst element in the element on the right\nSecond number is the element on the left\nThird number is the element in the vertical direction, down for even numbers and up for odd number" << std::endl;
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i){
        std::cout << i << " : " << element_right[i] << " , " << element_left[i] << " , " << element_vertical[i] << std::endl;
    }


    // if you have time improve the above line
    // ..................................................................


    // ..................................................................
    // preparing data to return
    // ..................................................................

    // create the array that contains the element number
    std::vector<int> elem_n( 2 * parms.num_element_in_x * parms.num_element_in_y );
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y ; ++i) {
        elem_n[i] = i;
    }

    // create the array that contains the elements vertices coordinates
    std::vector<std::vector<std::vector<double>>> ele_coords(2 * parms.num_element_in_x * parms.num_element_in_y, std::vector<std::vector<double>>(3, std::vector<double>(2)));
    for (int i = 0; i <  2 * parms.num_element_in_x * parms.num_element_in_y  ; ++i) {
        ele_coords[i][0][0] = allgridpoints_x[el_to_nod_1[i]];
        ele_coords[i][0][1] = allgridpoints_y[el_to_nod_1[i]];
        ele_coords[i][1][0] = allgridpoints_x[el_to_nod_2[i]];
        ele_coords[i][1][1] = allgridpoints_y[el_to_nod_2[i]];
        ele_coords[i][2][0] = allgridpoints_x[el_to_nod_3[i]];
        ele_coords[i][2][1] = allgridpoints_y[el_to_nod_3[i]];
    }

    // create the array that contains the elements at the boundaries of each element
    std::vector<std::vector<int>> ele_at_bondry( 2 * parms.num_element_in_x * parms.num_element_in_y , std::vector<int>(3) );
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y  ; ++i) {
        ele_at_bondry[i][0] = element_right[i];
        ele_at_bondry[i][1] = element_left[i];
        ele_at_bondry[i][2] = element_vertical[i];
    }

    // create the mesh data structure to be returned
    mesh grid;
    grid.element_number = elem_n;              // Element number
    grid.element_coordinates = ele_coords;     // First index represent element number. Second index represent vertices of the element (run between 0 and 2): 0 is the vertice at the square angle, 1 and 2 are the other vertices going counter clock wise. Third index represent x:0 and y:1 position
    grid.elements_at_boundary = ele_at_bondry; // Cointains the elements number of the elements at the boundaries. The first index in the element number of the current element. The second index runs over 0 and 2. The index 0 is the element at the right, 1 at the left and 2 in the vertical direction.
    grid.element_type = element_type;          // Contains zero or uno. Zero is for elements with square anglue down. One for squera anglue up.

    // ..................................................................
    // writing mesh data
    // ..................................................................

    // Create grid directory to store grid information
    std::string dirPath = "grid";
    std::string command = "rm -rf " + dirPath; // Remove directory and its contents
    int status = system(command.c_str());

    if (status == 0) {
        std::cout << "Directory cleaned successfully: " << dirPath << std::endl;
    } else {
        std::cerr << "Failed to clean directory: " << dirPath << std::endl;
        exit(EXIT_FAILURE);
    }

    // Now create the grid directory
    command = "mkdir -p " + dirPath;
    status = system(command.c_str());

    if (status == 0) {
        std::cout << "Directory created successfully: " << dirPath << std::endl;
    } else {
        std::cerr << "Failed to create directory: " << dirPath << std::endl;
        exit(EXIT_FAILURE);
    }

    // saving grid information for each element in grid directory
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {

        std::vector<std::string> lines(11);
        lines[0]="element_number=" + std::to_string(i);
        lines[1]="coordinate_1_x=" + std::to_string(ele_coords[i][0][0]);
        lines[2]="coordinate_1_y=" + std::to_string(ele_coords[i][0][1]);
        lines[3]="coordinate_2_x=" + std::to_string(ele_coords[i][1][0]);
        lines[4]="coordinate_2_y=" + std::to_string(ele_coords[i][1][1]);
        lines[5]="coordinate_3_x=" + std::to_string(ele_coords[i][2][0]);
        lines[6]="coordinate_3_y=" + std::to_string(ele_coords[i][2][1]);
        lines[7]="type=" + std::to_string(element_type[i]);
        lines[8]="right_element=" + std::to_string(ele_at_bondry[i][0]);
        lines[9]="left_element=" + std::to_string(ele_at_bondry[i][1]);
        lines[10]="vertical_element=" + std::to_string(ele_at_bondry[i][2]);

        writeToFile("grid/element" + std::to_string( i ) + ".txt", lines);
    }

    // ..................................................................
    // return mesh data
    // ..................................................................

    return grid;
}