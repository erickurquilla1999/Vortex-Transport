#include <iostream>
#include <cmath>
#include <cstdlib>

#include "Meshgeneration.H"

void generate_mesh(const parameters& parms){

    float grids_cord_x[parms.num_element_in_x+1];
    float grids_cord_y[parms.num_element_in_y+1];

    float element_size_x = parms.domain_x / parms.num_element_in_x;
    float element_size_y = parms.domain_y / parms.num_element_in_y;

    // generate x and y coordinates of the grid on the grid boundary
    for (int i = 0; i < parms.num_element_in_x + 1; ++i) {
        grids_cord_x[i] = -1 * parms.domain_x / 2 + element_size_x * i ;
    }
    for (int i = 0; i < parms.num_element_in_y + 1; ++i) {
        grids_cord_y[i] = -1 * parms.domain_y / 2 + element_size_y * i ;
    }

    float allgridpoints_x[ ( parms.num_element_in_x + 1 ) * ( parms.num_element_in_y + 1 ) ];
    float allgridpoints_y[ ( parms.num_element_in_x + 1 ) * ( parms.num_element_in_y + 1 ) ]; 

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

    float el_to_nod_1[ 2 * parms.num_element_in_x * parms.num_element_in_y ]; // node one is the square angle
    float el_to_nod_2[ 2 * parms.num_element_in_x * parms.num_element_in_y ]; // node two is the next to the node one counter clockwise
    float el_to_nod_3[ 2 * parms.num_element_in_x * parms.num_element_in_y ]; // node two is the next to the node two counter clockwise

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
    float element_type[2 * parms.num_element_in_x * parms.num_element_in_y];
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        element_type[i] = (i % 2 == 0 ? 0 : 1);
    }

    std::cout << "\nElement type information\n0 represents an element with the square angle down\n1 represents an element with the squere angle up" << std::endl;
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y; ++i) {
        std::cout << i << " : " << element_type[i] << std::endl;
    }

    // element boundary information
    float element_right[2 * parms.num_element_in_x * parms.num_element_in_y];
    float element_left[2 * parms.num_element_in_x * parms.num_element_in_y];
    float element_vertical[2 * parms.num_element_in_x * parms.num_element_in_y];
    
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

    // create grid directory to store grid information
    std::string dirPath = "grid";
    std::string command = "mkdir -p " + dirPath;
    int status = system(command.c_str());

    if (status == 0) {
        std::cout << "grid directory created successfully: " << dirPath << std::endl;
    } else {
        std::cerr << "Failed to create grid directory: " << dirPath << std::endl;
        exit(EXIT_FAILURE);
    }
}