#ifndef PARAMETERS_H_
#define PARAMETERS_H_

// input parameters struct
struct parameters{
    int p;
    float domain_x;
    float domain_y;
    int num_element_in_x;
    int num_element_in_y;
    float perturb_grid;

};

// this function reads the input parameter script and return an struct of the type parameters with all of them
parameters read_input_files(int argc, char* argv[]);

#endif // PARAMETERS_H_
