#ifndef PARAMETERS_H_
#define PARAMETERS_H_

struct parameters{
    int p;
    float domain_x;
    float domain_y;
};

parameters read_input_files(int argc, char* argv[]);

#endif // PARAMETERS_H_
