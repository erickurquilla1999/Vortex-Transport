#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include "Parameters.H"

// this function reads the input parameter script and return an struct of the type parameters with all of them
parameters read_input_files(int argc, char* argv[]) {
    // Check if input file is provided as a command-line argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Open input file
    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open input file '" << argv[1] << "'." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Define a map to store parameter names and values
    std::unordered_map<std::string, std::string> paramMap;

    // Read and parse input file
    std::string line;
    while (std::getline(inputFile, line)) {
        // Find the position of '=' in the line
        size_t pos = line.find('=');
        if (pos != std::string::npos) {
            // Extract parameter name and value and store in the map
            paramMap[line.substr(0, pos)] = line.substr(pos + 1);
        }
    }

    // Close input file
    inputFile.close();

    // Create the parameters struct that will be returned
    parameters params;

    // Extract parameters from the paramMap and assign them to the params struct
    params.p = std::stoi(paramMap["p"]);
    params.domain_x = std::stof(paramMap["domain_x"]);
    params.domain_y = std::stof(paramMap["domain_y"]);

    return params;
}