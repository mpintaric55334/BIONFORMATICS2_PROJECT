#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "math.h"
#include <vector>

std::vector<std::string> states = {"MM", "Ix", "Iy"};
std::vector<std::string> nucleotides = {"A", "C", "T", "G"};
std::string gap = "-";

std::map<std::string, double> read_pi(std::string pi_file_path){

    /*
    Function that initializes the pi(starting) transmission probability .

    Arguments:
        - pi_file_path: std::string => file path to the
        pi matrix storage file

    Returns:
        - pi: std::map<std::string, double> pi probability map
    */
    std::map<std::string, double> pi;
    std::ifstream input_file(pi_file_path);
    std::vector<double> probabilities;
    std::string line;
    while (std::getline(input_file, line)) {
        double probability = std::stod(line);
        probabilities.push_back(log(probability));
       
    }
    input_file.close();
    int n = probabilities.size();
    for (int i=0; i<n; i++){
        pi.insert(std::make_pair(states.at(i), probabilities.at(i)));
    }
    return pi;
}


std::map<std::string, double> read_A(std::string A_file_path){
    /*
    Function that initializes the A transmission probability map.

    Arguments:
        - A_file_path: std::string => file path to the
        A probabilities storage file

    Returns:
        - A: std::map<std::string, double> E probability map
    */
    std::map<std::string, double> A;
    std::ifstream input_file(A_file_path);
    std::vector<double> probabilities;
    std::vector<std::string> possible_values;
    for(const auto& state: states){
        for(const auto& state2: states){
            possible_values.push_back(state + state2);
        }
    }
    std::string line;
    while (std::getline(input_file, line)) {
        double probability = std::stod(line);
        probabilities.push_back(log(probability));
       
    }
    input_file.close();
    int n = probabilities.size();
    for (int i=0; i<n; i++){
        A.insert(std::make_pair(possible_values.at(i), probabilities.at(i)));
    }
    return A;
}


std::map<std::string, double> read_E(std::string E_file_path){
    /*
    Function that initializes the E transmission probability map.

    Arguments:
        - E_file_path: std::string => file path to the
        E probabilities storage file

    Returns:
        - E: std::map<std::string, double> E probability map
    */
    std::map<std::string, double> E;
    std::ifstream input_file(E_file_path);
    std::vector<double> probabilities;
    std::vector<std::string> possible_values;
    for (const auto& nucl: nucleotides){
        for (const auto& nucl2: nucleotides){
            possible_values.push_back(states.at(0) + nucl + nucl2);
        }
    }
    // state Ix
    for (const auto& nucl: nucleotides){
        possible_values.push_back(states.at(1) + nucl + gap);
    }
    // state Iy
    for (const auto& nucl: nucleotides){
        possible_values.push_back(states.at(2) + gap + nucl);
    }
    std::string line;
    while (std::getline(input_file, line)) {
        double probability = std::stod(line);
        probabilities.push_back(log(probability));
       
    }
    input_file.close();
    int n = probabilities.size();
    for (int i=0; i<n; i++){
        E.insert(std::make_pair(possible_values.at(i), probabilities.at(i)));
    }
    return E;
}
