#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "math.h"
#include <vector>
#include <algorithm>

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
double logSumExp(const std::vector<double>& logValues) {
    if(logValues.size()==0){
        return -1e9;
    }
    double maxLog = *std::max_element(logValues.begin(), logValues.end());
    double sumExp = 0.0;
    for (double logValue : logValues) {
        sumExp += std::exp(logValue - maxLog);
    }
    return maxLog + std::log(sumExp);
}
bool isSymbolEqualToObservation(int symbol, int observation) {
    return symbol == observation;
}

int transformObservationToInt(const std::string& A, const std::vector<std::string>& list) {
    // Iterate through the list and find the index of the given string A
    for (size_t i = 0; i < list.size(); ++i) {
        if (list[i] == A) {
            return static_cast<int>(i); // Return the index if found
        }
    }
    return -1; // Return -1 if the string A is not found in the list
}

// Function to read initial parameters from files
void readInitialParameters(std::vector<double>& Pi, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& E) {
    std::ifstream piFile("../transmission_values/estimate/pi.txt");
    std::ifstream aFile("../transmission_values/estimate/A.txt");
    std::ifstream eFile("../transmission_values/estimate/E.txt");

    if (!piFile || !aFile || !eFile) {
        std::cerr << "Error reading parameter files." << std::endl;
        exit(1);
    }

    // Read Pi
    double value;
    while (piFile >> value) {
        Pi.push_back(log(value));
    }

    // Read A
    while (aFile) {
        std::vector<double> row;
        for (size_t i = 0; i < Pi.size() && aFile >> value; ++i) {
            row.push_back(log(value));
        }
        if (!row.empty()) A.push_back(row);
    }

    // Read E
    size_t cnt = 0;
    size_t iterate = 16;
    while (eFile) {
        std::vector<double> row;
        if (cnt == 16) iterate = 4;
        for (size_t i = 0; i < iterate && eFile >> value; ++i) {
            row.push_back(log(value));
            cnt++;
        }
        if (!row.empty()) E.push_back(row);
    }
}


// PRINTAČA FUNKCIJA ZA PROVJERU GDJE POĐE PO ZLU
void printVectorOfVectors(const std::vector<std::vector<double>>& vecOfVecs) {
    // Loop through each vector in the vector of vectors
    for (const auto& vec : vecOfVecs) {
        // Loop through each element in the current vector and print it
        for (const auto& element : vec) {
            std::cout << element << " ";
        }
        // Print a new line after each inner vector
        //std::cout << std::endl;
    }
}

void writePiToFile(const std::vector<double>& vec, const std::string& filePath) {
    // Open the file in write mode
    std::ofstream outFile(filePath);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
        return;
    }

    // Write each element of the vector to the file, one per line
    for (double value : vec) {
        outFile << exp(value) << "\n";
    }

    // Close the file
    outFile.close();

    std::cout << "Pi written to " << filePath << " successfully." << std::endl;
}

// Function for writing matrix A to txt file
void write_A_ToFile(const std::vector<std::vector<double>>& vec, const std::string& filePath) {
    // Open the file in write mode
    std::ofstream outFile(filePath);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
        return;
    }

    // Write each element of the matrix to the file, one per line
    for (const auto& row : vec) {
        for (double val : row) {
            outFile << exp(val) << "\n";
        }
    }

    // Close the file
    outFile.close();

    std::cout << "A written to " << filePath << " successfully." << std::endl;
}

// Function for writing matrix E to txt file
void write_E_ToFile(const std::vector<std::vector<double>>& vec, const std::string& filePath) {
    // Open the file in write mode
    std::ofstream outFile(filePath);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << " for writing." << std::endl;
        return;
    }

    // Writing MM observations
    std::vector<double> mm_observation_probs = vec[0];
    double sum__mm = 0;
    for (int i = 0; i < 16; ++i) {
        sum__mm += exp(mm_observation_probs[i]);
    }
    for (int i = 0; i < 16; ++i) {
        outFile << exp(mm_observation_probs[i])/sum__mm << "\n";
    }

    // Writing Ix observations
    std::vector<double> Ix_observation_probs = vec[1];
    double sum_x = 0;
    for (int i = 16; i < 20; ++i) {
        sum_x += exp(Ix_observation_probs[i]);
    }

    for (int i = 16; i < 20; ++i) {
        outFile << exp(Ix_observation_probs[i])/sum_x << "\n";
    }

    // Writing Iy observations
    std::vector<double> Iy_observation_probs = vec[2];
    double sum_y = 0;
    for (int i = 20; i < 24; ++i) {
        sum_y += exp(Iy_observation_probs[i]);
    }
    for (int i = 20; i < 24; ++i) {
        outFile << exp(Iy_observation_probs[i])/sum_y << "\n";
    }

    // Close the file
    outFile.close();

    std::cout << "E written to " << filePath << " successfully." << std::endl;
}


// READING 2 sequences from a file and returning the pairs of nucleotides
std::vector<std::string> readAndPairSequences(const std::string& filePath) {
    std::vector<std::string> pairedSequences;

    // Open the file for reading
    std::ifstream inFile(filePath);

    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return pairedSequences; // Return empty vector on error
    }

    // Read the two lines (sequences) from the file
    std::string sequence1, sequence2;
    std::getline(inFile, sequence1);
    std::getline(inFile, sequence2);

    // Check if both sequences are of the same length
    if (sequence1.length() != sequence2.length()) {
        std::cerr << "Error: Sequences have different lengths." << std::endl;
        return pairedSequences; // Return empty vector on error
    }

    // Create pairs of corresponding characters
    for (size_t i = 0; i < sequence1.length(); ++i) {
        std::string pair = "";
        pair += sequence1[i];
        pair += sequence2[i];
        pairedSequences.push_back(pair);
    }

    return pairedSequences;
}


void normalize_E(std::vector<std::vector<double>> &E){
    double sum__mm = 0;
    for (int i = 0; i < 16; ++i) {
        sum__mm += exp(E[0][i]);
    }
    for (int i = 0; i < 16; ++i) {
        E[0][i] = exp(E[0][i])/sum__mm;
    }

    double sum_x = 0;
    for (int i = 16; i < 20; ++i) {
        sum_x += exp(E[1][i]);
    }

    for (int i = 16; i < 20; ++i) {
        E[1][i] = exp(E[1][i])/sum_x;
    }
    double sum_y = 0;
    for (int i = 20; i < 24; ++i) {
        sum_y += exp(E[2][i]);
    }
    for (int i = 20; i < 24; ++i) {
        E[2][i] = exp(E[2][i])/sum_y;
    }

}