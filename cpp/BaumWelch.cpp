#include <map>
#include <string>
#include "helper.cpp"
#include <vector>
#include <math.h>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <fstream>
const double INF = -1e10;
int MAXITER = 10;

std::vector<std::string> Pi_keys = {"MM", "Ix", "Iy"};
std::vector<std::string> A_keys = {"MM -> MM", "MM -> Ix", "MM -> Iy", "Ix -> MM", "Ix -> Ix", "Ix -> Iy", "Iy -> MM", "Iy -> Ix", "Iy -> Iy"};

std::vector<std::string> E_keys = {"MM -> AA", "MM -> AC", "MM -> AT", "MM -> AG", "MM -> CA", "MM -> CC", "MM -> CT", "MM -> CG",
                         "MM -> TA", "MM -> TC", "MM -> TT", "MM -> TG", "MM -> GA", "MM -> GC", "MM -> GT", "MM -> GG",
                         "Ix -> A-", "Ix -> C-", "Ix -> T-", "Ix -> G-",
                         "Iy -> -A", "Iy -> -C", "Iy -> -T", "Iy -> -G"};



// Baum-Welch algorithm implementation
void baumWelch(const std::vector<int>& observations, std::vector<double>& Pi, std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& E) {
    /*
    Function for baum welch algorithm for sequence observation

    Arguments:
        - const std::vector<int>& observations: sequence of observation pairs turned to integer
        - std::vector<double>& Pi: vector of pi probabilities
        - std::vector<std::vector<double>>& A: vector of A probabilities
        - std::vector<std::vector<double>>& E: vector of E probabilities
    */
    size_t N = Pi.size(); // Number of states
    size_t T = observations.size(); // Length of observation sequence

    // Forward and backward variables
    std::vector<std::vector<double>> alpha(T, std::vector<double>(N));
    std::vector<std::vector<double>> beta(T, std::vector<double>(N));

    // Step 1: Forward procedure
    for (size_t i = 0; i < N; ++i) {
        alpha[0][i] = Pi[i] + E[i][observations[0]];    
                        
    }

    for (size_t t = 1; t < T; ++t) {
        for (size_t i = 0; i < N; ++i) {
            std::vector<double> logContributions;
            for (size_t j = 0; j < N; ++j) {
                logContributions.push_back(alpha[t - 1][j] + A[j][i] + E[i][observations[t]]);
            }
            alpha[t][i] = logSumExp(logContributions);
        }
    }

    // Step 2: Backward procedure
    for (size_t i = 0; i < N; ++i) {
        beta[T - 1][i] = 0;
    }

    for (int t = T - 2; t >= 0; --t) {
        for (size_t i = 0; i < N; ++i) {
            std::vector<double> logContributions;
            for (size_t j = 0; j < N; ++j) {
                logContributions.push_back(A[i][j] + E[j][observations[t + 1]] + beta[t + 1][j]);
            }
            beta[t][i] = logSumExp(logContributions);
        }
    }

    // Step 3: Re-estimation
    std::vector<std::vector<double>> gamma(T, std::vector<double>(N));
    std::vector<std::vector<std::vector<double>>> xi(T - 1, std::vector<std::vector<double>>(N, std::vector<double>(N)));

    for (size_t t = 0; t < T - 1; ++t) {
        std::vector<double> logContributions;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                xi[t][i][j] = alpha[t][i] + A[i][j] + E[j][observations[t + 1]] + beta[t + 1][j];
                logContributions.push_back(xi[t][i][j]);
                    
            }
        }
        double logSum = logSumExp(logContributions);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {

                xi[t][i][j] -= logSum; 
                    
            }
        }
    }


    for (size_t t = 0; t < T; ++t) {
        std::vector<double> logContributions;
        for (size_t i = 0; i < N; ++i) {
            gamma[t][i] = alpha[t][i] + beta[t][i];
            logContributions.push_back(gamma[t][i]);
        }
        double logSum = logSumExp(logContributions);
        for (size_t i = 0; i < N; ++i) {
  
            gamma[t][i] -= logSum; 
        }
    }

    for (size_t i = 0; i < N; ++i) {
            
        Pi[i] = gamma[0][i];
            
    }

    // Update A
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            std::vector<double> numeratorlogContributions;
            std::vector<double> denominatorlogContributions;
            for (size_t t = 0; t < T - 1; ++t) {
                numeratorlogContributions.push_back(xi[t][i][j]);
                denominatorlogContributions.push_back(gamma[t][i]);
            }
            double numerator = logSumExp(numeratorlogContributions);
            double denominator = logSumExp(denominatorlogContributions);
            if (std::isfinite(denominator)) {
                A[i][j] = numerator - denominator; // Convert back to linear space
            } else {
                A[i][j] = INF;  //really large number to symbolize -inf
            }
                
        }
    }
    // Update E
    for (size_t i = 0; i < N; ++i) {
        for (size_t k = 0; k < E[i].size(); ++k) {
            std::vector<double> numeratorlogContributions;
            std::vector<double> denominatorlogContributions;
            for (size_t t = 0; t < T; ++t) {
                if (isSymbolEqualToObservation(k, observations[t])) {
                    numeratorlogContributions.push_back(gamma[t][i]);
                }
                denominatorlogContributions.push_back(gamma[t][i]);
            }
            double numerator = logSumExp(numeratorlogContributions);

            double denominator = logSumExp(denominatorlogContributions);
            if (std::isfinite(denominator)) {
                E[i][k] = numerator - denominator; // Convert back to linear space
            } else {
                E[i][k] = INF; // Handle case where denominator is zero (logDenominator = -inf)
            } 
                
        }
    }
}


int main() {
    std::vector<double> Pi;
    std::vector<std::vector<double>> A, E;
    std::vector<std::string> testNuclPairs = {"AA", "AT", "TT", "A-", "G-", "-C", "GG"}; //used for testing
    readInitialParameters(Pi, A, E);
    // inserting zeros into E (emmmission probabilities) so that E is a 3x24 matrix
    // Add 8 zeros to the back for Match (MM) state
    E[0].insert(E[0].end(), 8, 0.0);
    // Add 16 zeros to the front and 4 zeros to the back for Nucleotide and Gap (Ix) state
    E[1].insert(E[1].begin(), 16, 0.0);
    E[1].insert(E[1].end(), 4, 0.0);
    // Add 20 zeros to the front for Gap and Nucleotide (Iy) state
    E[2].insert(E[2].begin(), 20, 0.0);
    std::vector<std::string> possibleObservations = {"AA", "AC", "AT", "AG", "CA", "CC", "CT", "CG",
                                                "TA", "TC", "TT", "TG", "GA", "GC", "GT", "GG",
                                                "A-", "C-", "T-", "G-",
                                                "-A", "-C", "-T", "-G"};

    // Iterative expectation-maximization
    for (int iter = 0; iter < MAXITER; ++iter) { 
        for (const auto &entry : std::filesystem::directory_iterator("../data/train_data/baum_welch_train/")){
            std::string filePath_Pair = entry.path().string();
            std::vector<std::string> observationsNuclPairs = readAndPairSequences(filePath_Pair);
             // PRETVARANJE OBSERVACIJA U integere, odnosno indexe
            std::vector<int> observationsInts;
            // Use a for loop to append each element to the list
            for (int i = 0; i < observationsNuclPairs.size(); ++i) {    //testNuclPairs   //observationsNuclPairs
                int observacija = transformObservationToInt(observationsNuclPairs[i], possibleObservations);  //testNuclPairs  //observationsNuclPairs
                observationsInts.push_back(observacija);
            }
             // RUN BAUM WELCH
            baumWelch(observationsInts, Pi, A, E);
            normalize_E(E);
        }
    }
    
    // Writing updated parameters to files
    std::string filePath_Pi = "../transmission_values/trained_from_estimate/pi.txt";
    std::string filePath_A = "../transmission_values/trained_from_estimate/A.txt";
    std::string filePath_E = "../transmission_values/trained_from_estimate/E.txt";

    writePiToFile(Pi, filePath_Pi);
    write_A_ToFile(A, filePath_A);
    write_E_ToFile(E, filePath_E);
    

    return 0;
}