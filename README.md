This is a project done for Bionformatics 2 course at University of Zagreb, Faculty of Electrical Engineering and Computing.

The goal of the project is to make a Hidden Markov Model for pairwise sequence alignment.

To get started, either clone the repository or unzip the project files.

All the data is ready, parameters are optimized, and aligned sequences are stored, but if you wish to modify the process yourself, follow these steps.

1 => DATA PROCESSING:
    - run python file processing.py
        - data directories, and data size are already set, but if you want to add custom ones, change them in load_data function calling
    - run python file estimate.py
        - statistical estimates of HMM parameters are already written into files, and paths are set, but if you want to do custom ones, you can change it in this file
    - estimations are by default stored in transmission_values/estimate. Also, by hand, we created a random transmission prob files, stored in random subdirectory. Feel free to consult help_desctiption.txt file for description of prbability transmission files.

2 => PARAMETER OPTIMIZATION:
    - optimized parameters are already stored in transmission_values/trained_from_estimate. However, if you wish to optimize parameters of some other estimation, or random parameter initialization, inside helper.cpp, in readInitialParameters function, change the location of initial parameters, and in BaumWelch.cpp, change the data directory of output transmission which will best suit you. Default is initial parameters: transmission_values/estimate and trained parameters: transmission_values/trained_from_estimate. This combination of initial parameters and Baum Welch algorithm has shown to yield best alignments.
    - once you have chosen your parameter directories, compile Baum-Welch cpp using g++ -o baumwelch BaumWelch.cpp, then just run the exe file.

3 => VITERBI ALGORITHM:
    - to run the viterbi algorithm, go to run_viterbi.cpp, if you wish, change what parameter files you want to use in Model initialization, and if you want, change output directory of your alignments. Default is, parameters: transmission_values/trained_from_estimate, and output directory data/alignments_hmm_estim_train. You will notice several other alignments folder, containing alignments with other possible combination of parameters, such as alignments_hmm_estim, which stores alignments using only estimate parameters, or alignments_hmm_random, storing alignments generated wiht random model parameters, or alignments_hmm_rand_train, using alignments generated from optimized random parameters.
    - to run the algorithm, compile the file with g++ -o run_viterbi run_viterbi.cpp, and run the exe file

4 => NEEDLEMAN_WUNSCH
    - if you wish to compare alignments with needleman wunsch alignments, you can generate these alignments using needle_wunch.py. The file uses minineedle python library to align sequences. It stores the alignments in data/alignments_needle

5 => INFERENCE AND SCORING
    - to generate mean scores and score plots, run python file scoring.py

6 => SMALL EXAMPLE RUNNING
    - g++ -o small_example small_example.cpp, run exe file