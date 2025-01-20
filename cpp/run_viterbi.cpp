#include "viterbi.cpp"

/*Function to run viterbi algorithm with chosen parameters*/
int main(void){
    Model model = Model("../transmission_values/trained_from_estimate/pi.txt",
    "../transmission_values/trained_from_estimate/A.txt",
    "../transmission_values/trained_from_estimate/E.txt");
    std::string folderPath = "../data/data_test";
    int i = 0;
    for (const auto &entry : std::filesystem::directory_iterator(folderPath)){
        std::ifstream file(entry.path());
        std::string seq1;
        std::string seq2;
        getline(file, seq1);
        getline(file, seq2);
        std::pair<std::string, std::string> pair = std::make_pair(seq1, seq2);
        std::pair<std::string, std::string> aligned_pair = viterbi_algo(model, pair);

        std::string filePathWrite = "../data/alignments_hmm_estim_train/";
        filePathWrite += "aligned_no_" + std::to_string(i); 
        std::ofstream outFile(filePathWrite);
        if (outFile.is_open()) {
            outFile << aligned_pair.first << '\n';
            outFile << aligned_pair.second;
            outFile.close();
        } else {
            std::cerr << "Failed to open file for writing.\n";
        }
        std::cout << "Completed pair "<< i << std::endl;
        i++;
    }

}