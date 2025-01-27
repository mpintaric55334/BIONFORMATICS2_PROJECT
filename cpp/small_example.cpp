#include "viterbi.cpp"

/*Running small example*/
int main(void){
    Model model = Model("../transmission_values/trained_from_estimate/pi.txt",
    "../transmission_values/trained_from_estimate/A.txt",
    "../transmission_values/trained_from_estimate/E.txt");
    std::string folderPath = "../data/data_test";
    std::string seq1 = "ACAGGGACTTGAAAGCGAAAGTGAGACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCGC";
    std::string seq2 = "ATGGGTGCGAGAGCGTCGGTATTAAGCGGGGGACAATTAGATAGATGGGAAAAAATTC";
    std::pair<std::string, std::string> pair = std::make_pair(seq1, seq2);
    std::pair<std::string, std::string> aligned_pair = viterbi_algo(model, pair);
    std::cout << "Alignment" << std::endl << aligned_pair.first << std::endl << aligned_pair.second << std::endl;


}