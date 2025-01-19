#include <map>
#include <string>
#include "helper.cpp"

class Model {
    /*
    Model class of HMM, used for transmission values storage.
    */
    public:
        std::map <std::string, double> pi;
        std::map <std::string, double> A;
        std::map <std::string, double> E;

        Model(std::string pi_path, std::string A_path, std::string E_path){
            this->pi = read_pi(pi_path);
            this->A = read_A(A_path);
            this->E = read_E(E_path);
        }
};