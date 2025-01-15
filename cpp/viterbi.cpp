#include "Model.cpp"
#include <vector>
#include <map>
#include <string>
#include <math.h>
#include <algorithm>

const double INF = -1e23;

std::pair<double, int> max(std::vector<double> vec){
    double max_num = INF;
    int max_index = 0;
    for(int i=0; i<vec.size();i++){
        if(vec[i] > max_num){
            max_num = vec[i];
            max_index = i;
        }
    }
    return std::make_pair(max_num, max_index);
}

void viterbi_algo(Model &model, std::pair<std::string, std::string> &pair){
    std::string seq1 = pair.first;
    std::string seq2 = pair.second;
    long n = pair.first.size();
    long m = pair.second.size();

    double delta = model.A["MMIx"] + model.A["MMIy"] - log(2);  /*technically, MMIx and MMIy
    should be same, because if a sequence is above another shouldnt matter, but due
    to small no of sequence (computational limitations), MMIx and MMIy arent the same*/
    double epsilon = model.A["IxIx"] + model.A["IyIy"] - log(2); // same comment as above

    //initial initialization of V vectors, size m+1 and filled with INF
    std::vector<double> Vmm(m+1, INF); 
    std::vector<double> Vx(m+1, INF);
    std::vector<double> Vy(m+1, INF);

    //now, initialize tracker arrays, size (n+1) * (m+1) and filled with -1
    std::vector<std::vector<int>> tracker_mm(n+1, std::vector<int>(m+1, -1));
    std::vector<std::vector<int>> tracker_x(n+1, std::vector<int>(m+1, -1));
    std::vector<std::vector<int>> tracker_y(n+1, std::vector<int>(m+1, -1));
    // initialize 0
    Vmm[0] = 0;

    for (int i=1; i<=n; i++){
        std::vector<double> Vmm_next(m+1, INF); 
        std::vector<double> Vx_next(m+1, INF);
        std::vector<double> Vy_next(m+1, INF);
        for (int j=1; j<=m; j++){

            std::string nucl1 = std::to_string(seq1[i-1]);
            std::string nucl2 = std::to_string(seq2[j-1]);
            //Vmm part of pseudocode
            std::vector<double> Vm_values = {Vmm[j-1], Vx[j-1], Vy[j-1]};
            std::pair<double, int> values_pair = max(Vm_values);
            double max_value = values_pair.first;
            int max_index = values_pair.second;
            std::string search_string = "MM" + nucl1 + nucl2;
            Vmm_next[j] = model.E[search_string] + max_value;
            tracker_mm[i][j] = max_index;
            //Vx part of pseudocode
            std::vector<double> Vx_values = {Vmm[j] + delta, Vx[j] + epsilon};
            values_pair = max(Vx_values);
            max_value = values_pair.first;
            max_index = values_pair.second;
            search_string = "Ix" + nucl1 + "-";
            Vx_next[j] = model.E[search_string] + max_value;
            tracker_x[i][j] = max_index;
            //Vy part of pseudocode
            std::vector<double> Vy_values = {Vmm_next[j-1] + delta, Vy_next[j-1] + epsilon};
            values_pair = max(Vy_values);
            max_value = values_pair.first;
            max_index = values_pair.second;
            search_string = "Iy-" + nucl2;
            Vy_next[j] = model.E[search_string] + max_value;
            tracker_y[i][j] = max_index;
        }
        Vmm = Vmm_next;
        Vx = Vx_next;
        Vy = Vy_next;

    }

    // alignment construction
    std::string align_first = "";
    std::string align_second = "";
    std::vector<double> final_probs = {Vmm[m], Vx[m], Vy[m]};
    std::pair<double, int> max_pair = max(final_probs);
    int state = max_pair.second;
    std::vector<std::vector<int>> *tracker = nullptr;
    if(state == 0){
        tracker = &tracker_mm;
    }else if(state == 1){
        tracker = &tracker_x;
    }else if(state == 2){
        tracker = &tracker_y;
    }
    int i = n;
    int j = m;
    while (i > 0 && j > 0) {
        int state = (*tracker).at(i).at(j);

        if(state == 0){
            align_first.push_back(seq1[i-1]);
            align_second.push_back(seq2[j-1]);
            i -= 1;
            j -= 1;
            tracker = &tracker_mm;
        }else if(state == 1){
            align_first.push_back(seq1[i-1]);
            align_second.push_back('-');
            i -= 1; // check if its -i or -j
            tracker = &tracker_x;
        }else if(state == 2){
            align_first.push_back('-');
            align_second.push_back(seq2[j-1]);
            j -= 1; // check if its -i or -j
            tracker = &tracker_y;
        }
    }
    /*this part requires discussion*/
    if(i == 0 && j != 0){
        while (j != 0){
            align_first.push_back('-');
            align_second.push_back(seq2[j-1]);
            j-=1;
        }
    }
    if(j == 0 && i != 0){
        while (i != 0){
            align_first.push_back(seq1[i-1]);
            align_second.push_back('-');
            i-=1;
        }
    }
    std::reverse(align_first.begin(), align_first.end());
    std::reverse(align_second.begin(), align_second.end());
    std::cout << align_first << std::endl;
    std::cout << align_second << std::endl;

}


int main(void){
    std::string seq1 = "GATTACAAAAAAAAAAAAAAAAAAATCGCTCGT";
    std::string seq2 = "GTCGACGCAATATATATAAAAAAAAAAAAAAAAAAAAAAA";
    std::pair<std::string, std::string> pair = std::make_pair(seq1, seq2);
    Model model = Model("C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\transmission_values\\estimate\\pi.txt",
    "C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\transmission_values\\estimate\\A.txt",
    "C:\\Users\\Matija\\Desktop\\BIONFO2\\BIONFORMATICS2_PROJECT\\transmission_values\\estimate\\A.txt");
    viterbi_algo(model, pair);
}