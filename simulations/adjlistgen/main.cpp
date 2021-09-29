#include <iostream>
#include <iomanip>
#include "adjlistgen.hpp"

int main(int argc, char **argv) {
    uvector flatAdjList,offsetList;
    std::vector<float> weightList;
       
    WeightedListGen wl("/tmp/edge_list.txt");
    wl.load_weighted_adjacency_list();
    wl.copy_weighted_adjacency_list(flatAdjList,offsetList,weightList);
    for(int s = 0; s<offsetList.size()-1; s++){
      for(int neighbor=0; neighbor< (offsetList[s+1] - offsetList[s]); neighbor++){
	std::cout <<std::setprecision(4)<<std::fixed<< s << "\t"<< flatAdjList[ offsetList[s] + neighbor] << "\t" << weightList[offsetList[s] + neighbor] << "\n";
    }
  }
}