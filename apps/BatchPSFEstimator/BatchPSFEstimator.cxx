#include <DataModel.h>

#include <iostream>
#include <sstream>
#include <string>


int main(int argc, char* argv[]) {
  DataModel* model = new DataModel();

  if (argc < 2) {
    std::cout << "Usage: BatchPSFOptimizer <VPO settings file name>" << std::endl;
    return 1;
  }

  model->LoadSessionFile(std::string(argv[1]));
  model->Optimize();
  
  // Save the results to a different file with a modified name
  std::stringstream ss;
  ss << argv[1] << "-optimized.psfe";
  model->SaveSessionFile(ss.str());


  delete model;
}
