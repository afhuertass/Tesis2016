#include <arrayfire.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class ObstacleLoader{
 private:  
  int ***object;
  int Lx, Ly, Lz;
 public:
  ObstacleLoader(int Lx, int Ly, int Lz);
  ~ObstacleLoader();
  af::array loadFile(std::string file, int grid_sx , int grid_sy , int grid_sz  );
  af::array loadFile2(std::string file , int grid_sx, int grid_sy, int grid_sz);
};

