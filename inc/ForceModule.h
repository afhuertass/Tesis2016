
// Clase para calcular la fuerza sobre el objeto utilizando el metodo
// de intercambio de momentum asi que necesitamos:
// el arreglo f , las velocidades , los arreglos de las velocidades
// el arreglo de que puntos del lattice son fronteras y cuales son solidos.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
//#include "Vector.h" // vector matematico
#include <vector> // contenedor
#include "Vector.h" // vector fisico

#include <algorithm>

class ForceModule {
  public:
  ForceModule(std::string dir);
  ~ForceModule();
  void reLoadElements(std::string dir);
  vector3D calculateForce();
  
 private:
  double ****fs;
  int **V; // velocidades del lattice
  int ***w;
  int ***wb; 
  int Lx,Ly,Lz,q;
  std::vector<vector3D> Vels;
  int findInverse( vector3D v);
};
