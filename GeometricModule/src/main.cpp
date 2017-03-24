/*#include "meshIncludes.h"

#include "polyhedronContainer.h"
#include "Modelo.h"
#include "meshHelper.h"
*/
// #include "classes.h"
// #define CGAL_DISABLE_ROUNDING_MATH_CHECK 

#define BOOST_PARAMETER_MAX_ARITY 12
// remover en un futuro 


#include <inc/meshIncludes.h>
#include <inc/headers/Modelo.h>
#include <inc/headers/MeshHelper.h>
#include <inc/headers/PreGridProcessing.h>

int main(){

  Modelo casa("cubo.obj");
  Polyhedron_m p = casa.createPolyForMesh();
  //std::cout << "#size of Model " << sizeof(casa) << std::endl;
  //std::cout << "#size of Model " << sizeof(p) << std::endl;
  MeshHelper mh( p );
  //std::cout << "#size of MeshHelper " <<sizeof(mh) <<  std::endl;
  mh.saveMeshFile();

  PreGridProcessing preGrid(300,120,60,2.5);
  
  preGrid.GetSolidRegion(mh);
  
  preGrid.GetBoundaryRegion();
  return 0;
}
