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
#include <inc/headers/Model.h>
#include <inc/headers/MeshHelper.h>
#include <inc/headers/PreGridProcessing.h>

int main(){
  // it builds up the model 
  Modelo casa("cubo.obj");
  Polyhedron_m p = casa.createPolyForMesh();
  
  MeshHelper mh( p );
  
  mh.saveMeshFile();

  PreGridProcessing preGrid(300,120,60,2.5);
  
  preGrid.GetSolidRegion(mh);
  
  preGrid.GetBoundaryRegion();
  return 0;
}
