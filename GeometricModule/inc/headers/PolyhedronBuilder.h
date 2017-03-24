#include "meshIncludes.h"

template <class HDS>
class Polyhedron_builder : public CGAL::Modifier_base<HDS>{
public:
  std::vector<double> &coords;
  std::vector<aiFace> &faces;
  Polyhedron_builder( std::vector<double> &c , std::vector<aiFace> &f) : coords(c) , faces(f){}
  void operator()(HDS& hds) {
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds,false);
    
    B.begin_surface(coords.size()/3, faces.size()  );
    //  B.begin_surface( vertices , caras esperadas, h  );
    for(int i = 0; i < (int)coords.size() ; i+=3 ){
      B.add_vertex( Point( coords[i] , coords[i+1] , coords[i+2] ));
    }
    
    for ( int i = 0; i < (int) faces.size(); i+=1){
      B.begin_facet(); // Iteramos por cada cara agregando los indicesa el facet
      for ( int j = 0 ; j < faces[i].mNumIndices; j++) {
	B.add_vertex_to_facet( faces[i].mIndices[j] );
      } 
      B.end_facet(); // fin de la cara 
    }
    B.end_surface(); // construimos la superficie 
  }
  
  
};
