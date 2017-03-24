#include <inc/meshIncludes.h> 

#include <inc/headers/MeshHelper.h>
#include <assert.h>

using namespace CGAL::parameters;

MeshHelper::MeshHelper(Polyhedron_m &p){
  // en el constructor generamos el mesh
  
  //this->generateMesh(p);
  domain = new Mesh_domain(p);
  
  Mesh_domain doo(p);

}
MeshHelper::~MeshHelper(){
  //delete [] this->domain;
}
void MeshHelper::generateMesh(Polyhedron_m &p){
  Mesh_domain domain(p);
  std::cout << "#wtf?" << std::endl;
  Mesh_criteria criteria(facet_angle=25, facet_size=0.1, facet_distance=0.008,cell_radius_edge_ratio=2 , cell_size=1);
 std::cout << "#wtf 2 ?  mesh_criteria created" << std::endl;
  // criterios, para ensayos, probamos con estos y luego extendemos.

  this->mesh = CGAL::make_mesh_3<C3t3>( domain , criteria , no_perturb(), no_exude()) ;
   std::cout << "#wtf? 3  mesh created" << std::endl;

}
void MeshHelper::saveMeshFile(){
  
  std::ofstream medit_file("out.mesh");
  this->mesh.output_to_medit(medit_file);
  medit_file.close();

};
void MeshHelper::saveBoundary(){
  std::ofstream off_file("boundary.off");
  this->mesh.output_boundary_to_off(off_file);
  off_file.close();
  
}
bool MeshHelper::isPointIn(const Point_3_k  &query){
  // 
  /*Locate_type lt;
  int ii, jj;
  Tr &tr = mesh.triangulation();
  //Cell_handle infinite = tr.infite_cell();
  Cell_handle ch = tr.locate( query , lt , ii ,jj );
  
  return ( this->mesh.is_in_complex( ch  ) );
  // return ( this->domain->is_in_domain(query) );
  */
  typedef typename Mesh_domain::Is_in_domain Is_in_domain;
  typedef typename Mesh_domain::Subdomain Subdomain;
  typedef typename Mesh_domain::Subdomain_index Subdomain_index;
  Is_in_domain is_in_domain = domain->is_in_domain_object();
  Subdomain ss = is_in_domain( query );
  if ( ss ){
    return true;
  }else {
    return false;
  }
  
}
