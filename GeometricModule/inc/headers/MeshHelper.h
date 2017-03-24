class MeshHelper{
  // clase para generar Mesh3D a partir de polyhedro
 public:
  MeshHelper(Polyhedron_m &p);
  ~MeshHelper();
  void saveMeshFile();
  void saveBoundary();
  bool isPointIn(const Point &p );
 private:
  C3t3 mesh;
  Mesh_domain *domain;
  void generateMesh(Polyhedron_m &p);
  
  
};
