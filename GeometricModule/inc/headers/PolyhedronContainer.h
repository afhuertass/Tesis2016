
class PolyhedronContainer{
public:
  Polyhedron P; // para guardar el resultado.
  
  PolyhedronContainer(aiMesh* oMesh); // constructor, recibimos una mesh de la Scene.
  
  Polyhedron getPoly();
private:
  std::vector<double> coords;
  std::vector<aiFace> faces;
  
};
