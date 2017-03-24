class Modelo {
  // think el modelo va a tener los polyedros 
public:
  std::vector<Polyhedron> polys;
  Modelo(std::string path); // constructor, recibe la ruta para cargar el modelo 
  void getOff();
  Polyhedron_m createPolyForMesh();
private:
  void loadModelo(std::string path);
  void populatePolys(aiNode* nodo ,const aiScene* scene);
};

