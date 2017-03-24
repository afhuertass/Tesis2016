#include <inc/meshIncludes.h>

#include <inc/headers/PolyhedronContainer.h>

#include <inc/headers/Modelo.h>

Modelo::Modelo(std::string path){
  
  this->loadModelo(path);
  
}
void Modelo::loadModelo(std::string path){
   Assimp::Importer importer;
   const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs);
   if(!scene || scene->mFlags == AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) // if is Not Zero
     {
       std::cout << "ERROR::ASSIMP:: " << importer.GetErrorString() << std::endl;
       return;
     }
   // en esta parte el modelo esta cargado. y guardado en scene. debemos obteneer las mesh

   this->populatePolys(scene->mRootNode ,scene);
   
}
void Modelo::populatePolys(aiNode* node, const aiScene* scene){
  // esta cosa tiene que funcionar recursicamente.
  for ( int i = 0; i < node->mNumMeshes; i++){
    // iterar sobre todas las mesh del nodo actual
    aiMesh* mesh = scene->mMeshes[node->mMeshes[i] ];
    PolyhedronContainer polC(mesh);
    polys.push_back(  polC.getPoly()  ); // un arreglo de los pol
  }
  
  for( int i = 0; i < node->mNumChildren; i++){
    
    this->populatePolys(node->mChildren[i] , scene );
    
  }
  
}
void Modelo::getOff(){
  // sacar archivos .off
  std::ofstream os;
  os.open("out.off");
  
  //for(int i = 0; i < polys.size() ; i++){
  Polyhedron P = polys.at(0);
  
  os << P;
  //}
  //std::cout << "#polys:" << polys.size() << std::endl;
  os.close();
}
Polyhedron_m Modelo::createPolyForMesh(){
  Polyhedron_m p;
  this->getOff();
  std::ifstream input("out.off");
  input >> p;

  std::cout << "#polys:" << polys.size() << std::endl;
  return p;
}
