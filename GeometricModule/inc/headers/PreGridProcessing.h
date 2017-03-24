class PreGridProcessing{
  /* bien esta clase se encargara de poner la malla uno a , entonces 
     que reciba el meshhelper y de ahi lo convierte en 2 archivos de texto plano
que contienen: 
los nodos que son solidos. 
los nodos que son frontera.

en dos arreglos que seran cargados y usados a posterior por la simulacion 
necesitamos 3 distancias caracterizticas 
  */
 public:
  PreGridProcessing(int Lx, int Ly, int Lz, double range); // constructor, tamaño de las
  ~PreGridProcessing();
  void GetSolidRegion(MeshHelper & mh);
  void GetBoundaryRegion(); 
  // boundary region solo puede ser hallada 
 private:
  int Lx,Ly,Lz;
  int *** grid;
  int *** gridb;
  double range;
  

};
