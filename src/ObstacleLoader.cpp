// esta clase pretende cargar los obstaculos generados 
/*#include <arrayfire.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
*/
#include "ObstacleLoader.h"
#include <iostream>
ObstacleLoader::ObstacleLoader(int lx, int ly , int lz){
  this->Lx = lx;
  this->Ly = ly;
  this->Lz = lz;
  object = new int**[Lx];

  for ( int i = 0; i < Lx ; i++){
    object[i] = new int*[Ly];
  
    for(int j = 0; j < Ly ; j++){
      object[i][j] = new int[Lz];
  
    }
  }
  for( int i = 0; i < Lx ; i++)
    for(int j = 0 ; j < Ly ; j++)
      for( int k = 0  ; k < Lz ; k++)
	object[i][j][k] = 0;
}
ObstacleLoader::~ObstacleLoader(){
  
  int i,j;
  for ( i = 0; i < this->Lx; i++){
    for( j = 0 ; j < this->Ly; j++){
      delete [] object[i][j];
   
    }
    delete [] object[i];
  
  }
  
  delete [] object;
 
}
af::array ObstacleLoader::loadFile(std::string file , int grid_sx, int grid_sy, int grid_sz){
  
  // vamos a leer linea por linea, primero para contar las lineas

  std::ifstream archivo(file);
  
  std::string line;
  int size = 0, index = 0;
  float * obj;
  if(archivo.is_open()){
    while( getline(archivo,line ) ){
      if(line == "") continue;
      size +=1;
    }
    // crear el arreglo de ese tamaño
    std::cout << "#tamaño: " << size << std::endl;
    archivo.close();
    std::ifstream archivo2(file);
    
    obj = new float[size];
    while( std::getline( archivo2 ,line)){
      std::stringstream linestream(line);
      std::string data;
      float v1,v2,v3,v4;
      linestream >> v1 >> v2 >> v3 >> v4;
      //std::cout << v4 << std::endl;
      obj[index] = v4;
      index++;
    }
   
    
    af::array af_obj(grid_sx, grid_sy , grid_sz, obj);
    
    // float *h_obj = af_obj.host<float>();
   
    
    /* for(int i = 0; i < 300 ; i ++)
      for( int j = 0 ; j < 100 ; j++)
	for( int k = 0 ; k < 100 ; k++) {
	  if ( obj[i + 100*(j + 300*k)] == 1) continue;
	  std::cout << i << " " << j << " " << k << " " << obj[i + 100*(j + 300*k)] << std::endl;
	  //std::cout << i + 100*(j + 300*k) << std::endl;
	}
    */
    for( int id = 0 ; id < size ; id++){
      int z = id / (grid_sx*grid_sy );
      int idp = id - (z*grid_sx*grid_sy);
      int y = idp/grid_sx;
      int x = idp%grid_sx;
      if( obj[id] == 1 ) continue;
      //std::cout << x << " " << y << " " << z << " " << obj[id] << std::endl;
    }
    return af_obj;
    
  }else { //raise exception or just say something went wrong
    
    std::cout<< "# Imposible cargar archivo" << std::endl;
    return af::constant(0, 10); // just whatever
    
  } 
  //return 0.0;
}


af::array ObstacleLoader::loadFile2(std::string file , int Lx, int Ly, int Lz){
  
  std::ifstream file_h(file);
  
  if(! file_h.is_open()) {
    std::cout<< "# Imposible cargar archivo" << std::endl;
    return af::constant(0,10);
  }else { // proceder. 
    
    // obtener la primera linea:
    std::string line;
    int aa = 0 ;
    int *h_o;
    h_o = new int[Lx*Ly*Lz];
    int i = 0;
    while( std::getline(file_h,line ) ){
      
      std::stringstream linestream(line);
      int x,y,z, o;
      linestream >> x >> y >> z >> o;
      
      //object[x][y][z] = o;
      if (  o == 1 ) {
	//std::cout << "# Puntos no-esfera:"  << z << " " << y << " " << x<< std::endl;
      }
      h_o[i] = o;
      i++;
    }  
    for ( int k = 0 ; k < Lz ; k++)
      for( int j = 0 ; j < Ly ; j++)
      for( int ix = 0 ; ix < Lx ; ix++)

	  {
	    //if ( object[ix][j][k] == 1 ) continue;
	     //std::cout << ix << " " << j << " " << k << " " << object[ix][j][k] << std::endl;
	    //std::cout << ix << std::endl;
	    //h_o[i] = object[ix][j][k];
	    //std::cout << i <<  " " << object[ix][j][k]<<std::endl;
	    //i++;
	    //if ( object[ix][j][k] != 1 && object[ix][j][k] != 0 ) return af::constant(0,10); 
	  }
   
    af::array af_obj( Lx, Ly ,Lz , h_o);
    // af_print( af_obj(150,55,69) );
     
    int *h_a = af_obj.host<int>();
    i = -1;
   
    for(int z = 0 ; z < Lz ; z++) 
      for(int y = 0 ; y< Ly  ;y++)  
	for(int x = 0 ; x < Lx ; x ++)
	  
	  {
	    //int mm = h_a[x + Ly*(y + Lx*z) ];
	    i++;
	    int mm = h_a[ i  ];
	    if ( h_a[i] !=  0  ) continue;
	    //std::cout << x << " "<< y << " " << z << " " << mm << " " << h_o[i] <<std::endl;
	      
	  } 
    //std::cout << "# i:" << i << std::endl;
    delete [] h_a ;
    delete [] h_o;
    return af_obj;
  }
}

