#include "ForceModule.h"

#include <iostream>

ForceModule::ForceModule(std::string dir){
  
  // vamos a cargar los arreglos con la data de los archivos
  /* fs = dir/fs.txt
     v = dir/vs.txt
     w = dir/w.txt
     wb = dir/wb.txt
   */
  // primero con f
  //int lx, ly, lz, qz;}
  std::cout << "# feel like a bitch" << std::endl;
  std::ifstream file_fs(dir + "/fs.txt");
  std::cout << dir + "/fs.txt" << std::endl;
  std::string fLine;
  if ( file_fs.good()){
    
    std::getline(file_fs, fLine);
    std::stringstream fLineS(fLine);
    fLineS >> Lx >> Ly >> Lz >>  q ;
     
  }else {
    std::cout << " ERROR CARGANDO ARCHIVO CON FS" << std::endl;
    return ;
  }
  // debemos pedir la memoria del arreglo f basado en el tamaño Lx,Ly,Lz,q
  //std::cout << Lx << " " << Ly << " "<< Lz << " " << q << std::endl;
  // reservar espacio para los vectores
  this->Vels.reserve(q);
  fs = new double ***[Lx];
  for( int i = 0; i < Lx ;i++){
    fs[i] = new double **[Ly];
    for(int j = 0; j < Ly ; j++){
      fs[i][j] = new double *[Lz];
      for( int k = 0 ; k < Lz ; k++){
	fs[i][j][k] = new double[q];
      }
    }
  }
  V = new int*[3];
  for( int i = 0 ; i < 3 ; i++){
    V[i] = new int[q];
  }
  
  w = new int **[Lx];
  wb = new int **[Lx];
  for ( int i = 0 ; i < Lx; i++){
    w[i] = new int *[Ly];
    wb[i] = new int *[Ly];
    for( int j = 0 ; j < Ly ; j++){
      w[i][j] = new int[Lz];
      wb[i][j] = new int[Lz];
    }
  }
 
  while( std::getline(file_fs, fLine ) ){
    std::stringstream fLineS(fLine);
    int x,y,z,q;
    double f;
    fLineS >> x >> y >> z >> q >> f ;
    fs[x][y][z][q] = f; 
  }
  file_fs.close();
  // cargar vectores 
  std::ifstream file_v(dir + "/vs.txt");
  if( !file_v.good()){
    std::cout << "ERROR CARGANDO VELOCIDADES" << std::endl;
    return ;
  }
  // aqui el archivo es bueno 
  while(std::getline( file_v , fLine)){
    std::stringstream fLineS(fLine);
    double x,y,z;
    fLineS >> x >> y >> z;
    //V[vec_coo][vec_num] = vec_value;
    //std::cout << vec_coo << " " << vec_num << " " << vec_value << std::endl;
    vector3D v; v.cargue(x,y,z);
    this->Vels.push_back(v);
  }
  file_v.close();
  
  std::ifstream file_w(dir+"/ws.txt");
  if( ! file_w.good()){
    std::cout << "ERROR CARGANDO WS" << std::endl;
    return  ;
  } // cargar los ws 
  while( std::getline( file_w, fLine ) ){
    std::stringstream fLineS(fLine);
    int x,y,z, w_value;
    fLineS >> x >> y >> z >> w_value ;
    w[x][y][z] = w_value ;
    
  }
  file_w.close();
  
  std::ifstream file_wb(dir+"/wbs.txt") ;
  if( !file_wb.good()){
    std::cout << "ERROR CARGADO WBS" << std::endl;
    return ;
  }
  
  while( std::getline( file_wb , fLine )){
    std::stringstream fLineS(fLine);
    int x,y,z,wb_value;
    fLineS >> x >> y>> z >> wb_value;
    wb[x][y][z] = wb_value;
    
  }
  file_wb.close();
  std::cout << "# feel like a bitch" << std::endl;
}

ForceModule::~ForceModule(){
  
  int i,j,k;
 
  for( i = 0 ; i< this->Lx ; i++){
    for( j =0 ; j < this->Ly; j++){
      delete [] w[i][j];
      delete [] wb[i][j];
      for(k = 0 ; k < this->Lz;k++){
	delete [] fs[i][j][k];
      }
      delete [] fs[i][j];
      
    }
    delete [] w[i];
    delete [] wb[i];
    delete [] fs[i];
  }
  
  delete [] fs;
  delete [] w ;
  delete [] wb;
  
}
vector3D ForceModule::calculateForce(){
  vector3D force;
  force.cargue(0,0,0);
  int nodes;
  for( int i = 0 ; i < Lx ; i++)
    for( int j = 0 ; j < Ly; j++)
      for( int k = 0; k < Lz ; k++){

	//w = 0 , si es fluido
	// w = 1 , si es solido
	if( w[i][j][k] == 0){ // si no es nodo frontera
	  continue;
	}else { // si es uno, esta en la frontera
	  nodes++;
	  // calcular la fuerza
	  //std::cout << "Calculando nodo frontera:"  << i << " "<< j << " "<< k<<std::endl; 
	    std::cout << "Fluid? :" << w[i][j][k] << std::endl;
	 
	  for( int alpha = 1 ; alpha < q ; alpha++){
	    double df = 0;
	    vector3D ea = this->Vels[alpha];
	    vector3D ea_ = -1*ea;
	    int in = i+ea.x(), jn = j+ ea.y() , kn = k+ea.z();
	    int alpha_ = findInverse(ea_);
	    // std::cout << "   Nodos Vecinos:" << std::endl;
	    
	    //force += ( fs[i][j][k][alpha] - fs[i][j][k][alpha_]  )*ea;
	    if ( w[ in ][ jn ][ kn ] == 0 ) { // intercambia con nodos de 
	      // fluidos 
	      df = 1;
	      df *= ( fs[i][j][k][alpha] + fs[ in ][ jn ][ kn ][alpha_] ); //*ea;
	      //std::cout  << i << " " << j << " " << k << " "<< alpha <<" "   << in << " " << jn << " " << kn << " "<< alpha_ << std::endl;
	      //std::cout << w[i][j][k] << std::endl;
	     }
	    
	    //std::cout << "df:" << df << std::endl;
	    force += df*ea;
	    
	  }
	}
      }
  std::cout << "#   Fuerza:" << norma(force) << std::endl;
  // std::cout <<  force.x() ; // fuerza x << std::endl;
  std::cout << "#    Fuerza x:" << force.x() << std::endl;
  std::cout << "#    Fuerza y:" << force.y() << std::endl;
  std::cout << "#    Fuerza z:" << force.z() << std::endl;
  std::cout << "#    nodes: " << nodes << std::endl;
  return force;

}


int ForceModule::findInverse(vector3D v){
  
  for( int i = 0; i < this->Vels.size(); i++){
    if( Vels[i].x() == v.x() && Vels[i].y()==v.y() && Vels[i].z() == v.z()){
      return i;
    }
    
  }
  
  return 0;

}
