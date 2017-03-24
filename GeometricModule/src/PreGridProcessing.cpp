
#include <inc/meshIncludes.h>
#include <inc/headers/MeshHelper.h>
#include <inc/headers/PreGridProcessing.h>

#include <iostream>
#include <fstream>

PreGridProcessing::PreGridProcessing(int lx, int ly, int lz , double ran){
  // constructor 
  this->Lx= lx; this->Ly=ly;
  this->Lz = lz;
  this->range = ran;
  
  grid = new int**[Lx];
  gridb = new int**[Lx];
  for ( int i = 0; i < Lx ; i++){
    grid[i] = new int*[Ly];
    gridb[i] = new int*[Ly];
    for(int j = 0; j < Ly ; j++){
      grid[i][j] = new int[Lz];
      gridb[i][j] = new int[Lz];
    }
  }
}
PreGridProcessing::~PreGridProcessing(){
  
  int i,j;
  for ( i = 0; i < this->Lx; i++){
    for( j = 0 ; j < this->Ly; j++){
      delete [] grid[i][j];
      delete [] gridb[i][j];
    }
    delete [] grid[i];
    delete [] gridb[i];
  }

  delete [] grid;
  delete [] gridb;
}

void PreGridProcessing::GetSolidRegion(MeshHelper &mh) {
  // here comes the sun :v 
  std::cout << "Obteniendo region solida " << std::endl;
  double rx_min = 0; double rx_max = 15;
  double ry_min = 0; double ry_max = 6;
  double rz_min = -1.5; double rz_max = 1.5;
  
  double deltaX = (rx_max - rx_min)/this->Lx;
  double deltaY = (ry_max - ry_min)/this->Ly;
  double deltaZ = (rz_max - rz_min)/this->Lz;
  double x0,y0,z0 ;
  int ix,iy,iz;
  x0 = rx_min;
  y0 = ry_min;
  z0 = rz_min;
  for( int iz  = 0 ; iz < this->Lz; iz++ , z0+=deltaZ){
    for( int iy  = 0 ; iy < this->Ly; iy++ , y0+=deltaY){
      for( int ix = 0 ; ix < this->Lx ; ix++ , x0+=deltaX){
	
	Point_3_k punto(x0,y0 ,z0);
	//std:: cout << x0 << " " << y0 << " " << z0 << std::endl;
	if( mh.isPointIn( punto ) ){
	  
	  grid[ix][iy][iz] = 0; // es solido no habra velocidad
	  
	}else {
	  grid[ix][iy][iz] = 1; // fluido 
	}
	
      }
      x0 = rx_min;
    }
    y0 = ry_min;
  }
  std::ofstream file;
  std::ofstream w;
  file.open("watch.dat");
  w.open("ws.txt");
  //file << Lx << " " << Ly << " " << Lz << "\n";
  for( int iz  = 0 ; iz < this->Lz; iz++ ){
    for( int iy  = 0 ; iy < this->Ly; iy++){
      for( int ix = 0 ; ix < this->Lx ; ix++){
	if( grid[ix][iy][iz] == 1){
	  w << ix << " " << iy << " " << iz << " " << 0 << "\n";
	} else { // nodo solido  w = 1
	  w << ix << " " << iy << " " << iz << " " << 1 << "\n";
	}
	
	if ( grid[ix][iy][iz] == 1 ) continue;
	file << ix << " " <<  iy << " " << iz << " " << grid[ix][iy][iz] << "\n";
	
      }}}
  file.close();
  w.close();

}
void PreGridProcessing::GetBoundaryRegion(){
  // determinar la region 
  int i,j,k ;
  int i2,j2,k2;
  i2 = j2 = k2 = 0;
  std::cout << " Obteniendo Boundary region" << std::endl;
  for ( i = 0; i < this->Lx ; i++){
    for( j = 0 ; j < this->Ly; j++){
      for( k = 0; k < this->Lz ; k++){
	double sum = 0;
	// por cada  celda hay que checkear los vecinos. 
	if( grid[i][j][k] == 0 ){ // si es un nodo solido 
	  gridb[i][j][k] = 0; // no me interesa
	  continue; 
	} //
	  // me interesa los nodos fluidos y que estan en la frontera
	for( i2 = i-1 ; i2 <= i+1 ; i2++){
	  for(j2 = j-1; j2 <= j+1 ; j2++){
	    for(k2 = k-1 ; k2 <= k+1 ; k2++){
	      //std::cout << "isss: " << (i+i2+Lx)%Lx  << " " <<  (j+j2+Ly)%Ly << " " <<  (k+k2+Lz)%Lz<< std::endl; 
	      //std::cout << "is: " << i2  << " " << j2 << " " <<  k2 << std::endl; 
	      if( grid[ (i2+Lx)%Lx ][ (j2+Ly)%Ly ][ (k2+Lz)%Lz ] == 0 ){ // checkear es solido 
		sum++; 
		//std::cout << sum << std::endl;
	      }
	      
	    }
	  }
	}
	
	if( sum > 0  ){// el nodo solido, tiene por lo menos un vecino fluido 
	  gridb[i][j][k] = 1; // es un nodo frontera
	}else {
	  gridb[i][j][k] = 0; // no es frontera 
	}
	
      }
    }
  }
  std::ofstream file;
  file.open("wbs.txt");
  for(  i  = 0 ; i < this->Lx; i++ ){
    for(  j  = 0 ; j < this->Ly; j++){
      for( k = 0 ; k < this->Lz ; k++){
	//if ( gridb[i][j][k] == 0) continue;
	file << i << " " <<  j << " " << k << " " << gridb[i][j][k] << "\n";
      }}}
  file.close();


}
