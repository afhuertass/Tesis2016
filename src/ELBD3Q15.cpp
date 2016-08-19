
#include <arrayfire.h>
#include <ELB3DQ15.h>
#include <ObstacleLoader.h>
#include <ForceModule.h>
#include <vector>

// Modelo 3d para lattice bolztman utlizando el metodo entrpico 
using namespace af;
ELatticeBolztmannD3Q15::ELatticeBolztmannD3Q15(int lx,int ly, int lz){
  //constructor
  this->q = 15;
  this->Lx=lx; // tamaño del lattice
  this->Ly=ly;
  this->Lz=lz;
  fs = std::vector<array>(q);
  
  for( int i = 0 ; i < q ; i++){
    this->fs[i] = af::constant(0,Lx,Ly,Lz , f32);
  }
  w = af::constant( 0 , q , f32);
  alpha = af::constant( 1.1 , Lx, Ly,Lz , f32);
  
  this->ConditionsU = af::constant(1, Lx,Ly,Lz,f32);

  w(0) = 16; w(1)=w(2)=w(3)=w(4)=w(5)=w(6) = 8;
  w(7)=w(8)=w(9)=w(10)=w(11)=w(12)=w(13)=w(14)=1;
  w = 1/72.0*w;
  W[0] = 16.0/72;
  W[1] = W[2] = W[3] =   W[4] = W[5] = W[6] = 8.0/72;
  W[7] = W[8] = W[9] =   W[10] = W[11] = W[12] = W[13] = W[14] = 1.0/72; 
  
  V[0][0] = 0; V[0][1] = 1; V[0][2] = 0; V[0][3] = 0 ; V[0][4] = -1;
  V[0][5] = 0; V[0][6] = 0; V[0][7] = 1; V[0][8] = -1 ; V[0][9] = 1;
  V[0][10] = 1; V[0][11] = -1; V[0][12] = 1; V[0][13] = -1 ; V[0][14] = -1;
  
  V[1][0] = 0; V[1][1] = 0; V[1][2] = 1; V[1][3] = 0 ; V[1][4] = 0;
  V[1][5] = -1; V[1][6] = 0; V[1][7] = 1; V[1][8] = 1 ; V[1][9] = -1;
  V[1][10] = 1; V[1][11] = -1; V[1][12] = -1; V[1][13] = 1 ; V[1][14] = -1;

  V[2][0] = 0; V[2][1] = 0; V[2][2] = 0; V[2][3] = 1 ; V[2][4] = 0;
  V[2][5] = 0; V[2][6] = -1; V[2][7] = 1; V[2][8] = 1 ; V[2][9] = 1;
  V[2][10] = -1; V[2][11] = 1; V[2][12] = -1; V[2][13] = -1 ; V[2][14] = -1;

 
}
void ELatticeBolztmannD3Q15::Inicie(float r0 , float Ux0 , float Uy0 , float Uz0){
  /*array rh0 = af::constant( r0 , Lx, Ly, Lz , f32);
  array Uxs = af::constant( Ux0, Lx , Ly , Lz , f32);
  array Uys = af::constant( Uy0, Lx , Ly , Lz , f32);
  array Uzs = af::constant( Uz0, Lx , Ly , Lz , f32);
  */
  array rh0 = af::randu( Lx, Ly, Lz , f32)/100.0;
  array Uxs = af::randu( Lx , Ly , Lz , f32)/100.0;
  array Uys = af::randu( Lx , Ly , Lz , f32)/100.0;
  array Uzs = af::randu( Lx , Ly , Lz , f32)/100.0;
  
  //f = this->feq( rh0, Uxs , Uys, Uzs);
  SetConditions( Uxs , Uys, Uzs );
  int i = 0;
  for( i = 0; i< q ; i++ ) {
    this->fs[i] = this->feq2(i,rh0 , Uxs, Uys, Uzs ) ;
  }
  float Re = 100;
  float ulb = 0.03;
  float r = 25;
  float nulb = ulb*r/Re; // viscosidad de referencia 

  this->beta = 1/3.0*(nulb + 0.5); 
  
  /*double Re = 100 ; // numero de Reynolds
  double viscosidad = (r0*0.05*Lx)/Re; // rho*u*d/Re
  double cs2 = 1.0/3;
  double tau = viscosidad/cs2;
  
  double dx = 1.0/Lx;
  double dt = dx*dx;
  */
 
  std::cout << "#beta:" << beta << std::endl;
  std::cout << "# kinematic viscosity: " << 1.0/3*(1.0/beta - 1.0/2)<<std::endl; 
}
void ELatticeBolztmannD3Q15::SetConditions(array &Ux, array &Uy , array &Uz){
  
// poiseullie 
 // plano x = 0 la velocidad es maxima
  
 //  Ventiladores en la entrada
  //  Ux(0,span,span ) = 0.03; // this->UMAX;
  //  Uy(0,span,span ) = 0;
  //  Uz(0,span,span ) = 0;
  
  // plano y = 0
  // lid driven conditions : 
  Uy(span , 0 , span) = 0;

  // plano y = Ly , velocidad
  
  Uy(span , Ly-1 , span) = 0;

  // plano x = 0

  Ux(0,span,span ) = 0.0;
  Uy(0,span,span ) = 0;
  Uz(0,span,span ) = 0;
  

   // plano x = Lx
  Ux(Lx-1,span,span ) = 0.0;
  Uy(Lx-1,span,span ) = 0;
  Uz(Lx-1,span,span ) = 0;
  
  // plano z = 0
  Ux(span,span,0 ) = 0.0;
  Uy(span,span, 0) = 0;
  Uz(span,span,0 ) = 0;
 
   // plano z = Lz
  Ux(span,span,Lz-1 ) = 0.03;
  Uy(span,span, Lz-1) = 0;
  Uz(span,span, Lz-1) = 0;
  
}
void ELatticeBolztmannD3Q15::setArrayC(array c){
  this->ConditionsU = c;
}
array ELatticeBolztmannD3Q15::feq2(int i ,  array &rhos, array &Uxs , array &Uys , array &Uzs){ 
  // 0 x , 1 y , 2 , z 
  //std::cout << W[i] << std::endl;
  return rhos*W[i]*aux(i,0, Uxs)*aux(i,1,Uys)*aux(i,2,Uzs);
}

array ELatticeBolztmannD3Q15::rho(){
  array f = af::constant(0, Lx, Ly,Lz , f32 );
  
  //return (af::sum(f,3)); // suma a lo largo de la tercera componente
  for( int i = 0 ; i < q ; i++){
    f += this->fs[i];
  }
  return f;
}
array ELatticeBolztmannD3Q15::Jx(){
  array Jx = af::constant(0 , Lx, Ly,Lz, f32);
  
  Jx += fs.at(1) -fs.at(4) + fs.at(7)-fs.at(8)+fs.at(9)+fs.at(10)-fs.at(11)+fs.at(12)-fs.at(13)-fs.at(14) ;
  
  return Jx;
}
array ELatticeBolztmannD3Q15::Jy(){
  array Jy = af::constant(0 , Lx, Ly,Lz, f32);
  Jy += fs.at(2)-fs.at(5)+fs.at(7)+fs.at(8)-fs.at(9)+fs.at(10)-fs.at(11)-fs.at(12)+fs.at(13)-fs.at(14); 
  return Jy;
}
array ELatticeBolztmannD3Q15::Jz(){
 array Jz = af::constant(0 , Lx, Ly,Lz, f32);
 Jz+= fs.at(3)-fs.at(6)+fs.at(7)+fs.at(8)+fs.at(9)-fs.at(10)+fs.at(11)-fs.at(12)-fs.at(13)-fs.at(14); 
  return Jz;
}
array  ELatticeBolztmannD3Q15::aux(int i,int alpha, array &Ua){
  array u = af::constant(1, Lx,Ly,Lz, f32);
  array U2 = af::pow(Ua,2);
  array term1 =  (2*u-af::sqrt( 1+ 3*U2 ) );
  if ( V[alpha][i] == 0 ){
    return term1;
  }else {
    array term2 = af::pow((2*Ua+af::sqrt(1+3*U2 ) )/(u-Ua) , V[alpha][i]);
    return term1*term2;
  }
		  
 
}

void ELatticeBolztmannD3Q15::Colission(){
  float rex = 1/0.53; // tau = 0.53
  array rho = this->rho();
  array Ux = (this->Jx()/rho)*this->ConditionsU ;
  array Uy = (this->Jy()/rho)*this->ConditionsU ;
  array Uz = (this->Jz()/rho)*this->ConditionsU ;
  
  SetConditions( Ux , Uy, Uz );
  std::vector<array> feqs; // potencialmente lento 
  feqs = std::vector<array>(q);
  
  CalculateFeq( rho ,Ux, Uy , Uz , feqs);
  
  CalculateAlpha ( feqs  );

  for( int i = 0; i < q ; i++){
    fs[i] =  fs[i] + 2*beta*( feqs[i] - fs[i] ); 
  } 
  
}

void ELatticeBolztmannD3Q15::Adveccion(){
  for(int i = 0; i< q ; i++){
    array fp = this->fs[i];
   
    fs[i] = shift( fp , V[0][i] , V[1][i] , V[2][i] , 0 ) ;
   
  }
  
 
   
}

void ELatticeBolztmannD3Q15::SaveFs(std::string route_fs){
  
  std::ofstream file_fs(route_fs+"fs.txt");
  if( !file_fs.is_open() ) {
    std::cout << "ERROR ABRIENDO FS" << std::endl;
    return ;
  }
  file_fs << Lx << " " << Ly << " " << Lz << " " << q << "\n";
  float *fi;
  for( int alpha = 0 ; alpha < q ; alpha++){
    fi = this->fs[alpha].host<float>();
    for(int x = 0 ; x < Lx ; x ++)
      for(int y = 0 ; y< Ly  ;y++ ) 
	for(int z = 0 ; z < Lz ; z++){
	  float fval = fi[x + Ly*(y + Lx*z) ];
	  file_fs << x << " " << y << " " << z << " " << alpha << " " << fval <<"\n"; 
	  
	}
  }
  file_fs.close();
  // Guardando velocidades
  std::ofstream file_vs(route_fs+"vs.txt");
  if( !file_vs.is_open()){
    std::cout << "ERROR ABRIENDO VS" << std::endl;
    return;
  }
  
  for( int alpha = 0 ; alpha < q ; alpha++){
    file_vs << V[0][alpha] << " " << V[1][alpha]<< " " << V[2][alpha] << "\n";
  }
  file_vs.close();
  
  delete[] fi;
  
}

void ELatticeBolztmannD3Q15::saveVTK(std::string route){
  // guardar un legacy file para
  array rho = this->rho();
  array Ux = (this->Jx()/rho)*this->ConditionsU; 
  array Uy = (this->Jy()/rho)*this->ConditionsU;
  array Uz = (this->Jz()/rho)*this->ConditionsU;
  std::ofstream file_vtk( route );
  if( !file_vtk.is_open()){
    std::cout << "# error abriendo archivo " << route << std::endl; 
    return ;
  }
  file_vtk << "# vtk DataFile Version 3.0 \n";
  file_vtk << "vtk output \n ";
  file_vtk << "ASCII \n";
  file_vtk << "DATASET STRUCTURED_GRID \n";
  file_vtk << "DIMENSIONS " << Lx << " " << Ly << " " <<  Lz <<"\n";
  file_vtk << "POINTS " << Lx*Ly*Lz << " float \n";
  for(int k = 0; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++)
      for( int i = 0 ; i < Lx ; i++) 
	file_vtk << i << " " << j << " " << k << "\n";

  // hecha la grilla... a por los valores de las velocidades
  file_vtk << "POINT_DATA " << Lx*Ly*Lz << "\n";

  /// Ux field 
  float *h_ux = Ux.host<float>();
  file_vtk << "SCALARS Ux float \n";
  file_vtk << "LOOKUP_TABLE default \n";
  int index = 0;
  for( int k = 0 ; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++)
      for(int i = 0; i < Lx ; i++) { 
	file_vtk << h_ux[ index  ] << " \n";
	index++;
      }
 
  
  float *h_uy = Uy.host<float>();
  index = 0;
  file_vtk << "SCALARS Uy float \n";
  file_vtk << "LOOKUP_TABLE default \n";
  for( int k = 0 ; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++) 
      for( int i = 0 ; i < Lx ; i++){ 
	file_vtk << h_uy[ index ] << "\n";
	index++;
      }
  
  float *h_uz = Uz.host<float>();
  index = 0;
  file_vtk << "SCALARS Uz float \n";
  file_vtk << "LOOKUP_TABLE default \n";
  for( int k = 0 ; k < Lz ; k++)
    for( int i = 0 ; i < Lx ; i++)
      for( int j = 0 ; j < Ly ; j++)
	{
	  file_vtk << h_uz[ index ] << "\n";
	  index++;
      }
  

  file_vtk << "VECTORS VelocityField float \n";
  index = 0;
  for( int k = 0 ; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++)
      for( int i = 0 ; i < Lx ; i++)
      {
	file_vtk << h_ux[ index ] << " " << h_uy[index ] << " " << h_uz[ index ] << " \n";
	index++;
      }
  
  file_vtk.close();
  delete [] h_ux;
  delete [] h_uy;
  delete [] h_uz;
}
void ELatticeBolztmannD3Q15::CalculateAlpha(std::vector<array> &feq) {

  // paso 1 calcular la entropia de cada celda:

  array ent = af::constant( 0 , Lx , Ly , Lz , f32);
  array ent_neq = af::constant( 0 , Lx , Ly , Lz , f32);
  array delta_ent_derivative = af::constant( 0 , Lx , Ly , Lz , f32);
  array II = af::constant( 1 , Lx , Ly ,  Lz , f32);
  array tags = af::constant( 0 , Lx , Ly , Lz , f32);
  
  for ( int i = 0 ; i < q ; i++){
    ent += fs[i]*af::log(fs[i]/W[i]);
  }
  // pasos posibles calcular unos maximos para alpha
  // implementacion burda 
  for( int k = 0 ; k < 5 ; k++){

    
    for( int i = 0 ; i < q ; i++ ) {
      array h = fs[i] + alpha*(feq[i]-fs[i] );
      
      array t = af::log( h/W[i]);
      
      ent_neq += h*t;
      delta_ent_derivative += (feq[i] - fs[i])*(t + II);
     
    }
    

    // 1.1 en el alpha donde sea nan, el punto es que en dicha celda no se debe hacer newton

    alpha = 1.13*(   af::isNaN( ent_neq  )) + alpha*( !af::isNaN(ent_neq ) );
    //tags = ( af::abs(ent_neq - ent) < 1e-2  ) ;  // no requieren calcular alpha
    

    
    //alpha = 1.5*( (af::abs(ent_neq - ent)   )< 1e-2 ) + (alpha - ((ent_neq - ent)/delta_ent_derivative ))*(   af::abs(ent_neq - ent) > 1e-2  );
    // 
    af_print( alpha );
    float cc = af::sum<float>( alpha(span,span,span) >= 2  );
    if( cc >= Lx*Ly*Lz ){
      std::cout << " Todos son 2, iteracion: " << k  << std::endl;
      break; 
    }
    ent_neq = af::constant(0 , Lx , Ly , Lz , f32);
    delta_ent_derivative = af::constant(0 , Lx , Ly , Lz , f32);
   
  }
} 
void ELatticeBolztmannD3Q15::CalculateFeq( array & rhos, array & Uxs, array &Uys , array &Uzs , std::vector<array> &feqs){
  for( int i = 0 ; i < q ; i++){
    //
    feqs[i] = this->feq2( i, rhos,  Uxs , Uys , Uzs );
  }
  
  
}
array ELatticeBolztmannD3Q15::Falpha(std::vector<array> & feq){
  array Fa = af::constant( 0 , Lx , Ly , Lz , f32);
  
  for( int i = 0 ; i < q  ; i++){
    Fa += ( fs[i] + alpha*(feq[i] - fs[i]) )*af::log( (fs[i]+alpha*( feq[i]- fs[i] ) )/W[i] ) - fs[i]*af::log( fs[i]/W[i]); 
  }
  
  return Fa;
}

array ELatticeBolztmannD3Q15::DFalpha(std::vector<array> & feq){
  array DFa = af::constant(  0 , Lx , Ly , Lz , f32);
  array II = af::constant(1 , Lx , Ly , Lz , f32);
    for( int  i = 0 ; i < q ; i++){
      DFa += ( W[i]*II + af::log( fs[i] + alpha*(feq[i] - fs[i]) )  )*(feq[i] - fs[i] );  
  }
  
  return DFa;
}
