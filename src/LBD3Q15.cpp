
#include <arrayfire.h>
#include <LBD3Q15.h>
#include <ObstacleLoader.h>
#include <ForceModule.h>
#include <vector>
#include <cmath>
// Modelo 3d para lattice bolztman utlizando el metodo entrpico 
using namespace af;
LatticeBolztmannD3Q15::LatticeBolztmannD3Q15(int lx,int ly, int lz){
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

  this->ConditionsU = af::constant(0, Lx,Ly,Lz,f32);

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
  
  int vx[] = {0,1,0,0,-1,0,0,1,-1,1,1,-1,1,-1,-1  }; 
  int vy[] = {0,0,1,0,0,-1,0,1,1,-1,1,-1,-1,1,-1  };
  int vz[] = {0,0,0,1,0,0,-1,1,1,1,-1,1,-1,-1,-1 };
  
  vel_x = array( q , vx);
  vel_y = array (q , vy);
  vel_z = array(q , vz);
  
  rhos = array(Lx, Ly ,Lz , 1);
  Ux = array(Lx, Ly ,Lz , 1);
  Uy = array(Lx, Ly ,Lz , 1);
  Uz = array(Lx, Ly ,Lz , 1);
}
void LatticeBolztmannD3Q15::Inicie(float r0 , float Ux0 , float Uy0 , float Uz0){
  rhos = af::constant( r0 , Lx, Ly, Lz , f32);
  Ux = af::constant( Ux0, Lx , Ly , Lz , f32);
  Uy = af::constant( Uy0, Lx , Ly , Lz , f32);
  Uz = af::constant( Uz0, Lx , Ly , Lz , f32);

  //f = this->feq( rh0, Uxs , Uys, Uzs);
  int i = 0;
  for( i = 0; i< q ; i++ ) {
    // this->fs[i] = this->feq2(i,rh0 , Uxs, Uys, Uzs ) ;
    this->fs[i] = feq ( i );
  }
  float Re = 100;
  float ulb = 0.05;
  float r = Lx; // longitud caracteristica
  float nulb = ulb*r/Re; // viscosidad.
  this->beta = 1/(3.0*nulb + 0.5); 
  this->UMAX = 0.05;
  
}
void LatticeBolztmannD3Q15::setArrayC(array c){
  
  this->ConditionsU = c;
}
array LatticeBolztmannD3Q15::feq2(int i ,  array &rhos, array &Uxs , array &Uys , array &Uzs){ 
  // 0 x , 1 y , 2 , z 
  //std::cout << W[i] << std::endl;
  return rhos*W[i]*aux(i,0, Uxs)*aux(i,1,Uys)*aux(i,2,Uzs);
}
array LatticeBolztmannD3Q15::feq(int i ){ 
  /* SLOW IMPLEMENTATION
  array UpVi = Ux*V[0][i]+Uy*V[1][i] + Uz*V[2][i];
  array U2 = Ux*Ux + Uy*Uy + Uz*Uz;
  return rhos*W[i]*(1 + 3*UpVi + 9.0/2*af::pow( UpVi,2) - 3.0/2*U2 );
  */
  array UpVi = Ux*V[0][i]+Uy*V[1][i] + Uz*V[2][i];
  UpVi = af::flat( UpVi);
  array U2 = Ux*Ux + Uy*Uy + Uz*Uz;
  U2 = af::flat( U2);
  array feq_flat = ( W[i]*af::flat(rhos)*( 1 + 3*UpVi + 9.0/2*af::pow( UpVi,2) - 3.0/2*U2 )  );
  //return rhos*W[i]*(1 + 3*UpVi + 9.0/2*af::pow( UpVi,2) - 3.0/2*U2 );
  array feq( feq_flat , fs[0].dims() );
  return feq;
}
array LatticeBolztmannD3Q15::feq3(array &rhos, array &Uxs , array &Uys, array &Uzs ){
    
  //return (  af::tile( rhos ,  q )*af::tile( w , 1 , Lx , Ly , Lz ))*aux2(vel_x, Uxs)*aux2(vel_y,Uys)*aux2(vel_z,Uzs);
  // af::moddims( w , 1 , 1 , 1 ,q);
  array w_mod = af::tile( af::moddims( w , 1 , 1 , 1 ,q) , Lx , Ly , Lz);
  // completo: af::tile( rhos , 1 , 1, 1, q)*w_mod
  return  af::tile( rhos , 1 , 1, 1, q)*w_mod*aux2( vel_x , Uxs)*aux2( vel_y , Uys)*aux2(vel_z,Uzs) ;
} 
array LatticeBolztmannD3Q15::aux2(array &vels , array &Ua ){
  //array f2v = af::tile( (2*Ua+af::sqrt(1+3*(Ua*Ua) ))/(1-Ua) , 1 , 1 ,1 , q );
  array f2v = af::tile( (2*Ua+af::sqrt(1+3*(Ua*Ua) ))/(1-Ua) , 1 , 1 ,1 , q );
  array vels_mod = af::tile( af::moddims( vels , 1 , 1 , 1 ,q) , Lx , Ly , Lz);
  array f1v = af::tile( 2-af::sqrt( 1 + Ua*Ua) , 1, 1, 1, q  ) ;
  array res = af::pow( f2v ,  vels_mod  ); 
  return  f1v*res;

}
array  LatticeBolztmannD3Q15::aux(int i,int alpha, array &Ua){
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
array LatticeBolztmannD3Q15::rho(){
  array f = af::constant(0, Lx, Ly,Lz , f32 );
  
  //return (af::sum(f,3)); // suma a lo largo de la tercera componente
  for( int i = 0 ; i < q ; i++){
    f += this->fs[i];
  }
  return f;
}
array LatticeBolztmannD3Q15::Jx(){
  return ( fs.at(1) -fs.at(4) + fs.at(7)-fs.at(8)+fs.at(9)+fs.at(10)-fs.at(11)+fs.at(12)-fs.at(13)-fs.at(14) ) ;
  
}
array LatticeBolztmannD3Q15::Jy(){
  return ( fs.at(2)-fs.at(5)+fs.at(7)+fs.at(8)-fs.at(9)+fs.at(10)-fs.at(11)-fs.at(12)+fs.at(13)-fs.at(14) ); 
 
}
array LatticeBolztmannD3Q15::Jz(){
  return ( fs.at(3)-fs.at(6)+fs.at(7)+fs.at(8)+fs.at(9)-fs.at(10)+fs.at(11)-fs.at(12)-fs.at(13)-fs.at(14) ); 
 
}
void LatticeBolztmannD3Q15::Colission(){
  // tau = 0.53
  
  for( int i = 0; i < q ; i++){
    fs[i] =  fs[i] + beta*( feq( i ) - fs[i] ); 
  } 
  
}

void LatticeBolztmannD3Q15::Adveccion(){
   for( int i = 0 ; i < q ; i++ ){
    array fp = this->fs[i];
    fs[i] = shift( fp , V[0][i] , V[1][i] , V[2][i] , 0 ) ;
  }
  
  //af_print(this->fs[j] );
   
}
void LatticeBolztmannD3Q15::MacroscopicQuantities(){
  
  rhos = this->rho();
  Ux = (this->Jx()/rhos) ;
  Uy = (this->Jy()/rhos);
  Uz = (this->Jz()/rhos);
  
}

void LatticeBolztmannD3Q15::SaveFs(std::string route_fs){
  
  std::ofstream file_fs(route_fs+"/fs.txt");
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
  std::ofstream file_vs(route_fs+"/vs.txt");
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

void LatticeBolztmannD3Q15::saveVTK(std::string route){
  // guardar un legacy file para
  this->BounceBackBoundaries();
  this->MacroscopicQuantities();
  //this->SetConditions();
  float *h_ux = this->Ux.host<float>();
  float *h_uy = this->Uy.host<float>();
  float *h_uz = this->Uz.host<float>();
  float *h_rho = this->rhos.host<float>();
  

  
  int index=0;
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


  // Ux 
  /*
  index = 0 ;
  file_vtk << "SCALARS density float \n";
  file_vtk << "LOOKUP_TABLE default \n";
  for( int k = 0 ; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++)
      for( int i = 0 ; i < Lx ; i++)
	{ // h_nodes = 0 si es fluido , 1 si es 
	  file_vtk << h_rho[index] << "\n";
	  index++;
	}
  */
  /// Ux field 
  file_vtk << "VECTORS VelocityField float \n";
  index = 0;
  for( int k = 0 ; k < Lz ; k++)
    for( int j = 0 ; j < Ly ; j++)
      for( int i = 0 ; i < Lx ; i++)
      {
	file_vtk << h_ux[ index ] << " " << h_uy[ index ] << " " << h_uz[ index ] << " \n";
	index++;
      }
  
  file_vtk.close();
  delete [] h_ux;
  delete [] h_uy;
  delete [] h_uz;
  delete [] h_rho;
}
void LatticeBolztmannD3Q15::runTillStacionary(float eps){
  
  int i = 0;
  bool stacionary = false;
  this->BuildWalls();
  //this->BounceBackBoundaries();
  while( !stacionary ){

    this->Adveccion();
    this->BounceBackBoundaries();
    this->MacroscopicQuantities();
    this->Colission();
    i++;
    std::cout << "# step:" << i << std::endl;
    if ( i >= 100000 ) break;
  }
  // una vez terminada guardar
  this->saveVTK("./lid-driven-bg.vtk" );
}
void LatticeBolztmannD3Q15::runBounceBack(int steps){
  this->BuildWalls();
  for( int i =  0 ; i< steps; i++){
    this->Adveccion();
    this->BounceBackBoundaries();
    this->MacroscopicQuantities();
    this->Colission();
    std::cout << "# step:" << i << std::endl;
  }
  this->saveVTK("./plane.vtk");
}
void  LatticeBolztmannD3Q15::BounceBackBoundaries(){
  
  this->InterchangePopulations(1);
  this->InterchangePopulations(2);
  this->InterchangePopulations(3);
  this->InterchangePopulations(7);
  this->InterchangePopulations(8);
  this->InterchangePopulations(9);
  this->InterchangePopulations(10);

}
void  LatticeBolztmannD3Q15::BuildWalls(){
  // un arreglo de paredes... 
 
  // lid driven  
   this->ConditionsU(0,span,span) = 1; // plano x = 0
  this->ConditionsU(1,span,span ) = 1 ; // x = 1
  
  this->ConditionsU(Lx-1,span,span) = 1; // plano x = lx
  this->ConditionsU(Lx-2,span,span) = 1 ; // plano x =lx-1

  this->ConditionsU(span, span ,  0 ) = 1 ; // plano z = 0
  this->ConditionsU( span, span, 1 ) = 1; //z = 1

  this->ConditionsU(span, span , Lz - 1 ) = 1 ; // plano z = Lz
  this->ConditionsU(span, span, Lz-2)= 1; // plano z = Lz-1

  this->ConditionsU(span , 0 , span  ) = 2; // plano y = 0 ;
  this->ConditionsU(span, 1 , span ) = 2;
  
  this->ConditionsU(span , Ly-1 , span  ) = 1; // plano y = 0 ;
  this->ConditionsU(span, Ly-2 , span ) = 1;  //plano y = 1 
  
  
 
  // ventiladores en la punta
  /*
  this->ConditionsU( 0 , span , span ) = 2; // velocidiad algo
  this->ConditionsU( 0 , span, span) = 2;
  this->ConditionsU( 1 , span , span ) = 2; // velocidiad algo
  this->ConditionsU( 1 , span, span) = 2;
  */
}
int  LatticeBolztmannD3Q15::GetindexOpositeVector(int alpha){

  for ( int i = 1 ; i < q ; i++){
    
    if ( V[0][i] == -1*V[0][alpha] && V[1][i] == -1*V[1][alpha] && V[2][i] == -1*V[2][alpha] ) { 
      return i;
    }
  } 
}
void  LatticeBolztmannD3Q15::InterchangePopulations(int alpha ) {
  
  int alpha_ = this->GetindexOpositeVector(alpha);
  array buff = this->fs[alpha];
  buff = this->fs[alpha];
  // wall still conditions
  fs[alpha]( ConditionsU == 1 ) = fs[alpha_](ConditionsU == 1 ); // fs[1] = fs[7]
  fs[alpha_](ConditionsU == 1 ) = buff(ConditionsU == 1 ); // fs[7] = fs[1]

  // moving wall conditions
  fs[alpha]( ConditionsU == 2 )= fs[alpha_](ConditionsU==2) + 6*W[alpha]*rhos(ConditionsU==2)*V[0][alpha]*(UMAX);
  fs[alpha_](ConditionsU ==2) = buff(ConditionsU==2) + 6*W[alpha_]*rhos(ConditionsU==2)*V[0][alpha_]*(UMAX);
  
}
void LatticeBolztmannD3Q15::recalculateViscosity(float Re,float r){
  this->UMAX = 0.05;
 
  float ulb = UMAX;
  float l = r; // longitud caracteristica
  float nulb = ulb*l/Re; // viscosidad.
  this->beta = 1.0/(3*nulb + 0.5); 
  std::cout << "# tau:" << 1/beta << std::endl;
  std::cout << "# viscocidad:"  << nulb << std::endl;
  std::cout << "# UMAX : " << ulb << std::endl;
  std::cout << nulb << " " ;
  
}
void LatticeBolztmannD3Q15::runAndSaveFs(int steps, std::string rute){
  
  this->BuildWalls();
  af::timer timeit = timer::start();
  for( int ii = 0 ; ii < steps ; ii++){
    this->Adveccion();
    this->BounceBackBoundaries();
    this->MacroscopicQuantities();
    this->Colission();
    std::cout <<  "# step:" << ii << std::endl;  
  }
  double t = af::timer::stop(timeit);
  //std::cout << "#time(s):" << t << std::endl;
  // ahora calcular fuerzas 
  this->BounceBackBoundaries();
 
  this->SaveFs(rute); // guarda las ffunciones fi 
  //this->saveVTK("./wtf-fuerza.vtk" );
  
}
