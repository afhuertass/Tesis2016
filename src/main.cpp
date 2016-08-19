#include <LBD3Q15.h>
#include <ForceModule.h>
#include <ObstacleLoader.h>
#include <ELB3DQ15.h>
void calculateForce(){
  
  std::string ruta = "./simulation-runs/objectData/cube";
  ForceModule fm = ForceModule(ruta);

  fm.calculateForce();
}
void runSpheres(){
  int Lx, Ly , Lz, steps = 15000;
  Lx = 150 , Ly = 75 , Lz = 75;
  LatticeBolztmannD3Q15 fluid(Lx, Ly , Lz);
  fluid.Inicie( 1 , 0 , 0 , 0);
  ObstacleLoader obj( Lx , Ly , Lz );
  std::string esfera = "./simulation-runs/objectData/sphere/";
  
  fluid.setArrayC(  obj.loadFile( esfera + "/ws.txt", Lx , Ly , Lz) );
  
  // calcular la simulacion , calculas las fuerzas 
  //ForceModule fm = ForceModule(esfera);
  float Re = 100;
  for( Re = 4 ; Re <= 100 ; Re += 2){
    std::cout << "### Reynolds number " << Re << std::endl;
    fluid.Inicie( 1 , 0 , 0 , 0);
    fluid.recalculateViscosity( Re , 8.75);
    fluid.runAndSaveFs( steps , esfera);
    fluid.saveVTK("cubo-re100.vtk");
    ForceModule fm = ForceModule(esfera);
    //std::cout << Re  << " " ;
    fm.calculateForce();
    std::cout << std::endl;
      
    }
 
}
void runCubes() {
 int Lx, Ly , Lz, steps = 15000;
  Lx = 150 , Ly = 75 , Lz = 75;
  LatticeBolztmannD3Q15 fluid(Lx, Ly , Lz);
  fluid.Inicie( 1 , 0 , 0 , 0);
  ObstacleLoader obj( Lx , Ly , Lz );
  std::string esfera = "./simulation-runs/objectData/cube";
  
  fluid.setArrayC(  obj.loadFile( esfera + "/ws.txt", Lx , Ly , Lz) );
  
  // calcular la simulacion , calculas las fuerzas 
  //ForceModule fm = ForceModule(esfera);
  float Re = 50;
  //for( Re = 5 ; Re <= 50 ; Re += 5){
    std::cout << "### Reynolds number " << Re << std::endl;
    fluid.Inicie( 1 , 0 , 0 , 0);
    fluid.recalculateViscosity( Re , 12);
    std::cout << Re  << " " ;
    fluid.runAndSaveFs( steps , esfera  );
    ForceModule fm = ForceModule(esfera);
    //std::cout << Re  << " " ;
    fm.calculateForce();
    std::cout << std::endl;
      
    //}
  fluid.saveVTK("cube-re50.vtk");
}
void runPlane(){
  int Lx, Ly , Lz ;
  Lx = 300 ; Ly = 150 ; Lz = 150; 
  int steps = 15000;
  LatticeBolztmannD3Q15 fluid(Lx, Ly , Lz);
  fluid.Inicie( 1 , 0.0 , 0 , 0   );
  
  ObstacleLoader obl( Lx , Ly , Lz);
  
  std::string ruta = "./simulation-runs/objectData/plane/ws.txt";
  
  fluid.setArrayC ( obl.loadFile(ruta, Lx, Ly , Lz));
  
  fluid.runBounceBack( steps ); 
  
} 
void runLidDriven(){
  int Lx, Ly, Lz;
  Lx = 300;
  Ly = 120;
  Lz = 60;
  LatticeBolztmannD3Q15 lid(Lx, Ly , Lz);
  lid.Inicie(1 , 0.05 , 0 ,0 );
  ObstacleLoader obj( Lx , Ly , Lz );
  std::string ruta = "./simulation-runs/objectData/background/ws.txt";
  lid.setArrayC ( obj.loadFile(ruta, Lx, Ly , Lz));
  
  af::timer timeit = af::timer::start();
  lid.recalculateViscosity( 100 , Lx);
  lid.runTillStacionary( 1e-7 );
  double t = af::timer::stop( timeit );
  std::cout << "#Time (s): lid driven " << t << std::endl;
  
}
void spheresData() {
  
  
 
}
void EntropicCall(){
  int pasos = 1; 
  int Lx, Ly, Lz;
  Lx = Ly = Lz = 5;
  ELatticeBolztmannD3Q15 Efluid( Lx , Ly , Lz );
  Efluid.Inicie( 1 , 0.0 , 0 , 0 );
  af::timer timeit = timer::start();
  for ( int i = 0 ; i < pasos ; i++){
    
    Efluid.Adveccion();
    Efluid.Colission();
  }
  double t = af::timer::stop(timeit);
  double MUPS = (Lx*Ly*Lz)/1e6;
  MUPS = (MUPS*pasos)/t;
  std::cout << "#MUPs:" << MUPS << std::endl;
  
  //Efluid.saveVTK("./entropic-lid-driven.vtk");
} 


int main(int argc , char * argv[] ){
  
  int device = argc > 1 ? atoi(argv[1]) : 0 ;
  af::setDevice(device);

  //EntropicCall();
  
  //runKotaiiki();
  runLidDriven();
  //experiments();
  //runSpheres();
  // runPlane();
  //calculateForce();
}
