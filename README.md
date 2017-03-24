

Code i used for my thesis. it performs lattice bolztmann simulations on 3D.
Next a brief description of the files.

the inc/ folder holds the headers of the classes used.

inc/ELB3DQ15.h  - an incomplete implementation of an entropic lattice bolztmann.

inc/ForceModule.h a class used to compute hydrodynamic forces over a surface, using the momentum exchange method. for it to work a folder including the vectors of the lattice-bolztmann and the points of the surface, must be included.


inc/LBD3Q15.h this class performs the 3D simulation. To be created requieres to specify the size of simulation domain ( Lx, Ly , Lz ). it cointains a method to save a vtk file, than can be readed by ParaView.

inc/ObstacleReader.h  this class, takes a string pointing to a file and  turns dat file into a ArrayFire array object. this array will represent solid nodes in the simulation domain. The whole point of the thesis was to find a way to da this, to take a complex 3D model that can be designed in a program like Blender Autodesk , take this geometric information and introduce it in the simulation grid. This was correctly done... bit it still has many issues, the work flow i ended up using was extremly complicated and a single mistake lead to wastes of time and computation time. but in the end the goal of "Load 3d objects into LB simulations" and run those computations on GPUs, was achievied.

