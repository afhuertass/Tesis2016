
#pragma once

#define BOOST_PARAMETER_MAX_ARITY 12
// #define CGAL_MESH_3_VERBOSE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Triangulation_3.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>


#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;

typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

// definiciones mesh.cpp 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_m;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron_m, K> Mesh_domain;
typedef K::Point_3 Point_3_k;
typedef CGAL::Triangulation_3<K> Triangulation_3;

//typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Mesh_triangulation_3<
  Mesh_domain,
  CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
  CGAL::Parallel_tag                        // Tag to activate parallelism
  >::type Tr;
#else
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef C3t3::Cell_handle Cell_handle;
typedef Tr::Locate_type Locate_type;
typedef Triangulation_3::Point Point;


