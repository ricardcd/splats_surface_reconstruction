// Ricard Campos: Modification of CGAL s example Jet_fitting_3/Single_estimation.cpp to evaluate the computed function at any point

#include <CGAL/Cartesian.h>
#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
#include <CGAL/Random.h>

typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     Point_3;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form     My_Monge_form;

#include <fstream>
#include <vector>



int main(int argc, char *argv[])
{
  size_t d_fitting = 4;
  size_t d_monge = 4;
  const char* name_file_in = "../data/x2y2.pts";
  const char* name_file_out = "out.txt" ;
  int numPts = 1000 ;
  //check command line
  if (argc<5)
    {
      std::cout << " Usage : " << argv[0]
                << " <inputPoints.txt> <outputPoints.txt> <d_fitting> <d_monge> <num_points>" << std::endl
		<< "test with default arguments" << std::endl;
    }
  else {
    name_file_in = argv[1] ;
	name_file_out = argv[2] ;
    d_fitting = std::atoi( argv[3] ) ;
    d_monge = std::atoi( argv[4] ) ;
	numPts = std::atoi( argv[5] ) ;
  }

  //open the input file
  std::ifstream inFile(name_file_in);
  if ( !inFile )
    {
      std::cerr << "cannot open file for input\n";
      exit(-1);
    }
  
  //initalize the in_points container
  double x, y, z;
  std::vector<Point_3> in_points;
  while (inFile >> x) {
    inFile >> y >> z;
    Point_3 p(x,y,z);
    in_points.push_back(p);
  }
  inFile.close();
  std::cout << "Number of input points: " << in_points.size() << std::endl ;

  // fct parameters

  My_Monge_form monge_form;
  My_Monge_via_jet_fitting monge_fit;
  monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

  
   
  // Throw some random points and compute their projection on the monge
  std::vector< Point_3 > proj ;  
  double x_inf, x_sup, y_inf, y_sup, z_inf, z_sup ;
  //x_inf = -1.0419 ; 
  //x_sup = 1.2883 ;
  //y_inf = -9.4806 ;
  //y_sup = -7.0315 ;
  //z_inf = 6.3708 ;
  //z_sup = 8.8629 ;

  /*x_inf = -1.0 ; 
  x_sup = 1.0 ;
  y_inf = -1.0 ;
  y_sup = 1.0 ;
  z_inf = -1.0 ;
  z_sup = 1.0 ;*/
  
  x_inf = -4.2327 ; 
  y_inf = -1.1236 ;
  z_inf = 9.6333 ;
  x_sup = -1.9146 ;
  y_sup =  1.1236 ;
  z_sup = 11.8744 ;
  
  for ( int i = 0; i < numPts; i++ ) {
	double x = CGAL::default_random.get_double( x_inf, x_sup ) ;
	double y = CGAL::default_random.get_double( y_inf, y_sup ) ;
	double z = CGAL::default_random.get_double( z_inf, z_sup ) ;

	Point_3 curProj = monge_form.project( Point_3( x, y, z ) ) ;
	proj.push_back( curProj ) ;
  }

  // Debug: Project also ( 0,0,0 ) point
  //Point_3 curProj = monge_form.project( Point_3( 0.0, 0.0, 0.0 ) ) ;
  //proj.push_back( curProj ) ;


  //OUTPUT on std::cout
  CGAL::set_pretty_mode(std::cout);
  std::cout << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl
	    //<< monge_form;
		<< "origin = [ " << monge_form.origin().x() << " " << monge_form.origin().y() << " " << monge_form.origin().z() << " ] ;" << std::endl 
		<< "vecX = [ " << monge_form.maximal_principal_direction().x() << " " << monge_form.maximal_principal_direction().y() << " " << monge_form.maximal_principal_direction().z() << " ] ;" << std::endl 
		<< "vecY = [ " << monge_form.minimal_principal_direction().x() << " " << monge_form.minimal_principal_direction().y() << " " << monge_form.minimal_principal_direction().z() << " ] ;" << std::endl 
		<< "vecZ= [ " << monge_form.normal_direction().x() << " " << monge_form.normal_direction().y() << " " << monge_form.normal_direction().z() << " ] ;" << std::endl ;
  std::cout << monge_form << std::endl ;
  std::cout  << "condition_number : " << monge_fit.condition_number() << std::endl
	     << "pca_eigen_vals and associated pca_eigen_vecs :"  << std::endl;
  for (int i=0; i<3; i++)
    std::cout << monge_fit.pca_basis(i).first << std::endl
	      << monge_fit.pca_basis(i).second  << std::endl;


  // Print points on txt file
  std::ofstream out ; 
  out.open( name_file_out ) ;
  for ( int k = 0; k < proj.size(); k++ ) {
	  out << proj[k].x() << " " << proj[k].y() << " " << proj[k].z() << std::endl ;
  }
  out.close() ;

  return 0;
}
