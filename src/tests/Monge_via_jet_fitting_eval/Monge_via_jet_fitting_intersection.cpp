// Ricard Campos: Modification of CGAL s example Jet_fitting_3/Single_estimation.cpp to evaluate the computed function at any point

#include <CGAL/Cartesian.h>
#include "Monge_via_jet_fitting_eval/Monge_via_jet_fitting_eval.h"
#include <CGAL/Random.h>

typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     Point_3;
typedef Data_Kernel::Vector_3    Vector_3;
typedef Data_Kernel::Segment_3   Segment_3;
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
  std::cout << "Monge fitted " << std::endl ;

  // Create a segment
  double residuals[1] ;
  Point_3 start( -1, -1, 0.6 ) ;
  Point_3 end( 1, 1, 0.6 ) ;
  Segment_3 seg( start, end ) ;
  Point_3 itPoint ;


  double lb, ub ;
  lb = 0 ;
  ub = CGAL::sqrt( seg.squared_length() ) ;

  std::cout << "Lower bound = " << lb << std::endl ;
  std::cout << "Upper bound = " << ub << std::endl ;

  std::cout << "Computing intersection with a segment" << std::endl ;
  bool doIntersect = monge_form.intersection( start, seg.to_vector(), itPoint, residuals, lb, ub ) ;  


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

  std::cout << std::endl << std::endl << "IntersectionPoint = [ " << itPoint.x() << " " << itPoint.y() << " " << itPoint.z() << "] " << std::endl ;
  //std::cout << "Residuals = [ " << residuals[0] << " " << residuals[1] << " " << residuals[2] << "] " << std::endl ;
  std::cout << "Residuals = [ " << residuals[0] << "] " << std::endl ;
  

  return 0;
}
