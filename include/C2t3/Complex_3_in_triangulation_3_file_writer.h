
#ifndef CGAL_IO_COMPLEX_3_IN_TRIANGULATION_3_FILE_WRITER_H
#define CGAL_IO_COMPLEX_3_IN_TRIANGULATION_3_FILE_WRITER_H

#include <CGAL/IO/File_writer_OFF.h>

namespace CGAL { 
	namespace Mesh_3 {
	
	template < class Tr>
	typename Tr::size_type number_of_facets_on_surface(const Tr& T);

	template <class Triangulation>
	class Write_to_OFF_file 
	{
		CGAL::File_writer_OFF off;
		std::ostream& os;

		typedef Triangulation Tr;
		public:
		Write_to_OFF_file(std::ostream& os,  bool verbose) : off(verbose), os(os) {}

		bool write_header(const typename Tr::size_type number_of_vertices,
						  const typename Tr::size_type number_of_facets)
		{
			off.header().set_no_comments(true);
			off.write_header( os,
							  number_of_vertices,
							  0, // fake number of halfedges, not used.
							  number_of_facets);
			return os.good();
		}

		bool write_vertex(const typename Tr::Vertex_handle& v)
		{
			const typename Tr::Point& p = v->point();
			off.write_vertex(p.x(), p.y(), p.z());
			return os.good();
		}

		bool begin_facets()
		{
			off.write_facet_header();
			return os.good();
		}

		bool write_facet(const int index1,
					     const int index2,
						 const int index3)
		{
			off.write_facet_begin(3);
			off.write_facet_vertex_index(index1);
			off.write_facet_vertex_index(index2);
			off.write_facet_vertex_index(index3);
			off.write_facet_end();
			return os.good();
		}

		bool write_footer()
		{
			off.write_footer();
			return os.good();
		}

	}; // end class Write_to_OFF_file

	enum IO_option {	NO_OPTION = 0,
						IO_ORIENT_SURFACE = 1,
						IO_VERBOSE = 2 };
	
	} // End namespace Mesh_3


// Writes in the OFF file ALL the surface facets, regardless of their subdomain
template <class C3t3>
bool output_subdomain_surface_facets_to_off(std::ostream& os,
				   const C3t3& c3t3,
				   typename C3t3::Subdomain_index subdomain,
				   typename C3t3::Subdomain_index surface,
				   int options = 
				   Mesh_3::IO_ORIENT_SURFACE)
{
	using CGAL::Mesh_3::number_of_facets_on_surface;

	typedef typename C3t3::Triangulation Tr ;  
	typedef typename C3t3::Subdomain_index Subdomain_index ;
	typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
	typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
	typedef typename C3t3::Facets_in_complex_iterator Facets_in_complex_iterator ;
	typedef typename C3t3::Cells_in_complex_iterator Cells_in_complex_iterator ;
	typedef typename Tr::Facet_circulator Facet_circulator;
	typedef typename Tr::Facet Facet;
	typedef typename Tr::Edge Edge;
	typedef typename Tr::Cell_handle Cell_handle;
	typedef typename Tr::Vertex_handle Vertex_handle;
	typedef typename Tr::Point Point;
	typedef typename Tr::Geom_traits Gt;

	// Header.
	const Tr& tr = c3t3.triangulation();

	bool success = true;

	Mesh_3::Write_to_OFF_file<Tr> off( os, (options & Mesh_3::IO_VERBOSE) != 0 ) ;

	success &= off.write_header( tr.number_of_vertices(), c3t3.number_of_facets() ) ;
  
	CGAL_assertion(c3t3.number_of_facets() == number_of_facets_on_surface(tr));
  
	// Finite vertices coordinates.
	std::map<Vertex_handle, int> V;
	int inum = 0;
	for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
		 vit != tr.finite_vertices_end();
		 ++vit )
	{
	    V[vit] = inum++;
		success &= off.write_vertex(vit);
	}

	success &= off.begin_facets();

	if((options & Mesh_3::IO_ORIENT_SURFACE) == 0) {

		for( Finite_facets_iterator fit = tr.finite_facets_begin();
			fit != tr.finite_facets_end(); ++fit ) {

			const Cell_handle cell = fit->first;
			const int& index = fit->second;
			if ( cell->is_facet_on_surface(index)==true && 
				 c3t3.subdomain_index( fit->first ) == subdomain &&
				 ( ( c3t3.surface_index( *fit ).first == surface ) || ( c3t3.surface_index( *fit ).second == surface ) ) ) 
			{
				const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
				const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
				const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
				success &= off.write_facet(index1, index2, index3);
			}
		}
	}
	else // if facets must be oriented
	{		
		Finite_facets_iterator fit = tr.finite_facets_begin();
		std::set<Facet> oriented_set;
		std::stack<Facet> stack;

		typename Tr::size_type number_of_facets = c3t3.number_of_facets();

		CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

		while (oriented_set.size() != number_of_facets) 
		{
			while ( fit->first->is_facet_on_surface(fit->second) == false ||
					oriented_set.find(*fit) != oriented_set.end() ||
					// oriented_set.find(c3t3.opposite_facet(*fit)) !=	oriented_set.end() ) 
					// oriented_set.find(tr.mirror_facet(*fit)) !=	oriented_set.end() ) 
					oriented_set.find( std::make_pair(fit->first->neighbor(fit->second), fit->first->neighbor(fit->second)->index(fit->first) ) ) != oriented_set.end() ) 
			{
				++fit;
			}
			oriented_set.insert(*fit);
			stack.push(*fit);
			while(! stack.empty() )
			{
				Facet f = stack.top();
				stack.pop();
				for(int ih = 0 ; ih < 3 ; ++ih) {
					const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
					const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
					
					// Face status does not exist on a C3t3... so we will compute it
					// Is regular? (Count the number of incident facets)					
					Facet fn;
					Facet_circulator fc = tr.incident_facets( Edge( f.first, i1, i2 ) ) ;
					Facet_circulator fc_done = fc ;
					int numValidIncidentFacets = 0 ;
					do {
						// Check if the current facet is on the domain/surface index
						if ( fc->first->is_facet_on_surface(fc->second)==true && 
							 c3t3.subdomain_index( fc->first ) == subdomain &&
							 ( ( c3t3.surface_index( *fc ).first == surface ) || ( c3t3.surface_index( *fc ).second == surface ) ) ) {
							numValidIncidentFacets++ ;
							if ( *fc != f ) {
								// std::cout << "Neighbor!" << std::endl ;
								// Recover also the neighbor (used later)
								// fn = *fc ;
								fn = std::make_pair(fc->first->neighbor(fc->second), fc->first->neighbor(fc->second)->index(fc->first) ) ;
							}
						}
					} while( ++fc != fc_done ) ;
					
					// std::cout << "numValidIncidentFacets = " << numValidIncidentFacets << std::endl ;

					if ( numValidIncidentFacets == 2 ) { // i.e., REGULAR
						
						if (oriented_set.find(fn) == oriented_set.end()) {
							// if(oriented_set.find(tr.mirror_facet(fn)) == oriented_set.end())
							if( oriented_set.find( std::make_pair(fn.first->neighbor(fn.second), fn.first->neighbor(fn.second)->index(fn.first) ) ) == oriented_set.end() )
							{								
								oriented_set.insert(fn);
								stack.push(fn);
							}
							else {
								success = false; // non-orientable
							}
						}
					}
					else if( numValidIncidentFacets != 1 ) { // i.e., is not BOUNDARY
						success = false; // non manifold, thus non-orientable
					}
				} // end "for each neighbor of f"
			} // end "stack non empty"
		} // end "oriented_set not full"
    
		for( typename std::set<Facet>::const_iterator fit = 
			 oriented_set.begin();
			 fit != oriented_set.end();
			 ++fit)
		{
			const typename Tr::Cell_handle cell = fit->first;
			const int& index = fit->second;
			const int index1 = V[cell->vertex(tr.vertex_triple_index(index, 0))];
			const int index2 = V[cell->vertex(tr.vertex_triple_index(index, 1))];
			const int index3 = V[cell->vertex(tr.vertex_triple_index(index, 2))];
			success &= off.write_facet(index1, index2, index3);
			CGAL_assertion_code(++nb_facets);
		}

		CGAL_assertion(nb_facets == number_of_facets);
	} // end if(facets must be oriented)

	success &= off.write_footer();
	return success;
}

} // End namespace CGAL







#endif // CGAL_IO_COMPLEX_3_IN_TRIANGULATION_3_FILE_WRITER_H