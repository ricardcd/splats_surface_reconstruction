// Functions for post processing of Complex_2_in_triangulation_3 structures, after being used for Surface Meshing

#ifndef COMPLEX_2_IN_TRIANGULATION_3_POST_PROCESSING_H
#define COMPLEX_2_IN_TRIANGULATION_3_POST_PROCESSING_H


namespace CGAL {

namespace C2t3_PostProcessing { 

// Remove triangles forming a sliver tetrahedron, returns the number of facets removed
template < class C2t3 >
int removeSlivers( C2t3& c2t3 ) {

	typedef typename C2t3::Triangulation					Tr ;
	typedef typename Tr::Triangulation_data_structure		Tds ;
	typedef typename Tr::Geom_traits::Triangle_3			Triangle_3 ;
	typedef typename C2t3::Vertex_iterator					Vertex_iterator ;
	typedef typename C2t3::Edge_iterator					Edge_iterator ;
	typedef typename Tr::Cell_handle						Cell_handle ;
	typedef typename Tr::Facet								Facet ;
	typedef typename Tr::Cell_circulator					Cell_circulator ;
	typedef typename Tr::Facet_circulator					Facet_circulator ;


	// Run over all facets of the triangulation
	Tr tr = c2t3.triangulation() ;
	Tds tds = tr.tds() ;

	Edge_iterator eit ;
	int numRemoved = 0 ;
	for ( eit = c2t3.edges_begin(); eit != c2t3.edges_end(); ++eit ) {
		if ( c2t3.face_status( *eit ) == C2t3::SINGULAR ) {
			Cell_circulator cci = tds.incident_cells( *eit ) ;			
			Cell_circulator cc = ++cci ;			
			bool foundTetra = false ;			
			while ( !foundTetra && ( ++cc != cci ) ) {				
				// All 4 facets are in the 2D complex
				foundTetra = c2t3.is_in_complex( cc, 0 ) && 
							 c2t3.is_in_complex( cc, 1 ) && 
							 c2t3.is_in_complex( cc, 2 ) && 
							 c2t3.is_in_complex( cc, 3 ) ;
			}
			
			if ( foundTetra ) {
				bool removedSome = false ;
				int i = 0 ;
				int j = 1 ;
				while ( !removedSome && i < 3 ) {
					//std::cout << "  i = " << i << " | j = " << j << std::endl ;
					// Remove one of the regular edges
					if ( c2t3.face_status( cc, i, j ) == C2t3::REGULAR ) {						
						// Get the two facets associated to this edge
						//c2t3.remove_from_complex( cc, (i+2)%4 ) ;
						//c2t3.remove_from_complex( cc, (j+2)%4 ) ;

						// Get the proper indexes representing the 2 faces associated to the edge
						int k1, k2 ;
						if ( ((i+1)%4) != j ) 
							k1 = (i+1)%4 ;
						else
							k1 = (i-1)%4 ;
						if ( ((j+1)%4) != i ) 
							k2 = (j+1)%4 ;
						else
							k2 = (j-1)%4 ;

						// Compare the angle between normals of the triangles (should be compatibles)
						Triangle_3 t1 = tr.triangle( cc, (k1+2)%4 ) ;
						Triangle_3 t2 = tr.triangle( cc, (k2+2)%4 ) ;
						if ( t1.supporting_plane().orthogonal_vector() * t2.supporting_plane().orthogonal_vector() >= 0 ) {

							// Remove these faces from the 2D complex
							c2t3.remove_from_complex( cc, (k1+2)%4 ) ;
							c2t3.remove_from_complex( cc, (k2+2)%4 ) ;

							removedSome = true ;
							numRemoved += 2 ;

						}
					}
					j++ ;
					if ( j > 3 ) {
						i++ ;
						j = i + 1 ;
					}
				}
			}
		}
	}

	return numRemoved ;

}



// Removes ears in reconstruction, returns the number of facets deleted
template < class C2t3 >
int removeEars( C2t3& c2t3 ) {

	typedef typename C2t3::Triangulation					Tr ;
	typedef typename Tr::Triangulation_data_structure		Tds ;
	typedef typename Tr::Geom_traits::Triangle_3				Triangle_3 ;
	typedef typename C2t3::Vertex_iterator							Vertex_iterator ;
	typedef typename C2t3::Edge_iterator								Edge_iterator ;
	typedef typename Tr::Cell_handle									Cell_handle ;
	typedef typename Tr::Facet										Facet ;
	typedef typename Tr::Cell_circulator								Cell_circulator ;
	typedef typename Tr::Facet_circulator							Facet_circulator ;

	// Run over all facets of the triangulation
	Tr tr = c2t3.triangulation() ;
	Tds tds = tr.tds() ;

	Edge_iterator eit ;
	int numRemoved = 0 ;
	for ( eit = c2t3.edges_begin(); eit != c2t3.edges_end(); ++eit ) {
		if ( c2t3.face_status( *eit ) == C2t3::SINGULAR ) {
			// Non manifold edge found
			
			Facet_circulator fc0 = tds.incident_facets( *eit ) ;
			Facet_circulator fc = ++fc0 ;
			while ( ++fc != fc0 ) {
				if ( c2t3.face_status( fc->first, fc->second ) != C2t3::NOT_IN_COMPLEX ) {
					// Remove a facet if its 3 edges are not shared by another triangle
					bool invalidFacet = c2t3.face_status( fc->first, (fc->second+1)%4 ) != C2t3::REGULAR ||
										c2t3.face_status( fc->first, (fc->second+2)%4 ) != C2t3::REGULAR ||
										c2t3.face_status( fc->first, (fc->second+3)%4 ) != C2t3::REGULAR ;
					if ( invalidFacet ) {
						//c2t3.remove_from_complex( fc->first, fc->second ) ;
						c2t3.remove_from_complex( *fc ) ;
						numRemoved++ ;
					}
				}
			}		
			
		}
	}

	return numRemoved ;

}



// Removes all facets indicent to a non-manifold vertex, returns the number of deleted vertices
template < class C2t3 >
int removeNonManifoldVertices( C2t3& c2t3 ) {

	typedef typename C2t3::Triangulation					Tr ;
	typedef typename Tr::Triangulation_data_structure		Tds ;
	typedef typename Tr::Geom_traits::Triangle_3				Triangle_3 ;
	typedef typename C2t3::Vertex_iterator							Vertex_iterator ;
	typedef typename C2t3::Edge_iterator								Edge_iterator ;
	typedef typename Tr::Cell_handle									Cell_handle ;
	typedef typename Tr::Facet										Facet ;
	typedef typename Tr::Cell_circulator								Cell_circulator ;
	typedef typename Tr::Facet_circulator							Facet_circulator ;
	
	Tr tr = c2t3.triangulation() ;
	Tds tds = tr.tds() ;

	Vertex_iterator vit ;
	int numRemoved = 0 ;
	for ( vit = c2t3.vertices_begin(); vit != c2t3.vertices_end(); ++vit ) {
		if ( c2t3.face_status( vit ) == C2t3::SINGULAR ) {
			
			typename std::vector< Facet > incidentFacets ;
			typename std::vector< Facet >::iterator fit ;
			tds.incident_facets( vit, std::back_inserter( incidentFacets ) ) ;			
			for ( fit = incidentFacets.begin(); fit != incidentFacets.end(); ++fit ) {
				if ( c2t3.face_status( *fit ) != C2t3::NOT_IN_COMPLEX ) {
					c2t3.remove_from_complex( *fit ) ;
					numRemoved++ ;
				}
			}
		}
	}

	return numRemoved ;
}



// Returns true if embedded surface is manifold, false otherwise
template < class C2t3 >
bool manifoldnessCheck( C2t3& c2t3, int& numNonManifoldEdges, int& numNonManifoldVertices ) {
	
	typedef typename C2t3::Vertex_iterator							Vertex_iterator ;
	typedef typename C2t3::Edge_iterator							Edge_iterator ;	

	// Count non manifold edges
	numNonManifoldEdges = 0 ;
	Edge_iterator eit ;
	for ( eit = c2t3.edges_begin(); eit != c2t3.edges_end(); ++eit ) {
		if ( c2t3.face_status( *eit ) == C2t3::SINGULAR ) {
			numNonManifoldEdges++ ;
		}
	}
	
	numNonManifoldVertices = 0 ;
	Vertex_iterator vit ;
	for ( vit = c2t3.vertices_begin(); vit != c2t3.vertices_end(); ++vit ) {
		if ( c2t3.face_status( vit ) == C2t3::SINGULAR ) {
			numNonManifoldVertices++ ;
		}
	}

	return ( numNonManifoldEdges == 0 && numNonManifoldVertices == 0 ) ;

}

} // End namespace C2t3_PostProcessing

} // End namespace CGAL


#endif // COMPLEX_2_IN_TRIANGULATION_3_POST_PROCESSING_H