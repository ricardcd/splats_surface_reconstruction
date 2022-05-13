class Construct_intersection { 
				
		typedef boost::optional<typename Tree::Object_and_primitive_id> AABB_intersection;
		const Self& self ;
		
	public: 
		
		typedef typename Tree::Object_and_primitive_id	Obj_prim_id ;
		
		Construct_intersection( const Self& self ) : self( self ) {}

		// Intersection between the surface and a segment
		Intersection operator()( const Segment_3& segment ) const 
		{
			// Return the already computed intersection
			// (It should have been computed in Do_intersect_surface)
			// std::cout << "Recovering the computed intersection..." << std::endl ;
			return self.m_computedIntersection ;
			// return computeIntersection( segment ) ;
		}

		Intersection computeIntersection( const Segment_3& segment ) const {
			MesherIntersectionCounterSegments++ ;
			MesherIntersectionCounter++ ;

			//std::cout << "Creating the new segment" << std::endl ;
			Segment_3 newSegment = minLengthSegment( segment, self.minLengthQuerySegment() ) ;
			//std::cout << "New segment created: " << newSegment << std::endl ;
			//std::cout << "Old segment        : " << segment << std::endl ;
			//if( newSegment == segment )
				//std::cout << "They are equal        : " << segment << std::endl ;

			// Check if there is any intersection
			if ( self.tree()->do_intersect( newSegment ) ) {
				//std::cout << "Some intersection" << std::endl ;
				foundMesherIntersectionCounter++ ;
				
				// Get the list of intersections
				std::vector< Obj_prim_id > obj_prim_vect ;
				std::vector< Point_3 >	point_vector ; 
				std::vector< Circle_3 >	circle_vector ; 
				self.tree()->all_intersections( newSegment, std::back_inserter( obj_prim_vect ) ) ;
				typename std::vector< Obj_prim_id >::iterator it ;

				if ( obj_prim_vect.size() == 1 )
					numSegmentsMesherWithOneIntersection++ ;

				// Debug: See what happens if only an intersection is considered if there are more than one primitive interescted...
				//if ( obj_prim_vect.size() > 2 ) { 

				// Get the intersection points and its corresponding Circle_3 primitives
				//for( int i = 0; i < obj_prim_vect.size(); i++ ) {
				for ( it = obj_prim_vect.begin(); it!=obj_prim_vect.end(); ++it ) {
					const Point_3* p = CGAL::object_cast< Point_3 >( &( it->first ) ) ;
					//const Point_3* p = CGAL::object_cast< Point_3 >( &(obj_prim_vect.at(i).first) ) ;
					
					if ( p ) {
						point_vector.push_back( *p ) ;
						Circle_3 c( *(it->second) ) ;
						circle_vector.push_back( c ) ;
					}
				}
				
				// Compute the intersection point				
				Point_3 intPoint = computeIntersectionPoint( self, 
															 newSegment, 
															 point_vector, 
															 circle_vector ) ;

				if ( intPoint.x() == 99999999) { 
					// std::cout << "NULL point! self.holeFillingRadius() = " << self.holeFillingRadius() << std::endl ;

					/* Hole filling, if needed/possible */
					if ( self.holeFillingRadius() > 0 ) {
						Point_3 chp = closeHole( self, newSegment, segment ) ;
						// std::cout << "Hole filling point: " << chp << std::endl ;
						if ( intPoint.x() == 99999999) { 
							return Intersection() ;
						}
						else {
							return Intersection( chp, 
												 self.index_from_surface_patch_index( self.make_surface_index( AABB_primitive_id() ) ),
												 2 ) ;
						}	
					}
					else {
						// std::cout << "NULL Intersection point!: " << std::endl ;
						return Intersection() ;
					}
				}
				//else
				// Resulting intersection point has to be on the original query segment
				else {
					// if( segment.has_on( intPoint ) )
					FT dist = CGAL::sqrt( CGAL::squared_distance( intPoint, newSegment.source() ) ) ;
					FT orSegLength = CGAL::sqrt( segment.squared_length() ) ;
					FT hDifSeg = ( CGAL::sqrt( newSegment.squared_length() ) - orSegLength ) / 2 ;
					if ( segment == newSegment || ( ( dist - hDifSeg > 0 ) && ( dist - hDifSeg < orSegLength ) ) ) {
						// std::cout << "Intersection on the original query segment, returning Intersection object" << std::endl ;
						return Intersection( intPoint, 
											 self.index_from_surface_patch_index( self.make_surface_index( AABB_primitive_id() ) ),
											 2 ) ;
					}
					else {
						// std::cout << "Intersection is NOT on the original query segment, returning Intersection object" << std::endl ;
						return Intersection() ;
					}
				}				
			}
			else {
				if ( self.holeFillingRadius() <= 0 ) {
					// std::cout << "Query segment does not intersect any of the (planar) splats" << std::endl ;
					return Intersection() ;
				}
				else {
					/* Hole filling, if needed/possible */
					Point_3 chp = closeHole( self, newSegment, segment ) ;				
					return Intersection( chp, 
										 self.index_from_surface_patch_index( self.make_surface_index( AABB_primitive_id() ) ),
										 2 ) ;
				}
			}			
		}
		
		// Given the intersected circles, compute a unique intersection point depending on circles' monge and intersection type required
		Point_3 computeIntersectionPoint( const Surface_3& surface, const Segment_3& segment, std::vector< Point_3 >& point_vector, std::vector< Circle_3 >& circle_vector ) const {
			
			// std::cout << "Computing intersection point..." << std::endl ;

			if ( point_vector.size() != circle_vector.size() ) {
				std::cerr << "Error,  point_vector.size() != circle_vector.size()!!!!!" << std::endl ;
			}

			Point_3 ptS = segment.vertex( 0 ) ;
			Point_3 ptE = segment.vertex( 1 ) ;
		
			// Compute the intersection points with the monges (if any)
			std::vector< Point_3 >::iterator pIt = point_vector.begin() ; 
			std::vector< Circle_3 >::iterator cIt = circle_vector.begin() ; 
			std::vector< int > indicesToDelete ;
			int i = 0 ;
			// --- Debug (Start) ---
			//std::cout << "Size = " << point_vector.size() << std::endl ;
			// --- Debug  (End)  ---
			for ( ; pIt != point_vector.end(); ++pIt, ++cIt, i++ ) {
					
				if ( (*cIt).has_monge() ) {

					// --- Debug (Start) ---
					//std::cout << (*cIt).monge() << std::endl ;
					// --- Debug  (End)  ---

					// Compute lower and upper bound
					typename Kernel::Vector_3 vSub1 = ( *pIt - ptS ) ;
					typename Kernel::Vector_3 vSub2 = ( *pIt - ptE ) ;
							
					double lb = -CGAL::sqrt( vSub1.squared_length() ) ;
					double ub = CGAL::sqrt( vSub2.squared_length() ) ;

					//std::cout << "Lower Bound = " << lb << std::endl ;
					//std::cout << "Upper Bound = " << ub << std::endl ;

					// Compute the intersection with the monge, if any
					Point_3 result ;
					double residuals[1] ;	
					// Degree 1, or > than 0 (may need optimization)
					if ( (*cIt).monge().intersection( *pIt, segment.to_vector(), result, residuals, lb, ub, 1E-5 ) ) {						
						// --- Debug (Start) ---
						// std::cout << "Point Before Optim = " << *pIt << std::endl ;
						// --- Debug  (End)  ---
							
						// --- Debug (Start) ---
						// std::cout << "Point After Optim = " << *pIt << std::endl ;
						//std::cout << "Point After Optim = " << point_vector[i] << std::endl ;
						//std::cout << "Intersection Detected, residuals = " << residuals[0] << std::endl ;
						// --- Debug  (End)  ---

						// Check that the new intersection falls inside the predefined radius
						// if ( CGAL::squared_distance( result, cIt->monge().origin() ) < cIt->squared_radius() ) {
							*pIt = result ;							
						/*}
						else {
							std::cout << "Outside radius..." << std::endl ;
							indicesToDelete.push_back( i ) ;
						}*/
					}
					else {
						// --- Debug (Start) ---
						// std::cout << "NO Intersection Detected = " << residuals[0] << std::endl ;
						// --- Debug  (End)  ---
						indicesToDelete.push_back( i ) ;						
					}				

					// --- Debug (Start) ---
					//char temp = cin.get() ;
					// --- Debug  (End)  ---

				}
			}
			
			// --- Debug (Start) ---
			//std::cout << "Number of indices to delete = " << indicesToDelete.size() << std::endl ;
			// std::cout << "number of points = " << point_vector.size() << std::endl ;
			// std::cout << "number of circles = " << circle_vector.size() << std::endl ;
			// --- Debug  (End)  ---

			// Delete intersection points falling outside the segment				
			if ( indicesToDelete.size() == point_vector.size() ) {
				// std::cout << "All candidates were deleted!" << std::endl ;
				return NULL_POINT ;
			}
			if ( indicesToDelete.size() > 0 ) {
				// Sort indices in descending order
				std::sort( indicesToDelete.begin(), indicesToDelete.end() ) ;
				std::reverse( indicesToDelete.begin(), indicesToDelete.end() ) ;
				for ( i = 0; i < indicesToDelete.size(); i++ ) {
					point_vector.erase( point_vector.begin() + indicesToDelete[i] ) ;
					circle_vector.erase( circle_vector.begin() + indicesToDelete[i] ) ;
				}
			}
			
			// Compute the intersection with the surface
			Point_3 intPoint, intPoint2 ;
			double ransacThres ;
			switch ( surface.intersectionType() ) {
				case AABB_Tree_Splats_Mesher_Oracle::Centroid :
					intPoint = IntersectionCentroid( point_vector ) ;
					break ;
				case AABB_Tree_Splats_Mesher_Oracle::WeightedCentroid :
					intPoint = IntersectionWeightedCentroid( segment, 
															 point_vector, 
															 circle_vector, 
															 surface.distanceSigma() ) ;
						
					// --- Debug (Start) ---
					//intPoint2 = IntersectionCentroid( point_vector ) ;
					//std::cout << "Weighted =     " << intPoint << std::endl ;
					//std::cout << "Non Weighted = " << intPoint2 << std::endl ;
					//if ( intPoint == intPoint2 ) 
					//	std::cout << "They are equal!" << std::endl ; 
					//else
					//	std::cout << "They are different" << std::endl ; 
					//std::cout << "Press a key to continue..." << std::endl ; 
					//std::cin.get() ;
					// --- Debug  (End)  ---

					break ;
				case AABB_Tree_Splats_Mesher_Oracle::MixtureOfGaussiansMaximum :
					intPoint = IntersectionMixtureOfGaussiansMaximum(	segment, 
																		point_vector, 
																		circle_vector, 
																		surface.distanceSigma(),
																		surface.mogSigma(),
																		surface.nBins() ) ;						
					// --- Debug (Start) ---
					//std::cout << "Intersection Point = " << intPoint << std::endl ;
					//std::cout << "Press a key to continue..." << std::endl ; 
					//std::cin.get() ;
					// --- Debug  (End)  ---
					break ;					
				case AABB_Tree_Splats_Mesher_Oracle::RansacCentroid :
					// Fixed threshold
					//intPoint = IntersectionRansacCentroid( point_vector, surface.ransacDistThres() ) ;
					// Threshold depending on the segment length
					if ( segment.squared_length() < surface.smallestRansacSegment()*surface.smallestRansacSegment() ) {
						intPoint = IntersectionCentroid( point_vector ) ;
					}
					else {
						ransacThres = segment.squared_length()*surface.ransacDistThres()*surface.ransacDistThres() ; // Squared ransac threshold (ransac compares squared distances)					
						intPoint = IntersectionRansacCentroid( point_vector, ransacThres ) ;
					}
					break ;
				case AABB_Tree_Splats_Mesher_Oracle::RansacWeightedCentroid :					
					if ( segment.squared_length() < surface.smallestRansacSegment()*surface.smallestRansacSegment() ) {
						//std::cout << "Very small SEGMENT = " << segment.squared_length() << std::endl ;
						//ransacThres = segment.squared_length() ;
						intPoint = IntersectionWeightedCentroid( segment, 
																 point_vector, 
																 circle_vector, 
																 surface.distanceSigma() ) ;
					}
					else {
						ransacThres = segment.squared_length()*surface.ransacDistThres()*surface.ransacDistThres() ; // Squared ransac threshold (ransac compares squared distances)
						// --- Debug (Start) ---
						/*if ( ransacThres < 0.001*0.001 )
							std::cout << "Very small ransac value = " << ransacThres << std::endl ;*/
						// --- Debug  (End)  ---
						intPoint = IntersectionRansacWeightedCentroid(	segment,
																		point_vector, 
																		circle_vector, 
																		surface.distanceSigma(), 
																		ransacThres ) ;
					}
					break ;
				case AABB_Tree_Splats_Mesher_Oracle::RansacMixtureOfGaussiansMaximum :
					ransacThres = segment.squared_length()*surface.ransacDistThres()*surface.ransacDistThres() ; // Squared ransac threshold (ransac compares squared distances)
					intPoint = IntersectionRansacMixtureOfGaussiansMaximum( segment, 
																			point_vector, 
																			circle_vector, 
																			surface.distanceSigma(),
																			surface.mogSigma(),
																			surface.nBins(),
																			ransacThres ) ;		
					break ;
				case AABB_Tree_Splats_Mesher_Oracle::RansacMaximumScore :
					// Threshold depending on the segment length
					ransacThres = segment.squared_length()*surface.ransacDistThres()*surface.ransacDistThres() ; // Squared ransac threshold (ransac compares squared distances)
					intPoint = IntersectionRansacMaximumScore( point_vector, circle_vector, ransacThres ) ;
					break ;					
				default:
					std::cerr << "Unknown intersection type!" << std::endl ;
					return NULL_POINT ;
			}

			// Optimize the position along the segment according to photoconsistency
			if ( surface.optimType() >= 0 ) {
					
				// Compute the mean normal of intersected circles
				Vector_3 normal = meanNormal( circle_vector ) ;

				// Filter views according to normal and projection in the image plane
				std::vector< int > compatibleIndexes = AuxLib::compatibleImagesIndexes<Kernel>( surface.dataset(), intPoint, normal, surface.angleCompatibilityThres() ) ;
								
				if ( compatibleIndexes.size() >= 2 ) {
						
					// Set the optimization parameters					
					std::pair< double, double > limits = AuxLib::limitsOfDisplacementValues< Kernel >( segment, intPoint ) ;
					
					PhotoConsOptimizer< Kernel > pcOptim ;
					pcOptim.setStartPoint( intPoint ) ;
					pcOptim.setDirVector( AuxLib::normalize< Kernel >( segment.to_vector() ) ) ;
					pcOptim.setDataset( surface.dataset() ) ;
					pcOptim.setCompatibleViewsIndexes( compatibleIndexes ) ;
					pcOptim.setCorrelRadius( surface.correlRadius() ) ;
					pcOptim.setCorrelThres( surface.correlThres() ) ;
					pcOptim.setLimits( limits.first, limits.second ) ;
					pcOptim.setOptimType( surface.optimType() ) ;

					// --> Normal as the direction from the point to the best compatible camera
					BoostMatrix iv = AuxLib::invView( surface.dataset()->getView( compatibleIndexes[0] ) ) ;
					Point_3 camPos( iv( 0, 3 ), iv( 1, 3 ), iv( 2, 3 ) ) ;
					pcOptim.setInitialNormal( AuxLib::normalize< Kernel >( Vector_3( intPoint, camPos ) ) ) ;
					
					// --> Normal as the mean of the normals of the intersected triangles
					// pcOptim.setInitialNormal( normal ) ; // Only used in optim type == DoubleOptim3DOF
										
					//std::cout << "[AABBTriangleSoupPhotoCons] Point before optimization = " << intPoint << std::endl ;
					Point_3 oldIntPoint = intPoint ;
					bool optimResult = pcOptim.PhotoConsOptimizer< Kernel >::optimizePhotoCons() ;

					// Any change?
					//if( intPoint != Point_3() ) {
					if( optimResult ) {
						MesherSuccesfullyOptimizedPts++ ;
						//std::cout << "[AABBTriangleSoupPhotoCons] (Successful) Point after optimization  = " << intPoint << std::endl ;
						return pcOptim.getPoint() ;
					}
					//std::cout << "[AABBTriangleSoupPhotoCons] (Failed) Point after optimization  = " << intPoint << std::endl ;
				
				}
					
			}
			else {
				// std::cout << "Intersection point computed: " << intPoint << std::endl ;
				return intPoint ;
			}

		}
		
		// Intersection between the surface and a ray
		Intersection operator()(const Ray_3& ray) const
		{
			MesherIntersectionCounterRays++ ;
			MesherIntersectionCounter++ ;
			//surface.increment_number_of_intersection_tests() ;
			
			// Change the query to reuse the machinery already coded for segments
			Segment_3 querySegment = Segment_3( ray.source(), 
												ray.source()+( ray.to_vector()*999999 ) ) ; // Create a huge segment following the query ray
			return (*this)( querySegment ) ;
		}
		
		// Intersection between the surface and a line
		Intersection operator()(const Line_3& line) const
		{
			MesherIntersectionCounterLines++ ;
			MesherIntersectionCounter++ ;
			////surface.increment_number_of_intersection_tests() ;
			
			// Change the query to reuse the machinery already coded for segments
			Point_3 pointOnLine = line.point() ;
			Segment_3 querySegment = Segment_3( pointOnLine + ( line.to_vector()*-999999 ), 
												pointOnLine+( line.to_vector()*999999 ) ) ; // Create a huge segment following the query line
			return (*this)( querySegment ) ;
		}


		Point_3 closeHole( const Surface_3& surface, const Segment_3& newSegment, const Segment_3& segment ) const {
			// Vector direction
			Vector_3 segDir = AuxLib::normalize< Kernel >( newSegment.to_vector() ) ;

			// Get orthonormal base
			Plane_3 plane( newSegment.source(), segDir ) ;
			Vector_3 u = AuxLib::normalize< Kernel >( plane.base1() ) * surface.holeFillingRadius() ;
			Vector_3 v = AuxLib::normalize< Kernel >( plane.base2() ) * surface.holeFillingRadius() ;
					
			// Get N segments parallel to the query segment around a circle
			std::vector< Segment_3 > holeFillingQuerySegments ;
			double partPi = ( 2*CGAL_PI ) / static_cast<double>( surface.holeFillingNumQueries() ) ;
			FT radius = surface.holeFillingRadius() ;
			for ( double j = 0; j < ( 2*CGAL_PI ); j=j+partPi ) {
				// Starting point (plane at query segments' start point)
				Point_3 s = newSegment.source() + radius * cos( j ) * u + radius * sin( j ) * v ;
				// Starting point (plane at query segments' end point)
				Point_3 e = newSegment.target() + radius * cos( j ) * u + radius * sin( j ) * v ;
				// Save the segment
				holeFillingQuerySegments.push_back( Segment_3( s, e ) ) ;
			}

			// For each segment, compute the intersection with other circles
			std::vector< Segment_3 >::iterator itHfS ;
			std::vector< Circle_3 > circlesHf ;
			std::vector< Circle_3 >::iterator found ;
			std::vector< Point_3 > pointsHf ;
			int numHfSegIntersect = 0 ;
			//std::cout << "Original Segment: " << segment << std::endl ;
			for ( itHfS = holeFillingQuerySegments.begin(); itHfS != holeFillingQuerySegments.end(); ++itHfS ) {
				//std::cout << "Hole filling query Segment: " << *itHfS << std::endl ;
				if ( surface.tree()->do_intersect( *itHfS ) ) {
					numHfSegIntersect++ ;
					//std::cout << "Segment intersect the splats" << std::endl ;
							
					std::vector< Obj_prim_id > obj_prim_vect ;							
					surface.tree()->all_intersections( *itHfS, std::back_inserter( obj_prim_vect ) ) ;
					typename std::vector< Obj_prim_id >::iterator it ;

					// Insert circles not already in the list
					for ( it = obj_prim_vect.begin(); it != obj_prim_vect.end(); ++it ) {

						// Check if the current circle is in the list
						found = std::find( circlesHf.begin(), circlesHf.end(), *(it->second) ) ;
						if ( found == circlesHf.end() ) {
							// Not in the list, check if the intersection with the original segment is a point
														
							// Compute the intersection points with the query segment by computing the intersection with the plane the disc represents
							Object o = CGAL::intersection( (*(it->second)).supporting_plane(), newSegment ) ;
							Point_3 pInt ;
							if ( CGAL::assign( pInt, o ) ) {						
								pointsHf.push_back( pInt ) ;
								circlesHf.push_back( *(it->second) ) ;								
							}
						}
					}

				}
			}
			//std::cout << "Hole filling: num circles = " << circlesHf.size() << std::endl ;

			// At least half of the segment queries must interect the surface
			if ( numHfSegIntersect > ( surface.holeFillingNumQueries() / 2 ) ) { 
				
				// Compute the intersection point			
				//std::cout << "Hole filling: num points = " << pointsHf.size() << std::endl ; 
				//std::cout << "Hole filling: num circles = " << circlesHf.size() << std::endl ; cin.get() ;
				Point_3 intPoint = computeIntersectionPoint( surface, 
															 newSegment, 
															 pointsHf, 
															 circlesHf ) ;

				//return Object() ;
				if ( intPoint.x() == 99999999 )						
					return NULL_POINT ;
				//else
				// Resulting intersection point has to be on the original query segment
				else {
					//if( segment.has_on( intPoint ) ) { 
					// if( segment.has_on( intPoint ) )
					FT dist = CGAL::sqrt( CGAL::squared_distance( intPoint, newSegment.source() ) ) ;
					FT orSegLength = CGAL::sqrt( segment.squared_length() ) ;
					FT hDifSeg = ( CGAL::sqrt( newSegment.squared_length() ) - orSegLength ) / 2 ;
					if ( segment == newSegment || ( ( dist - hDifSeg > 0 ) && ( dist - hDifSeg < orSegLength ) ) )
						return intPoint ;				
					}

			}
			else {
				return NULL_POINT ;
			}

		}


		/* Auxiliar functions */

		Segment_3 minLengthSegment( const Segment_3& segment, FT minLength ) const {

			if ( minLength <= 0 ) 
				return segment ;

			double curLength = CGAL::sqrt( segment.squared_length() ) ;
			if ( curLength >= minLength ) 
				return segment ;
			else {
				double rest = minLength - curLength ;
				double halfRest = rest / 2 ;

				// Project the new points (source/target) along the segment
				Vector_3 normDir = AuxLib::normalize<Kernel>( segment.to_vector() ) ;
				Point_3 newSource = segment.source() + ( -normDir * halfRest ) ;
				Point_3 newTarget = segment.target() + (  normDir * halfRest ) ;

				return Segment_3( newSource, newTarget ) ;
			}

		}

		// Compute mean normal from triangle set
		Vector_3 meanNormal( const std::vector< Circle_3 >& circles ) const {

			typename std::vector< Circle_3 >::const_iterator it ;
			double x_acc = 0.0, y_acc = 0.0, z_acc = 0.0 ;
			
			for ( it = circles.begin(); it != circles.end(); ++it ) {
			//for ( int i = 0; i < circles.size(); i++ ) {
				
				Vector_3 normal = it->supporting_plane().orthogonal_vector() ;
				//Vector_3 normal = circles.at( i ).supporting_plane().orthogonal_vector() ; 
				Vector_3 normalUnit = AuxLib::normalize<Kernel>( normal ) ;
					
				x_acc += normal.x() ;
				y_acc += normal.y() ;
				z_acc += normal.z() ;

			}
			Vector_3 vect_acc( x_acc, y_acc, z_acc ) ;
			vect_acc = vect_acc / static_cast< double >( circles.size() ) ;

			//std::cout<< "[AABB_Tree_Splats_Mesher_Oracle] Mean Normal before normalize = " << vect_acc << std::endl ;
	
			// Normalize vector (if needed)
			vect_acc = AuxLib::normalize< Kernel >( vect_acc ) ;		
		
			//std::cout<< "[AABB_Tree_Splats_Mesher_Oracle] Mean Normal result = " << vect_acc << std::endl ;
	
			return vect_acc ;

		}

		/* Intersection types */

		// The intersection is the centroid
		Point_3 IntersectionCentroid( const std::vector< Point_3 >& points ) const {
			return CGAL::centroid( points.begin(), points.end() ) ;
		}

		// The intersection is a weighted centroid
		// If a score != -1 is present in the circles, it will be used to weight the gaussians
		Point_3 IntersectionWeightedCentroid( const Segment_3& segment,
											 const std::vector< Point_3 >& points, 
											 const std::vector< Circle_3 >& circles, 
											 const double& sigma ) const {
			typename std::vector< Point_3 >::const_iterator itP ;
			typename std::vector< Circle_3 >::const_iterator itC ;
			FT weightedDispSum = 0 ;
			FT totalWeight = 0 ;

			Point_3 ps = segment.source() ;
			Point_3 pe = segment.target() ;

			Vector_3 dir = segment.to_vector() ;
			// Unitize
			dir = dir / CGAL::sqrt( dir.squared_length() ) ;

			itC = circles.begin() ;
			for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
				// Displacement of the point along the segment
				FT disp = CGAL::sqrt( CGAL::squared_distance( ps, *itP ) ) ;

				// Compute the distance from the intersection to the center of the Circle_3
				FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

				// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
				//FT weight = centeredGaussian1D( 0, distCent, sigma ) ;
				//FT weight = centeredGaussian1D( distCent, 0, sigma ) ;
				FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
				FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;
				//FT weight = weighting( distCent, curSigma ) ;

				if ( itC->score() > 0 ) 
					weight = weight * itC->score() ;

				//FT weight = distCent ;

				// --- Debug (Start) --- 
				//std::cout << "distCent = " << distCent << std::endl ;
				//std::cout << "radius = " << CGAL::sqrt( itC->squared_radius() ) << std::endl ;
				//std::cout << "sigma = " << sigma << std::endl ;
				//std::cout << "curSigma = " << curSigma << std::endl ;
				//std::cout << "weight = " << weight << std::endl ;
				//std::cin.get() ;
				// --- Debug  (End)  ---

				// Weight the current displacement
				weightedDispSum += disp * weight ;

				// Increase the total weight
				totalWeight += weight ;
			}

			FT meanDisp = weightedDispSum / totalWeight ;

			//if ( CGAL::sqrt( CGAL::squared_distance( ps, meanPoint ) ) < CGAL::sqrt( segment.squared_length() ) ) {
			if ( meanDisp < CGAL::sqrt( segment.squared_length() ) ) {
				Point_3 meanPoint = ps + ( meanDisp * dir ) ;
				return meanPoint ;
			}
			else {			

				//itC = circles.begin() ;
				//for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
				//	FT disp = CGAL::sqrt( CGAL::squared_distance( ps, *itP ) ) ;
				//	FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;
				//	FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
				//	FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;
				//	if ( itC->score() > 0 ) 
				//		weight = weight * itC->score() ;
				//	
				//	// --- Debug (Start) --- 
				//	std::cout << "distCent = " << distCent << std::endl ;
				//	std::cout << "radius = " << CGAL::sqrt( itC->squared_radius() ) << std::endl ;
				//	std::cout << "sigma = " << sigma << std::endl ;
				//	std::cout << "curSigma = " << curSigma << std::endl ;
				//	std::cout << "weight = " << weight << std::endl ;
				//	//std::cin.get() ;
				//}
				//std::cout << "Number of intersection points = " << points.size() << std::endl ;
				//std::cout << "Number of intersected circles = " << circles.size() << std::endl ;
				//std::cout << "weightedDispSum = " << weightedDispSum << std::endl ;
				//std::cout << "totalWeight = " << totalWeight << std::endl ;
				//std::cout << "meanDisp = "<< meanDisp << std::endl ;
				//std::cout << "seg length = "<< CGAL::sqrt( segment.squared_length() ) << std::endl ;
				//std::cin.get() ;				
				return NULL_POINT ;
			}
		}

		// The intersection is the maximum of a mixture of gaussians model
		// If a score != -1 is present in the circles, it will be used to weight the gaussians
		Point_3 IntersectionMixtureOfGaussiansMaximum(	const Segment_3& segment,
														const std::vector< Point_3 >& points, 
														const std::vector< Circle_3 >& circles, 
														const double& distanceSigma,
														const double& mogSigma,
														const int& nBins ) const {
			typename std::vector< Point_3 >::const_iterator itP ;
			typename std::vector< Circle_3 >::const_iterator itC ;
			std::vector< FT > centers ;
			std::vector< FT > weights ;
			
			Point_3 ps = segment.source() ;
			Point_3 pe = segment.target() ;

			Vector_3 dir = segment.to_vector() ;
			// Unitize
			dir = dir / CGAL::sqrt( dir.squared_length() ) ;

			//FT totalSqWeightSum = 0 ;
			FT totalWeightSum = 0 ;

			itC = circles.begin() ;
			for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
				// Displacement of the point along the segment
				FT disp = CGAL::sqrt( CGAL::squared_distance( ps, *itP ) ) ;
				centers.push_back( disp ) ; // Store it on the centers vector

				// Compute the distance from the intersection to the center of the Circle_3
				FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

				// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
				// FT weight = centeredGaussian1D( distCent, 0, distanceSigma ) ;
				FT curSigma = CGAL::sqrt( itC->squared_radius() ) * distanceSigma ; // Now distancesSigma is just a multiplicative factor!!
				FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;				
				//FT weight = weighting( distCent, curSigma ) ;				

				if ( itC->score() > 0 ) 
					weight = weight * itC->score() ;

				weights.push_back( weight ) ;

				// Increase the total weight
				totalWeightSum += weight ;
				//FT weightSq = weight * weight ;
				//totalSqWeightSum += weightSq ;
			}
			//// Norm2 of the weights vector
			//FT normWeights = CGAL::sqrt( totalSqWeightSum ) ;

			// Make the weights vector sum up to 1
			typename std::vector< FT >::iterator itW ;
			for ( itW = weights.begin(); itW != weights.end(); ++itW ) {
				*itW = *itW / totalWeightSum ;
			}			
			
			// --- Debug (Start) ---
			/*double sumWeights = 0 ;
			for ( itW = weights.begin(); itW != weights.end(); ++itW ) {
				sumWeights += *itW ;				
			}	
			std::cout << "SumWeights = " << sumWeights << std::endl ;

			CImg< double > mogEvalImage( nBins, 1 ) ;
			std::vector< FT >::iterator itCent ;
			std::vector< CImg< double > > singleGaussianImages ;
			std::vector< CImg< double > >::iterator itI ;
			for ( int i = 0; i < centers.size(); i++ ) {
				CImg< double > newImg( nBins, 1 ) ;
				singleGaussianImages.push_back( newImg ) ;
			}*/
			// --- Debug  (End)  ---

			// Evaluate the MoG at nBins locations along the segment and get the maximum
			double max = 0 ;
			double maxDisp = 0 ;
			double segmentLength = CGAL::sqrt( CGAL::squared_distance( ps, pe ) ) ;
			double delta = segmentLength / nBins ;
			for( int i = 0; i < nBins; i++ ) {
				double disp = i*delta ;
				//double eval = evalMoG( disp, centers, weights, mogSigma ) ;
				double eval = evalMoG( disp, centers, weights, segmentLength*mogSigma ) ; // Now mogSigma is a multiplicative factor!!

				if ( eval > max ) {
					max = eval ;
					maxDisp = disp ;
				}

				// --- Debug (Start) ---
				//std::cout << "Eval = " << eval << std::endl ;
				/*mogEvalImage.atXY( i, 0 ) = eval ;
								
				itI = singleGaussianImages.begin() ;						
				for ( itCent = centers.begin(); itCent != centers.end(); ++itCent, ++itI ) {
					itI->atXY( i, 0 ) = centeredGaussian1D( disp, *itCent, mogSigma ) ;
				}			*/
				// --- Debug  (End)  ---			
			}

			// --- Debug (Start) ---
			//CImg< double > mogGraph( 500, 400, 1, 3, 0 ) ;
			//CImg< double > mogGraphMOG( 500, 400, 1, 3, 0 ) ;
			//CImgDisplay mogDisplay( mogGraph, "MoG Evaluation along the segment" ) ;
			//const unsigned char red[] = { 255,0,0 } ;
			//const unsigned char green[] = { 0, 255,0 } ;
			//const unsigned char blue[] = { 0, 0, 255 } ;
			//
			//// Draw centers of the gaussians
			//for ( itCent = centers.begin(); itCent != centers.end(); ++itCent ) {
			//	int centIm = ( *itCent / segmentLength ) * mogGraph.width() ;

			//	//mogGraph.draw_point( centIm, 399, blue ) ;
			//	mogGraph.draw_line( centIm, 0, centIm, mogGraph.width()-1, blue ) ;
			//}
			//
			//// Draw single gaussians
			//for ( itI = singleGaussianImages.begin(); itI != singleGaussianImages.end(); ++itI ) {
			//	mogGraph.draw_graph( *itI, green, 1, 1, 0, 0, 0 ) ;
			//}

			//// Draw mixture
			//mogGraph.draw_graph( mogEvalImage, red, 1, 1, 0, 0, 0 ) ;
			//mogGraphMOG.draw_graph( mogEvalImage, red, 1, 1, 0, 0, 0 ) ;
			//
			//mogGraph.save_bmp( "temp.bmp" ) ;
			//mogGraphMOG.save_bmp( "temp_mog.bmp" ) ;
			//while( !mogDisplay.is_key() ) {
			//	mogDisplay.wait() ;
			//}
			// --- Debug  (End)  ---

			Point_3 maxPoint = ps + maxDisp * dir ;

			return maxPoint ;
		}

		// The intersection is a centroid computed using RANSAC
		Point_3 IntersectionRansacCentroid( const std::vector< Point_3 >& points, 
									  const double& ransacDistThres ) const {
			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				bool exec =	ransac::ransacFitCentroid< Kernel >( points, 
																 ransacDistThres, 
																 centroidRansac,
																 inliers ) ;

				if ( !exec ) {
					//std::cerr << "Ransac problem" << std::endl ;
					//return CGAL::centroid( points.begin(), points.end() ) ;
					return NULL_POINT ;
				}

				FT xSum = 0, ySum = 0, zSum = 0 ;
				std::vector< int >::iterator inIt ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					xSum += points[ *inIt ].x() ;
					ySum += points[ *inIt ].y() ;
					zSum += points[ *inIt ].z() ;
				}
				return Point_3( xSum / inliers.size(), ySum / inliers.size(), zSum / inliers.size() ) ;

			}
			else {
				if (  points.size() == 2 && CGAL::squared_distance( points[0], points[1] ) < ransacDistThres )
					return CGAL::centroid( points.begin(), points.end() ) ;
				else				
					return NULL_POINT ;
					//return CGAL::centroid( points.begin(), points.end() ) ;
			}
		}

		// The intersection is a centroid computed using RANSAC
		Point_3 IntersectionRansacWeightedCentroid( const Segment_3& segment,
													const std::vector< Point_3 >& points, 
													const std::vector< Circle_3 >& circles, 
													const double& sigma,
													const double& ransacDistThres ) const {

			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				// Compute weights
				std::vector< Point_3 >::const_iterator itP ;
				std::vector< Circle_3 >::const_iterator itC ;
				std::vector< double > weights ;
				itC = circles.begin() ;				
				for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
					// Compute the distance from the intersection to the center of the Circle_3
					FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

					// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
					//FT weight = centeredGaussian1D( 0, distCent, sigma ) ;
					//FT weight = centeredGaussian1D( distCent, 0, sigma ) ;
					FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
					FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;	
					//FT weight = weighting( distCent, curSigma ) ;	

					if ( itC->score() > 0 )				
						weight = weight * itC->score() ; // Score also with weight

					//std::cout << "[IntersectionRansacWeightedCentroid] Radius = " << CGAL::sqrt( itC->squared_radius() ) << "   distCent = " << distCent << "   Sigma = " << curSigma << "   Weight = " << weight << std::endl ;

					weights.push_back( weight ) ;
				}
				
				bool exec =	ransac::ransacFitWeightedCentroid< Kernel >( points,
																		 weights,
																		 ransacDistThres,
																		 centroidRansac,
																		 inliers ) ;

				if ( !exec ) {
					//std::cerr << MesherIntersectionCounter << " Ransac problem..." << std::endl ;
					//return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
					return NULL_POINT ;
				}

				FT xSum = 0, ySum = 0, zSum = 0 ;
				std::vector< int >::iterator inIt ;
				std::vector< Point_3 > inlierPoints ;
				std::vector< Circle_3 > inlierCircles ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					inlierPoints.push_back( points[ *inIt ] ) ;					
					inlierCircles.push_back( circles[ *inIt ] ) ;
				}
				return IntersectionWeightedCentroid( segment, inlierPoints, inlierCircles, sigma ) ;

				// return centroidRansac ;

			}
			else {
				//std::cerr << MesherIntersectionCounter << " Not enough points (" << points.size() << ")" << std::endl ;
				if (  points.size() == 2 && CGAL::squared_distance( points[0], points[1] ) < ransacDistThres )
					return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
				else	
					//return points[0];
					return NULL_POINT ;
					//return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;					
			}
		}

		Point_3 IntersectionRansacMixtureOfGaussiansMaximum(	const Segment_3& segment,
																const std::vector< Point_3 >& points, 
																const std::vector< Circle_3 >& circles, 
																const double& distanceSigma,
																const double& mogSigma,
																const int& nBins,
																const double& distThres ) const {
			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				bool exec =	ransac::ransacFitCentroid< Kernel >( points, 
																 distThres, 
																 centroidRansac,
																 inliers ) ;

				if ( !exec ) {
					//std::cerr << "Ransac problem" << std::endl ;
					//return CGAL::centroid( points.begin(), points.end() ) ;
					//return IntersectionMixtureOfGaussiansMaximum( segment, points, circles, distanceSigma, mogSigma, nBins ) ;	
					return NULL_POINT ;
				}

				FT xSum = 0, ySum = 0, zSum = 0 ;
				std::vector< int >::iterator inIt ;
				std::vector< Point_3 > inlierPoints ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					inlierPoints.push_back( points[ *inIt ] ) ;					
				}
				return IntersectionMixtureOfGaussiansMaximum( segment, inlierPoints, circles, distanceSigma, mogSigma, nBins ) ;	

			}
			else {
				return IntersectionMixtureOfGaussiansMaximum( segment, points, circles, distanceSigma, mogSigma, nBins ) ;	
				//return NULL_POINT ;
			}
		}



		// The intersection is the inlier point (RANSAC) corresponding to the patch with largest score
		Point_3 IntersectionRansacMaximumScore( const std::vector< Point_3 >& points, 
											const std::vector< Circle_3 >& circles,
											const double& distThres ) const {
			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				bool exec =	ransac::ransacFitCentroid< Kernel >( points, 
																 distThres, 
																 centroidRansac,
																 inliers ) ;

				if ( !exec ) {
					//std::cerr << "Ransac problem" << std::endl ;
					// All are inliers
					for ( int i = 0; i < points.size(); i++ ) {
						inliers.push_back( i ) ;
					}
				}

				std::vector< int >::iterator inIt ;
				double maxScore = 0 ;
				int maxScoreInd = -1 ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					if ( circles[*inIt].score() > maxScore ) {
						maxScoreInd = *inIt ;
						maxScore = circles[*inIt].score() ;
					}
				}
				return points[ maxScoreInd ] ;
			}
			else {
				if ( points.size() > 1 ) {
					if ( circles[0].score() > circles[1].score() ) {
						return points[0] ;
					}
					else {
						return points[1] ;
					}
				}
				else {
					return points[0] ;
				}
			}
		}



		Point_3 IntersectionRansacWeightedCentroidScore( const Segment_3& segment,
														 const std::vector< Point_3 >& points, 
														 const std::vector< Circle_3 >& circles, 
														 const double& sigma,
														 const double& distThres ) const {
			Point_3 centroidRansac ;
			std::vector< int > inliers ;
			if ( points.size() > 2 ) {
				// Compute weights
				std::vector< Point_3 >::const_iterator itP ;
				std::vector< Circle_3 >::const_iterator itC ;
				std::vector< double > weights ;
				itC = circles.begin() ;				
				for ( itP = points.begin(); itP != points.end(); ++itP, ++itC ) {
					// Compute the distance from the intersection to the center of the Circle_3
					FT distCent = CGAL::sqrt( CGAL::squared_distance( *itP, itC->center() ) ) ;

					// Compute the weight of the intersection as the distance from the center of the intersected Circle_3 primitive
					//FT weight = centeredGaussian1D( 0, distCent, sigma ) ;
					//FT weight = centeredGaussian1D( distCent, 0, sigma ) ;
					FT curSigma = CGAL::sqrt( itC->squared_radius() ) * sigma ; // Now sigma is a multiplicative factor!!
					FT weight = centeredGaussian1D( distCent, 0, curSigma ) ;	
					//FT weight = weighting( distCent, curSigma ) ;	
					
									
					//std::cout << "[IntersectionRansacWeightedCentroid] Radius = " << CGAL::sqrt( itC->squared_radius() ) << "   distCent = " << distCent << "   Sigma = " << curSigma << "   Weight = " << weight << std::endl ;

					weights.push_back( weight ) ;
				}
				
				bool exec =	ransac::ransacFitWeightedCentroid< Kernel >( points,
																		 weights,
																		 distThres,
																		 centroidRansac,
																		 inliers ) ;

				if ( !exec ) {
					std::cerr << MesherIntersectionCounter << " Ransac problem..." << std::endl ;
					return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
				}

				FT xSum = 0, ySum = 0, zSum = 0 ;
				std::vector< int >::iterator inIt ;
				std::vector< Point_3 > inlierPoints ;
				for ( inIt = inliers.begin(); inIt != inliers.end(); ++inIt ) {
					inlierPoints.push_back( points[ *inIt ] ) ;					
				}
				return IntersectionWeightedCentroid( segment, inlierPoints, circles, sigma ) ;

				// return centroidRansac ;

			}
			else {
				//std::cerr << MesherIntersectionCounter << " Not enough points (" << points.size() << ")" << std::endl ;
				return IntersectionWeightedCentroid( segment, points, circles, sigma ) ;
			}
		
		}

	} ; // end class Construct_intersection