/* 
	The functions developed in this file are a conversion to CGAL from the ones in Wild Magic 5.9 library 
	Web page: http://www.geometrictools.com/
*/


#ifndef SQUARED_DISTANCE_SEGMENT_TRIANGLE_3_H
#define SQUARED_DISTANCE_SEGMENT_TRIANGLE_3_H

#include <CGAL/basic.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Point_3.h>
#include <ClosestPoint/closest_point_triangle_3.h>
#include <iostream> // Debug only, should be removed

namespace CGAL {

namespace internal {


template <typename K>
typename K::FT line_segment_distance_with_disp( const Line_3<K> l, 
												const Point_3<K> lineOrigin, 
												const Segment_3<K> s, 
												typename K::FT& disp )
{
	typedef typename K::FT FT ;

	Vector_3<K> dummyEnd( s.end()[0], s.end()[1], s.end()[2] ) ;
	Vector_3<K> dummyStart( s.start()[0], s.start()[1], s.start()[2] ) ;
	Vector_3<K> midPtVec = FT( 0.5 ) * ( dummyStart + dummyEnd )  ;	// Midpoint of the segment	
	Point_3<K> segCenter( midPtVec[0], midPtVec[1], midPtVec[2] ) ;

    Vector_3<K> diff = lineOrigin - segCenter ;
	Vector_3<K> segDir = s.to_vector() / CGAL::sqrt( s.squared_length() ) ;
	FT a01 = ( -l.to_vector() ) * segDir ;
    FT b0 = diff * l.to_vector() ;
    FT c = diff.squared_length() ;
    FT det = CGAL::abs( FT( 1 ) - a01*a01 ) ;
    FT b1, s0, s1, sqrDist, extDet ;

    if ( det >= 1e-8 )
    {
        // The line and segment are not parallel.
        b1 = ( -diff ) * segDir ;
        s1 = a01*b0 - b1 ;
        extDet = ( CGAL::sqrt( s.squared_length() ) * FT( 0.5 ) ) * det ;

        if (s1 >= -extDet)
        {
            if (s1 <= extDet)
            {
                // Two interior points are closest, one on the line and one
                // on the segment.
                FT invDet = FT( 1 ) / det ;
                s0 = ( a01*b1 - b0 )*invDet ;
                s1 *= invDet ;
                sqrDist = s0*(s0 + a01*s1 + FT(2)*b0) +
                    s1*(a01*s0 + s1 + FT(2)*b1) + c;
            }
            else
            {
                // The endpoint e1 of the segment and an interior point of
                // the line are closest.
                s1 = CGAL::sqrt( s.squared_length() ) * FT(0.5) ;
                s0 = -(a01*s1 + b0);
				sqrDist = -s0*s0 + s1*(s1 + FT(2)*b1) + c;				
            }
        }
        else
        {
            // The end point e0 of the segment and an interior point of the
            // line are closest.
            s1 = - ( CGAL::sqrt( s.squared_length() ) * FT(0.5) ) ;
            s0 = -(a01*s1 + b0);
            sqrDist = -s0*s0 + s1*(s1 + FT(2.0)*b1) + c;
        }
    }
    else
    {
        // The line and segment are parallel.  Choose the closest pair so that
        // one point is at segment center.
        s1 = FT( 0 ) ;
        s0 = -b0 ;
        sqrDist = b0*s0 + c ;
    }

	// --- Debug (start) ---
	/*Point_3<K> closestOnLine = lineOrigin + s0*l.to_vector() ;
    Point_3<K> closestOnSegment = segCenter + s1*segDir ;

	std::cout << " Closest on line = " << closestOnLine << std::endl ;
	std::cout << " Closest on segment = " << closestOnSegment << std::endl ;*/
	
	// Overwrite!
	// sqrDist = squared_distance( closestOnLine, closestOnSegment ) ;
	// --- Debug (end) ---


    disp = s0;
    
    // Account for numerical round-off errors.
    if ( sqrDist < FT( 0 ) )
    {
        sqrDist = FT( 0 );
    }
    return sqrDist ;
}


template <class K>
void generate_complement_basis (	Vector_3<K>& u, 
									Vector_3<K>& v,
									const Vector_3<K>& w)
{
	typedef typename K::FT FT ;
    FT invLength ;

	FT u0, u1, u2, v0, v1, v2 ;

    if ( CGAL::abs( w[0] ) >= CGAL::abs( w[1] ) )
    {
        // W.x or W.z is the largest magnitude component, swap them
        invLength = 1 / CGAL::sqrt( w[0]*w[0] + w[2]*w[2] ) ;
        u0 = -w[2]*invLength ;
        u1 = FT( 0 ) ;
        u2 = +w[0]*invLength ;
        v0 = w[1]*u2;
        v1 = w[2]*u0 - w[0]*u2;
        v2 = -w[1]*u0;
    }
    else
    {
        // W.y or W.z is the largest magnitude component, swap them
        invLength = 1/ CGAL::sqrt( w[1]*w[1] + w[2]*w[2] ) ;
        u0 = FT( 0 ) ;
        u1 = +w[2]*invLength;
        u2 = -w[1]*invLength;
        v0 = w[1]*u2 - w[2]*u1;
        v1 = -w[0]*u2;
        v2 = w[0]*u1;
    }

	u = Vector_3<K>( u0, u1, u2 ) ;
	v = Vector_3<K>( v0, v1, v2 ) ;

}


template <class K>
typename K::FT	
squared_distance_aux(	const Line_3<K>& l,
						const Triangle_3<K>& t,
						const Point_3<K> origin,
						typename K::FT& disp ) {

	typedef typename K::FT FT ;

	// Test if line intersects triangle.  If so, the squared distance is zero.
    Vector_3<K> edge0 = t.vertex(1) - t.vertex(0) ;
    Vector_3<K> edge1 = t.vertex(2) - t.vertex(0) ;
    Vector_3<K> normal = CGAL::cross_product( edge0, edge1 ) ;
    Vector_3<K> line_vec = l.to_vector();
    FT NdD = normal * line_vec;
    if ( CGAL::abs( NdD ) > 1e-08 )
    {
        // The line and triangle are not parallel, so the line intersects
        // the plane of the triangle.
        Vector_3<K> diff = origin - t.vertex(0) ;
        Vector_3<K> U, V ;
        generate_complement_basis( U, V, l.to_vector() ) ;
        FT UdE0 = U * edge0 ;
        FT UdE1 = U * edge1 ;
        FT UdDiff = U * diff ;
        FT VdE0 = V * edge0 ;
        FT VdE1 = V * edge1 ;
        FT VdDiff = V * diff ;
        FT invDet = FT(1) / ( UdE0*VdE1 - UdE1*VdE0 ) ;

        // Barycentric coordinates for the point of intersection.
        FT b1 = (VdE1*UdDiff - UdE1*VdDiff)*invDet;
        FT b2 = (UdE0*VdDiff - VdE0*UdDiff)*invDet;
        FT b0 = FT(1) - b1 - b2;

        if ( b0 >= FT(0) && b1 >= FT(0) && b2 >= FT(0) )
        {
            // Line parameter for the point of intersection.
            FT DdE0 = l.to_vector() * edge0 ;
            FT DdE1 = l.to_vector() * edge1 ;
            FT DdDiff = l.to_vector() * diff ;
            disp = b1*DdE0 + b2*DdE1 - DdDiff ;
			
			return FT( 0 ) ;
        }
    }

    // Either (1) the line is not parallel to the triangle and the point of
    // intersection of the line and the plane of the triangle is outside the
    // triangle or (2) the line and triangle are parallel.  Regardless, the
    // closest point on the triangle is on an edge of the triangle.  Compare
    // the line to all three edges of the triangle.
    FT sqrDist = 99999999999.9 ;
	FT dispTmp = 0.0 ;
    for (int i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
    {
        Segment_3<K> segment( t.vertex(i0), t.vertex(i1) ) ;

		
		FT testSqrDist = CGAL::squared_distance( l, segment ) ;

        FT sqrDistTmp = line_segment_distance_with_disp( l, origin, segment, dispTmp ) ;
		if (sqrDistTmp < sqrDist) {
            sqrDist = sqrDistTmp ;
            disp = dispTmp ;            
        }
    }

    return sqrDist;

}



} // End namespace internal



template <class K>
typename K::FT
squared_distance(	const Segment_3<K> & s,
					const Triangle_3<K> & t ) {

	// Get segment's midpoint
	typedef typename K::FT FT ;

	// Vector_3<K> dir = s.to_vector() ; // Direction vector 
	Vector_3<K> dir = s.to_vector() / CGAL::sqrt( s.squared_length() ) ; // Normalized direction vector
	Vector_3<K> dummyEnd( s.end()[0], s.end()[1], s.end()[2] ) ;
	Vector_3<K> dummyStart( s.start()[0], s.start()[1], s.start()[2] ) ;
	Vector_3<K> midPtVec = FT( 0.5 ) * ( dummyStart + dummyEnd )  ;	// Midpoint of the segment	
	Point_3<K> midPt( midPtVec[0], midPtVec[1], midPtVec[2] ) ;
	FT halfLength = FT(0.5)*CGAL::sqrt( s.squared_length() ) ; // Half the length of the segment

	typename K::Line_3 line( midPt, dir ) ;

	FT disp ;
	
	FT sqDist = internal::squared_distance_aux( line, t, midPt, disp ) ;
	// std::cout << "Distance to line = " << sqDist << std::endl ;

	if ( disp >= -halfLength ) {
		
		if ( disp > halfLength ) {
	
			// Closest point is the s.end()
			// std::cout << "Closest point is the s.end()" << std::endl ;
			Point_3<K> cp = CGAL::closest_point_3<K>( s.end(), t ) ;
			sqDist = CGAL::squared_distance( s.end(), cp ) ;
		}
	}
	else {
		// Closest point is s.start()
		// std::cout << "Closest point is the s.start()" << std::endl ;
		Point_3<K> cp = CGAL::closest_point_3<K>( s.start(), t ) ;
		sqDist = CGAL::squared_distance( s.start(), cp ) ;
	}
		
	return sqDist ;

}

template <class K>
typename K::FT
squared_distance(	const Triangle_3<K> & t,
					const Segment_3<K> & s	) {
	return squared_distance( s, t ) ;
}


} // End namespace CGAL


#endif // SQUARED_DISTANCE_SEGMENT_TRIANGLE_3_H