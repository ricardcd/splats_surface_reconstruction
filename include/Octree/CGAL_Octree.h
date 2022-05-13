// Octree code, closely based on the source example in flipcode:
// Web: http://www.flipcode.com/archives/Octree_Implementation.shtml
// From: Paul Nettle [midnight@FluidStudios.com]
//
// Author: Ricard Campos (rcampos@eia.udg.edu)

#include <vector>

namespace CGAL_Octree {


// Cubic bounding volume (NOT A BOUNDING BOX!)
template < class K >
class CubicBounds
{
public:
    // Some typedefs
    typedef typename K::FT FT ;
    typedef typename K::Point_3 Point_3 ;
    typedef typename K::Vector_3 Vector_3 ;


    // Constructors
    CubicBounds() : m_center(0.0,0.0,0.0), m_radius(0.0) { }

    CubicBounds( const Point_3 &center, const FT radius ) {
        m_center = center ;
        m_radius = radius ;
    }

    CubicBounds( const std::vector< Point_3 >& points ) {

        // Determine min/max of the given set of points
        FT min_x = points[0].x() ;
        FT min_y = points[0].y() ;
        FT min_z = points[0].z() ;

        FT max_x = points[0].x() ;
        FT max_y = points[0].y() ;
        FT max_z = points[0].z() ;


        typename std::vector<Point_3>::const_iterator it ;
        for ( it = points.begin();
              it != points.end();
              ++it )
        {
            if ( it->x() < min_x ) min_x = it->x() ;
            if ( it->y() < min_y ) min_y = it->y() ;
            if ( it->z() < min_z ) min_z = it->z() ;
            if ( it->x() > max_x ) max_x = it->x() ;
            if ( it->y() > max_y ) max_y = it->y() ;
            if ( it->z() > max_z ) max_z = it->z() ;
        }

        Point_3 minPt( min_x, min_y, min_z ) ;
        Point_3 maxPt( max_x, max_y, max_z ) ;

        // The radius of the volume (dimensions in each direction)
        Vector_3 radius = maxPt - minPt ;

        // Find the center of this space
        m_center = minPt + radius * 0.5 ;

        // We want a CUBIC space. By this, I mean we want a bounding cube, not
        // just a bounding box. We already have the center, we just need a
        // radius that contains the entire volume. To do this, we find the
        // maxumum value of the radius' X/Y/Z components and use that

        FT r = radius.x() ;
        if ( r < radius.y() ) r = radius.y() ;
        if ( r < radius.z() ) r = radius.z() ;
        m_radius = r ;
    }

    Point_3 center() const { return m_center ; }
    double radius() const { return m_radius ; }

    void setCenter( Point_3 center ) { m_center = center ; }
    void setRadius( FT radius ) { m_radius = radius ; }

private:
    Point_3     m_center ;
    double      m_radius ;
} ;



template < class K >
class  CGAL_Octree {
public:
    // Some typedefs
    typedef typename K::FT FT ;
    typedef typename K::Point_3 Point_3 ;
    typedef typename K::Vector_3 Vector_3 ;


    // Constructor/destructor
    CGAL_Octree() : m_pointCount(0), m_points(), m_center(0.0, 0.0, 0.0), m_radius(0.0) { memset(m_child, 0, sizeof(m_child)); }
    ~CGAL_Octree() {
        // delete [] m_child ;
    }

    std::vector< Point_3 > const points() const { return m_points ; }
    unsigned int numPoints() const { return m_points.size() ; }

    const bool build( const std::vector< Point_3 >& points,
                      const unsigned int maxDepth,
                      const CubicBounds< K > &bounds,
                      const unsigned int currentDepth = 0 )
    {
        //std::cout << "numPoints = " << points.size() << " maxDepth = " << maxDepth << " / currentDepth = " << currentDepth << std::endl ;

        // We are in a leaf because we recursed too deep into the tree
        //    (currentDepth >= maximumDepth)
        //
        //    NOTE: We specifically use ">=" for the depth comparison so that we
        //          can set the maximumDepth depth to 0 if we want a tree with
        //          no depth.
        if ( points.size() <= 1 || currentDepth >= maxDepth )
        {
                // Just store the points in the node, making it a leaf
                m_pointCount = (int)points.size() ;
                m_points = points ;
                return true ;
        }

        // We'll need this (see comment further down in this source)
        unsigned int childPointCounts[8];
        std::vector<unsigned int> codes( points.size(), 0 ) ;

        // Get the number of points
        int count = (int)points.size() ;

        // Classify each point to a child node
        for ( int i = 0; i < count; i++ )
        {
            // Current point
            Point_3 p = points[i] ;

            // Center of this node
            const Point_3 c = bounds.center() ;

            // Here, we need to know which child each point belongs to. To
            // do this, we build an index into the _child[] array using the
            // relative position of the point to the center of the current
            // node
            codes[i] = 0 ;
            if ( p.x() > c.x() ) codes[i] |= 1 ;
            if ( p.y() > c.y() ) codes[i] |= 2 ;
            if ( p.z() > c.z() ) codes[i] |= 4 ;

            // We'll need to keep track of how many points get stuck in each
            // child so we'll just keep track of it here, since we have the
            // information handy.
            childPointCounts[ codes[i] ]++ ;
        }

        // TODO: improve the definition of this unmodifyed table...
        std::vector< Vector_3 > boundsOffsetTable ;
        boundsOffsetTable.push_back( Vector_3( -0.5, -0.5, -0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( +0.5, -0.5, -0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( -0.5, +0.5, -0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( +0.5, +0.5, -0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( -0.5, -0.5, +0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( +0.5, -0.5, +0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( -0.5, +0.5, +0.5 ) ) ;
        boundsOffsetTable.push_back( Vector_3( +0.5, +0.5, +0.5 ) ) ;

        // Recursively call build() for each of the 8 children
        for (int i = 0; i < 8; i++)
        {
            // Don't bother going any further if there aren't any points for
            // this child
            if (!childPointCounts[i]) continue ;

            // Allocate the child
            m_child[i] = new CGAL_Octree<K>() ;

            // Allocate a list of points that were coded JUST for this child
            // only
            std::vector< Point_3 > newList ;

            // Go through the input list of points and copy over the points
            // that were coded for this child
            for ( int j = 0; j < count; j++ )
            {
                if ( codes[j] == i )
                {
                    newList.push_back( points[j] ) ;
                }
            }

            // At this point, we have a list of points that will belong to
            // this child node. We'll want to remove any point with a
            // duplicate 'n' in them...
            //
            // [We won't need to reallocate the list, since it's temporary]
//            int     newCount = 0;
//            for (j = 0; j < childPointCounts[i]; j++)
//            {
//                    // Remove duplicates from newList
//                    // ...
//                    // Keep track of the new count in 'newCount'
//            }

            // Generate a new bounding volume -- We do this with a touch of
            // trickery...
            //
            // We use a table of offsets. These offsets determine where a
            // node is, relative to it's parent. So, for example, if want to
            // generate the bottom-left-rear (-x, -y, -z) child for a node,
            // we use (-1, -1, -1).
            //
            // However, since the radius of a child is always half of its
            // parent's, we use a table of 0.5, rather than 1.0.
            //
            // These values are stored the following table. Note that this
            // won't compile because it assumes Points are structs, but you
            // get the idea.

            // Calculate our offset from the center of the parent's node to
            // the center of the child's node
            Vector_3 offset = boundsOffsetTable[i] * bounds.radius() ;

            // Create a new Bounds, with the center offset and half the radius
            CubicBounds<K> newBounds ;
            newBounds.setRadius( bounds.radius() * 0.5 ) ;
            newBounds.setCenter( bounds.center() + offset ) ;

            // Recurse
            //std::cout << "newList.size() = " << newList.size() << std::endl ;
            m_child[i]->build( newList, maxDepth, newBounds, currentDepth+1 ) ;
        }

        return true;
    }


    // Gets the mean point in each leaf. Do not forget to empty the meanPts vector before using this function!
    const bool getMeanPoints( std::vector<Point_3>& meanPts ) const
    {
        // If I'm a node, recursively traverse my children
        if (!m_pointCount)
        {
            for (unsigned int i = 0; i < 8; i++)
            {
                // We store incomplete trees (i.e. we're not guaranteed
                // that a node has all 8 children)
                if (!m_child[i]) continue;

                m_child[i]->getMeanPoints( meanPts ) ;
            }
        }
        else {
            // We are in a leaf, do stuff
            FT mx = 0, my = 0, mz = 0 ;
            typename std::vector<Point_3>::const_iterator it ;
            for ( it = m_points.begin(); it != m_points.end(); ++it ) {
                mx = mx + it->x() ;
                my = my + it->y() ;
                mz = mz + it->z() ;
            }
            int numPts = (int)m_points.size() ;
            mx = mx / numPts ;
            my = my / numPts ;
            mz = mz / numPts ;

            Point_3 meanPt( mx, my, mz ) ;
            meanPts.push_back( meanPt ) ;
        }

        return true;
    }


    // TODO: Generic Octree traversal
//    const bool Octree::traverse(callback proc, void *data) const
//    {
//            // Call the callback for this node (if the callback returns false, then
//            // stop traversing.

//            if (!proc(*this, data)) return false;

//            // If I'm a node, recursively traverse my children

//            if (!_pointCount)
//            {
//                    for (unsigned int i = 0; i < 8; i++)
//                    {
//                            // We store incomplete trees (i.e. we're not guaranteed
//                            // that a node has all 8 children)

//                            if (!_child[i]) continue;

//                            if (!_child[i]->traverse(proc, data)) return false;
//                    }
//            }

//            return true;
//    }

private:
    CGAL_Octree<K>          *m_child[8] ;
    unsigned int            m_pointCount ;
    std::vector< Point_3 >  m_points ;
    Point_3                 m_center ;
    double                  m_radius ;

};

}
