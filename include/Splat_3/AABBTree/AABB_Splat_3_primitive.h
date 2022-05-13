
#ifndef AABB_SPLAT_3_PRIMITIVE_H
#define AABB_SPLAT_3_PRIMITIVE_H

#include "Splat_3/Splat_3.h"

namespace CGAL {

    template <class GeomTraits, class Iterator>
    class AABB_Splat_3_primitive
    {

    public:
        
		// types
        typedef Iterator Id; // Id type
        typedef typename GeomTraits::Point_3 Point; // point type
        typedef typename CGAL::Splat_3<GeomTraits> Datum; // datum type

    private:
        
		// member data
        Id m_it; // iterator
        Datum m_datum; // a Splat_3
		        
    public:
		
		// constructor
        AABB_Splat_3_primitive() {}
        
		AABB_Splat_3_primitive(Id it)
			: m_it(it)
        {
            // m_datum = *it ; // copy splat
			m_datum = Datum( *it ) ; // copy splat
        }

        AABB_Splat_3_primitive(const AABB_Splat_3_primitive& primitive) {
			m_datum = primitive.datum() ;
            m_it = primitive.id() ;
        }

    public:
        Id& id() { return m_it ; }
        const Id& id() const { return m_it ; }
        Datum& datum() { return m_datum ; }
        const Datum& datum() const { return m_datum ; }

        /// Returns a point on the primitive
        Point reference_point() const { return m_datum.center() ; }
    };

}  // end namespace CGAL


#endif // AABB_SPLAT_3_PRIMITIVE_H

