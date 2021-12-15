/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include <cgal_lib.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Gray_level_image_3.h>

namespace OpenMEEG {

    typedef CGAL::Labeled_mesh_domain_3<K> ImplicitDomain;
    typedef CGAL::Mesh_triangulation_3<ImplicitDomain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    // Criteria
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
    typedef CGAL::Gray_level_image_3<Geom_traits::FT, Geom_traits::Point_3> Gray_level_image;

    // To avoid verbose function and named parameters call
    using namespace CGAL::parameters;

    Mesh cgal_mesh_3Dlabeled_image(const std::string& input_filename,const double radius_bound,const double distance_bound) {
        // Mesh criteria
        CGAL::Image_3 image;
        image.read(input_filename.c_str());
        std::cout << "Input image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << std::endl;

        // Domain
        ImplicitDomain domain = ImplicitDomain::create_labeled_image_mesh_domain(image);

        // Mesh criteria
        Mesh_criteria criteria(facet_angle=30,facet_size=radius_bound,facet_distance=distance_bound);

        // Meshing
        const C3t3& c3t3 = CGAL::make_mesh_3<C3t3>(domain,criteria,no_exude(),no_perturb());

        // Output
        return CGAL_to_OM(c3t3);
    }

    Mesh cgal_mesh_3Dlevelset_image(const std::string& input_filename,const double levelset_value,const bool positive_inside, const double radius_bound,const double distance_bound)
    {
        // Mesh criteria
        double value_outside = 1.;
        const Gray_level_image image(input_filename.c_str(),levelset_value,positive_inside,value_outside);
        std::cout << "Input INR image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << "\n Positive values are " << (positive_inside?"Inside":"Outside") << std::endl;

        // Domain
        // definition of the surface, with 10^-8 as relative precision
        ImplicitDomain domain = ImplicitDomain::create_gray_image_mesh_domain(image,levelset_value,1e-3);

        // Mesh criteria
        Mesh_criteria criteria(facet_angle=30,facet_size=radius_bound,facet_distance=distance_bound);

        // Meshing
        const C3t3& c3t3 = CGAL::make_mesh_3<C3t3>(domain,criteria,no_exude(),no_perturb());

        // Output
        return CGAL_to_OM(c3t3);
    }
}
