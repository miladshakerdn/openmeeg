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

#pragma once

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <string>

#include <mesh.h>
#include <geometry.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>

#include <CGAL/config.h>
#include <CGAL/version.h>

namespace OpenMEEG {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef K::Point_3 Point_3;
    typedef K::FT FT;
    typedef Polyhedron::HalfedgeDS HDS;

    template <typename CGAL_MESH>
    Mesh CGAL_to_OM(const CGAL_MESH& mesh) {

        typedef typename CGAL_MESH::Triangulation Triangulation;
        typedef typename Triangulation::Point     Point;
        typedef typename CGAL_MESH::Vertex_handle VertexHandle;
        typedef typename CGAL_MESH::Facets_in_complex_iterator ExternalFacetsIterator;

        Mesh m(mesh.triangulation().number_of_vertices(),mesh.number_of_facets());

        std::set<VertexHandle> vertex_set;
        for (ExternalFacetsIterator fit=mesh.facets_in_complex_begin(); fit!=mesh.facets_in_complex_end(); ++fit) {
            unsigned j = 0;
            for (unsigned i=0; i<4; ++i)
                if (i!=fit->second) {
                    vertex_set.insert(fit->first->vertex(i));
                    j += 1;
                }
            if (j!=3)
                throw "Wrong complex";
        }

        Vertices vertices;
        std::map<VertexHandle,unsigned> mapping;
        unsigned index = 0;
        for (const auto& vit : vertex_set) {
            const Point& point = vit->point();
            const Vertex v(CGAL::to_double(point.x()),CGAL::to_double(point.y()),CGAL::to_double(point.z()));
            vertices.push_back(v);
            mapping[vit] = index++;
        }
        m.geometry().add_vertices(vertices);

        for (ExternalFacetsIterator fit=mesh.facets_in_complex_begin(); fit!=mesh.facets_in_complex_end(); ++fit) {
            TriangleIndices indices;
            for (unsigned i=0,j=0; i<4; ++i)
                if (i!=fit->second)
                    indices[j++] = mapping[fit->first->vertex(i)];
            m.add_triangle(indices);
        }

        m.update(true);
        m.correct_global_orientation();
        return m;
    }

    Mesh cgal_mesh_function(const double sphere_radius,const double hemisphere,const double radius_bound,const double distance_bound);

    #if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,6,0)
    Mesh cgal_refine(const Mesh& m_in,const double radius_bound,const double distance_bound,const std::string& sizing_field);

    Mesh cgal_decimate(const Mesh& m_in,const unsigned nb_points);
    #endif

    #ifdef CGAL_ImageIO
    Mesh cgal_mesh_3Dlabeled_image(const std::string& input_filename,const double radius_bound,const double distance_bound);
    Mesh cgal_mesh_3Dlevelset_image(const std::string& input_filename,const double levelset_value,const bool positive_inside,const double radius_bound,const double distance_bound);
    #endif
}
