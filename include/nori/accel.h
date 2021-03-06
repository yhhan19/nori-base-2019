/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
	typedef struct node {
		node *child[8];
		uint32_t *index, t, len;
	} node;
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

	node * alloc();

	void getBoundingBoxes(BoundingBox3f *b, BoundingBox3f bbox) const;

	void insert(node * x, BoundingBox3f bbox, uint32_t index, int depth);

	void stat(node * x);

	/// Build the acceleration data structure (currently a no-op)
    void build();

	bool search(node * x, BoundingBox3f bbox, Ray3f & ray, Intersection &its, bool shadowRay, uint32_t &f) const;

	bool search_(node * x, BoundingBox3f bbox, Ray3f & ray, Intersection & its, bool shadowRay, uint32_t & f) const;

	bool bf(Ray3f & ray, Intersection & its, bool shadowRay, uint32_t & f) const;

	//bool search_(node * x, BoundingBox3f bbox, Ray3f & ray, Intersection & its, bool shadowRay, uint32_t & f) const;

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
	node *root = nullptr;
};

NORI_NAMESPACE_END
