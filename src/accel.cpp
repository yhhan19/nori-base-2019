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

#include <nori/accel.h>
#include <nori/timer.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

int innerN = 0, leafN = 0, triN = 0;

Accel::node *Accel::alloc() {
	node *x = new node;
	for (int i = 0; i < 8; i++) 
		x->child[i] = nullptr;
	x->len = 10;
	x->index = (uint32_t *)malloc(sizeof(uint32_t) * x->len);
	for (uint32_t i = 0; i < x->len; i++)
		x->index[i] = 0;
	x->t = 0;
	return x;
}

void Accel::getBoundingBoxes(BoundingBox3f *b, BoundingBox3f bbox) const {
	float L[3], R[3], M[3];
	for (int i = 0; i < 3; i++) {
		L[i] = bbox.min[i];
		R[i] = bbox.max[i];
		M[i] = (bbox.min[i] + bbox.max[i]) / 2;
	}
	b[0] = BoundingBox3f(Point3f(L[0], L[1], L[2]), Point3f(M[0], M[1], M[2]));
	b[1] = BoundingBox3f(Point3f(M[0], L[1], L[2]), Point3f(R[0], M[1], M[2]));
	b[2] = BoundingBox3f(Point3f(L[0], M[1], L[2]), Point3f(M[0], R[1], M[2]));
	b[3] = BoundingBox3f(Point3f(M[0], M[1], L[2]), Point3f(R[0], R[1], M[2]));
	b[4] = BoundingBox3f(Point3f(L[0], L[1], M[2]), Point3f(M[0], M[1], R[2]));
	b[5] = BoundingBox3f(Point3f(M[0], L[1], M[2]), Point3f(R[0], M[1], R[2]));
	b[6] = BoundingBox3f(Point3f(L[0], M[1], M[2]), Point3f(M[0], R[1], R[2]));
	b[7] = BoundingBox3f(Point3f(M[0], M[1], M[2]), Point3f(R[0], R[1], R[2]));
}

void Accel::insert(node *x, BoundingBox3f bbox, uint32_t index, int depth) {
	if (depth == 16) {
		if (x->t == x->len) {
			x->len *= 2;
			uint32_t *temp = (uint32_t *)malloc(sizeof(uint32_t) * x->len);
			for (uint32_t i = 0; i < x->t; i++) {
				temp[i] = x->index[i];
			}
			delete x->index;
			x->index = temp;
		}
		x->index[x->t++] = index;
	}
	else {
		if (x->child[0] == nullptr) {
			if (x->t < x->len) {
				x->index[x->t++] = index;
			}
			else {
				BoundingBox3f b[8]; getBoundingBoxes(b, bbox);
				for (int i = 0; i < 8; i++)
					x->child[i] = alloc();
				for (uint32_t i = 0; i < x->t; i++)
					for (int j = 0; j < 8; j++) 
						if (m_mesh->boundingBoxIntersect(x->index[i], b[j])) 
							insert(x->child[j], b[j], x->index[i], depth + 1);
				for (int i = 0; i < 8; i++) 
					if (m_mesh->boundingBoxIntersect(index, b[i]))
						insert(x->child[i], b[i], index, depth + 1);
			}
		}
		else {
			BoundingBox3f b[8]; getBoundingBoxes(b, bbox);
			for (int i = 0; i < 8; i++) 
				if (m_mesh->boundingBoxIntersect(index, b[i]))
					insert(x->child[i], b[i], index, depth + 1);
		}
	}
}

void Accel::stat(node *x) {
	if (x->child[0] == nullptr) {
		leafN++;
		triN += x->t;
	}
	else {
		innerN++;
		for (int i = 0; i < 8; i++)
			stat(x->child[i]);
	}
}

void Accel::build() {
	Timer timer;
	root = alloc();
	for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) 
		insert(root, m_bbox, idx, 0);
	stat(root);
	cout << "octree node size: " << sizeof(node) << "." << endl;
	cout << leafN << " leaves, " << innerN << " inner nodes, " << (float)triN / leafN << " triangles per leaf." << endl;
	cout << "done. (took " << timer.elapsedString() << ")" << endl;
}

bool Accel::search(node *x, BoundingBox3f bbox, Ray3f &ray, Intersection &its, bool shadowRay, uint32_t &f) const {
	bool ret = false;
	if (x->child[0] == nullptr) {
		for (uint32_t i = 0; i < x->t; i++) {
			float u, v, t;
			if (m_mesh->rayIntersect(x->index[i], ray, u, v, t)) {
				if (shadowRay) return true;
				ray.maxt = its.t = t;
				its.uv = Point2f(u, v);
				its.mesh = m_mesh;
				f = x->index[i];
				ret = true;
			}
		}
	}
	else {
		float near[8], far;
		int index[8], t = 0;
		BoundingBox3f b[8]; getBoundingBoxes(b, bbox);
		for (int i = 0; i < 8; i++) {
			if (b[i].rayIntersect(ray, near[i], far)) 
				if (ray.mint <= far && near[i] <= ray.maxt)
					index[t++] = i;
		}
		for (int i = 0; i < t; i ++) 
			for (int j = i + 1; j < t; j ++)
				if (near[index[i]] > near[index[j]]) {
					int temp = index[i];
					index[i] = index[j];
					index[j] = temp;
				}
		for (int i = 0; i < t; i ++) {
			if (near[index[i]] <= ray.maxt) {
				bool res = search(x->child[index[i]], b[index[i]], ray, its, shadowRay, f);
				if (shadowRay && res) return true;
				ret = ret || res;
			}
		}
	}
	return ret;
}

bool Accel::search_(node *x, BoundingBox3f bbox, Ray3f &ray, Intersection &its, bool shadowRay, uint32_t &f) const {
	bool ret = false;
	if (x->child[0] == nullptr) {
		for (uint32_t i = 0; i < x->t; i++) {
			float u, v, t;
			if (m_mesh->rayIntersect(x->index[i], ray, u, v, t)) {
				if (shadowRay) return true;
				ray.maxt = its.t = t;
				its.uv = Point2f(u, v);
				its.mesh = m_mesh;
				f = x->index[i];
				ret = true;
			}
		}
	}
	else {
		BoundingBox3f b[8]; getBoundingBoxes(b, bbox);
		for (int i = 0; i < 8; i++) {
			if (b[i].rayIntersect(ray)) {
				bool res = search_(x->child[i], b[i], ray, its, shadowRay, f);
				if (shadowRay && res) return true;
				ret = ret || res;
			}
		}
	}
	return ret;
}

bool Accel::bf(Ray3f &ray, Intersection &its, bool shadowRay, uint32_t &f) const {
	bool foundIntersection = false;
	/* Brute force search through all triangles */
	for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
		float u, v, t;
		if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
			/* An intersection was found! Can terminate
			   immediately if this is a shadow ray query */
			if (shadowRay)
				return true;
			ray.maxt = its.t = t;
			its.uv = Point2f(u, v);
			its.mesh = m_mesh;
			f = idx;
			foundIntersection = true;
		}
	}
	return foundIntersection;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
	//foundIntersection = bf(ray, its, shadowRay, f);
	//foundIntersection = search_(root, m_bbox, ray, its, shadowRay, f);
	foundIntersection = search(root, m_bbox, ray, its, shadowRay, f);

	if (shadowRay && foundIntersection)
		return true;

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

