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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
	dpdf = DiscretePDF(0);
	s = 0;
	for (uint32_t idx = 0; idx < getTriangleCount(); ++idx) {
		dpdf.append(surfaceArea(idx));
		s += surfaceArea(idx);
	}
	dpdf.normalize();
	int N = 5;
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			for (int k = 1; k < N; k++) {
				float fi = i / (float)N, fj = j / (float)N, fk = k / (float)N;
				Mesh::samplePoint sp = sample(Point3f(fi, fj, fk));
				cout << sp.p.toString() << endl;
			}
			cout << endl;
		}
		cout << endl;
	}
}

Mesh::samplePoint Mesh::sampleTriangle(int index, Point2f p) {
	uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
	const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
	float a = 1 - sqrt(1 - p.coeff(0)), b = p.coeff(1) * sqrt(1 - p.coeff(0)), c = 1 - a - b;
	samplePoint sp;
	float x = p0.x() * a + p1.x() * b + p2.x() * c, y = p0.y() * a + p1.y() * b + p2.y() * c, z = p0.z() * a + p1.z() * b + p2.z() * c;
	sp.p = Point3f(x, y, z);
	Vector3f v = (p1 - p0).cross(p2 - p0);
	float n = v.norm();
	sp.n = Normal3f(v.x() / n, v.y() / n, v.z() / n);
	sp.pdf = surfaceArea(index) / s;
	return sp;
}

Mesh::samplePoint Mesh::sample(Point3f p) {
	uint32_t index = dpdf.sample(p.coeff(0));
	return sampleTriangle(index, Point2f(p.coeff(1), p.coeff(2)));
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}


bool Mesh::axis(Vector3f ex, const Vector3f *u, const Vector3f *v, Vector3f a) const {
	float r = 0, min = 0, max = 0;
	for (int k = 0; k < 3; k++)
		r += ex.dot(u[k]) * fabs(u[k].dot(a));
	for (int k = 0; k < 3; k++) {
		float proj = v[k].dot(a);
		if (k == 0 || proj > max) max = proj;
		if (k == 0 || proj < min) min = proj;
	}
	if (max < -r || min > r) return false;
	return true;
}

bool Mesh::boundingBoxIntersect(uint32_t index, const BoundingBox3f &box) const {
	uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
	const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
	Vector3f c = (box.min + box.max) / 2;
	Vector3f ex = (box.max - box.min) / 2;
	Vector3f v[3], e[3], u[3];
	v[0] = p0 - c; v[1] = p1 - c; v[2] = p2 - c;
	e[0] = v[1] - v[0]; e[1] = v[2] - v[1]; e[2] = v[0] - v[2];
	u[0] = Vector3f(1.0f, 0.0f, 0.0f);
	u[1] = Vector3f(0.0f, 1.0f, 0.0f);
	u[2] = Vector3f(0.0f, 0.0f, 1.0f);
	for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++) 
			if (!axis(ex, u, v, u[i].cross(e[j])))
				return false;
	for (int i = 0; i < 3; i++) 
		if (!axis(ex, u, v, u[i]))
			return false;
	if (!axis(ex, u, v, e[0].cross(e[1])))
		return false;
	return true;
 }

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
