#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
	Point3f p;
	Color3f e;

	SimpleIntegrator(const PropertyList &props) {
		p = props.getPoint("position");
		e = props.getColor("energy");
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);
		Normal3f n = its.shFrame.n;
		Vector3f v = p - its.p;
		if (scene->rayIntersect(Ray3f(its.p, v)) || v.dot(n) < 0) 
			return Color3f(0.0f);
		float C = v.dot(n) / sqrt(v.dot(v)) / v.dot(v) / (4 * M_PI * M_PI);
		return  Color3f(e.r() * C, e.g() * C, e.b() * C);
	}

	std::string toString() const {
		return "SimpleIntegrator[]";
	}
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
