#include "ray.h"

scalar_t Ray::env_ior = 1.0;

Ray::Ray() {
	ior = env_ior;
	energy = 1.0;
	time = 0;
}

Ray::Ray(const Vector3 &origin, const Vector3 &dir) {
	this->origin = origin;
	this->dir = dir;
	ior = env_ior;
	energy = 1.0;
	time = 0;
}

// TODO: remove this hack
void Ray::transform(const Matrix4x4 &xform) {
	//Matrix4x4 inv_transp_xform = xform.inverse();
	//inv_transp_xform.transpose();
	
	/*origin.transform(xform);
	dir.transform((Matrix3x3)xform);*/

	Vector3 d = dir - origin;
	d.transform((Matrix3x3)xform);
	origin.transform(xform);
	
	dir = d + origin;
}

Ray Ray::transformed(const Matrix4x4 &xform) const {
	Ray foo = *this;
	foo.transform(xform);
	return foo;
}

void Ray::enter(scalar_t new_ior) {
	ior = new_ior;
	ior_stack.push(ior);
}

void Ray::leave() {
	if(ior_stack.empty()) {
		//std::cerr << "empty ior stack?\n";
		return;
	}
	ior_stack.pop();
	ior = ior_stack.empty() ? env_ior : ior_stack.top();
}
