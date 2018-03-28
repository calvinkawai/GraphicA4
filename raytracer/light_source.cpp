/***********************************************************

	Starter code for Assignment 3

	Implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"
#include <iostream>
using namespace std;

void PointLight::shade(Ray3D& ray) {
	// TODO: implement this function to fill in values for ray.col
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray
	// is available.  So be sure that traverseScene() is called on the ray
	// before this function.

  Vector3D N = ray.intersection.normal;
  Vector3D V = -ray.dir;
  Vector3D L = get_position() - ray.intersection.point;
  N.normalize();
  V.normalize();
  L.normalize();

  Vector3D R = 2.0 * (L.dot(N)) * N - L;
  R.normalize();

  // phong shading
  float NL = max(0.0,N.dot(L));
  float VR = max(0.0, pow(V.dot(R), ray.intersection.mat->specular_exp));

  Color ambient = ray.intersection.mat->ambient * col_ambient;
  Color diffuse = ray.intersection.mat->diffuse * (NL * col_diffuse);
  Color specular = Color(0.0, 0.0, 0.0);

  if (NL > 0) {
    specular = ray.intersection.mat->specular * (VR * col_specular);
  }

  ray.col = ambient + diffuse + specular;
  ray.col.clamp();
}
