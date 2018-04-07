/***********************************************************

   Starter code for Assignment 3

   Implements scene_object.h

***********************************************************/

#include <cmath>
#include "scene_object.h"
#include <iostream>
using namespace std;

bool UnitSquare::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
                           const Matrix4x4& modelToWorld) {
        // TODO: implement intersection code for UnitSquare, which is
        // defined on the xy-plane, with vertices (0.5, 0.5, 0),
        // (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
        // (0, 0, 1).
        //
        // Your goal here is to fill ray.intersection with correct values
        // should an intersection occur.  This includes intersection.point,
        // intersection.normal, intersection.none, intersection.t_value.
        //
        // HINT: Remember to first transform the ray into object space
        // to simplify the intersection test.

        Ray3D rayInObject;
        rayInObject.origin = worldToModel * ray.origin;
        rayInObject.dir = (worldToModel * ray.dir);

        Vector3D normal = Vector3D(0.0,0.0,1.0);
        Point3D centre = Point3D(0.0,0.0,0.0);

        // t = (l0 - p0)n / (ln)
        Vector3D L = centre - rayInObject.origin;

        if (rayInObject.dir.dot(normal) > 1e-6) {
                return false;
        }

        float t = L.dot(normal)/(rayInObject.dir.dot(normal));
        Point3D intersectP = rayInObject.origin + t * rayInObject.dir;

        if (intersectP[0] < 0.5 && intersectP[0] > -0.5 && intersectP[1] < 0.5 && intersectP[1] > -0.5) {
                if (ray.intersection.none || t < ray.intersection.t_value) {
                        ray.intersection.point = modelToWorld * intersectP;

                        ray.intersection.normal = transNorm(worldToModel, normal);
                        ray.intersection.normal.normalize();

                        ray.intersection.none = false;
                        ray.intersection.t_value = t;

                        //-------------- UV coordinate ------------//
                        ray.intersection.u = intersectP[0] + 0.5;
                        ray.intersection.v = 0.5 - intersectP[1];

                        //std::cout << ray.intersection.u << " and " << ray.intersection.v << '\n';

                        return true;
                }
        }

        return false;
}

bool UnitSphere::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
                           const Matrix4x4& modelToWorld) {
        // TODO: implement intersection code for UnitSphere, which is centred
        // on the origin.
        //
        // Your goal here is to fill ray.intersection with correct values
        // should an intersection occur.  This includes intersection.point,
        // intersection.normal, intersection.none, intersection.t_value.
        //
        // HINT: Remember to first transform the ray into object space
        // to simplify the intersection test.

        Ray3D rayInObject;
        rayInObject.origin = worldToModel * ray.origin;
        rayInObject.dir = worldToModel * ray.dir;

        Point3D centre = Point3D(0,0,0);
        Vector3D L = rayInObject.origin - centre;

        // analytic solution, use implitic function
        float t0, t1, t;
        float a = rayInObject.dir.dot(rayInObject.dir);
        float b = 2 * rayInObject.dir.dot(L);
        float c = L.dot(L) - 1;

        float delta = b * b - 4 * a * c;

        // calulate the intersection points
        if (delta < 0) {
                return false;
        } else if (delta == 0) {
                t0 = -0.5 * b / a;
                t1 = t0;
        } else {
                t0 = 0.5 * (-b + sqrt(delta)) / a;
                t1 = 0.5 * (-b - sqrt(delta)) / a;
        }

        if (t0 > t1) swap(t0, t1);

        // check if intersection point is valid
        if (t0 > 0) {
                t = t0;
                // cout << "Solutions " << t0 << endl;
        } else if (t0 < 0 && t1 > 1) {
                t = t1;
        } else {
                return false;
        }

        // fill Intersection
        if (ray.intersection.none || t < ray.intersection.t_value) {
                Point3D intersectP = (rayInObject.origin + (t * rayInObject.dir));
                ray.intersection.point = modelToWorld * intersectP;

                Vector3D normal = intersectP - centre;
                ray.intersection.normal = transNorm(worldToModel, normal);
                ray.intersection.normal.normalize();
                ray.intersection.none = false;
                ray.intersection.t_value = t;

                //------------ UV coordinate -------------------//
                Vector3D n = ray.intersection.normal;
                ray.intersection.u = (atan2(n[0], n[2])/(2*M_PI)) + 0.5;
                ray.intersection.v = n[1] * 0.5 + 0.5;

                return true;
        }

        return false;
}

// ---------------------- Cube intersection ---------------------------- //
bool UnitCube::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
                         const Matrix4x4& modelToWorld) {
        // TODO: implement intersection code for UnitCube, which is
        // defined on the xyz-plane, with vertices (0, 0, 0),
        // (1, 0, 0), (1, 1, 0), (1, 0, 1), (0,1,0), (0,1,1), (1,1,1), (0,0,1)
        //
        // Your goal here is to fill ray.intersection with correct values
        // should an intersection occur.  This includes intersection.point,
        // intersection.normal, intersection.none, intersection.t_value.

        Ray3D rayInObject;
        rayInObject.origin = worldToModel * ray.origin;
        rayInObject.dir = (worldToModel * ray.dir);

        double tmin = (0 - rayInObject.origin[0])/rayInObject.dir[0];
        double tmax = (1 - rayInObject.origin[0])/rayInObject.dir[0];

        if (tmin > tmax) swap(tmin, tmax);

        double tymin = (0 - rayInObject.origin[1])/rayInObject.dir[1];
        double tymax = (1 - rayInObject.origin[1])/rayInObject.dir[1];

        if (tymin > tymax) swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax)) {
                return false;
        }

        if (tymin > tmin) {
                tmin = tymin;
        }

        if (tymax < tmax) {
                tmax = tymax;
        }

        double tzmin = (0 - rayInObject.origin[2])/rayInObject.dir[2];
        double tzmax = (1 - rayInObject.origin[2])/rayInObject.dir[2];

        if (tzmin > tzmax) swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax)) {
                return false;
        }

        if (tzmin > tmin) {
                tmin = tzmin;
        }

        if (tzmax < tmax) {
                tmax = tzmax;
        }

        double t = tmin;
        if (t < 0) {
                return false;
        }
        if (ray.intersection.none || t < ray.intersection.t_value) {

                Point3D intersectP = (rayInObject.origin + (t * rayInObject.dir));
                ray.intersection.point = modelToWorld * intersectP;

                Vector3D normal;
                if ((1 - intersectP[0]) < 1e-15) {
                        normal = Vector3D(1.0,0.0,0.0);
                        ray.intersection.u = intersectP[1]/4.0 + 0.25;
                        ray.intersection.v = 0.34 - intersectP[2]/3.0;

                } else if ((0-intersectP[0]) < 1e-15) {
                        normal = Vector3D(-1.0,0.0,0.0);
                        ray.intersection.u = intersectP[1]/4.0 + 0.25;
                        ray.intersection.v = intersectP[2]/3.0 + 0.66;

                } else if ((1 - intersectP[1]) < 1e-15) {
                        normal = Vector3D(0.0,1.0,0.0);
                        ray.intersection.u = 0.5 + intersectP[2]/4.0;
                        ray.intersection.v = 0.67 - intersectP[0]/3.0;
                } else if ((0-intersectP[1]) < 1e-15) {
                        normal = Vector3D(0.0,-1.0,0.0);
                        ray.intersection.u = 0.25 - intersectP[2]/4.0;
                        ray.intersection.v = 0.67 - intersectP[0]/3.0;
                } else if ((1 - intersectP[2]) < 1e-15) {
                        normal = Vector3D(0.0,0.0,1.0);
                        ray.intersection.u = 1.0 - intersectP[1]/4.0;
                        ray.intersection.v = 0.67 - intersectP[0]/3.0;
                } else if ((0-intersectP[2]) < 1e-15) {
                        normal = Vector3D(0.0,0.0,-1.0);
                        ray.intersection.u = 0.25 + intersectP[1]/4.0;
                        ray.intersection.v = 0.67 - intersectP[0]/3.0;

                } else {
                        std::cout << intersectP << '\n';
                        std::cout << "Error: wrong normal on cube" << '\n';
                }


                ray.intersection.normal = transNorm(worldToModel, normal);
                //ray.intersection.normal.normalize();
                ray.intersection.none = false;
                ray.intersection.t_value = t;


                //------------ UV coordinate -------------------//
                return true;
        }

        return false;
}

// ------------------- Cube ends her ----------------------- //

void SceneNode::rotate(char axis, double angle) {
        Matrix4x4 rotation;
        double toRadian = 2*M_PI/360.0;
        int i;

        for (i = 0; i < 2; i++) {
                switch(axis) {
                case 'x':
                        rotation[0][0] = 1;
                        rotation[1][1] = cos(angle*toRadian);
                        rotation[1][2] = -sin(angle*toRadian);
                        rotation[2][1] = sin(angle*toRadian);
                        rotation[2][2] = cos(angle*toRadian);
                        rotation[3][3] = 1;
                        break;
                case 'y':
                        rotation[0][0] = cos(angle*toRadian);
                        rotation[0][2] = sin(angle*toRadian);
                        rotation[1][1] = 1;
                        rotation[2][0] = -sin(angle*toRadian);
                        rotation[2][2] = cos(angle*toRadian);
                        rotation[3][3] = 1;
                        break;
                case 'z':
                        rotation[0][0] = cos(angle*toRadian);
                        rotation[0][1] = -sin(angle*toRadian);
                        rotation[1][0] = sin(angle*toRadian);
                        rotation[1][1] = cos(angle*toRadian);
                        rotation[2][2] = 1;
                        rotation[3][3] = 1;
                        break;
                }
                if (i == 0) {
                        this->trans = this->trans*rotation;
                        angle = -angle;
                }
                else {
                        this->invtrans = rotation*this->invtrans;
                }
        }
}

void SceneNode::translate(Vector3D trans) {
        Matrix4x4 translation;

        translation[0][3] = trans[0];
        translation[1][3] = trans[1];
        translation[2][3] = trans[2];
        this->trans = this->trans*translation;
        translation[0][3] = -trans[0];
        translation[1][3] = -trans[1];
        translation[2][3] = -trans[2];
        this->invtrans = translation*this->invtrans;
}

void SceneNode::scale(Point3D origin, double factor[3] ) {
        Matrix4x4 scale;

        scale[0][0] = factor[0];
        scale[0][3] = origin[0] - factor[0] * origin[0];
        scale[1][1] = factor[1];
        scale[1][3] = origin[1] - factor[1] * origin[1];
        scale[2][2] = factor[2];
        scale[2][3] = origin[2] - factor[2] * origin[2];
        this->trans = this->trans*scale;
        scale[0][0] = 1/factor[0];
        scale[0][3] = origin[0] - 1/factor[0] * origin[0];
        scale[1][1] = 1/factor[1];
        scale[1][3] = origin[1] - 1/factor[1] * origin[1];
        scale[2][2] = 1/factor[2];
        scale[2][3] = origin[2] - 1/factor[2] * origin[2];
        this->invtrans = scale*this->invtrans;
}
