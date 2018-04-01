/***********************************************************

	Starter code for Assignment 3

	Implementations of functions in raytracer.h,
	and the main function which specifies the scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

void Raytracer::traverseScene(Scene& scene, Ray3D& ray)  {
	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
}

void Raytracer::computeTransforms(Scene& scene) {
	// right now this method might seem redundant. But if you decide to implement
	// scene graph this is where you would propagate transformations to child nodes

	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		node->modelToWorld = node->trans;
		node->worldToModel = node->invtrans;
	}
}

void Raytracer::computeShading(Ray3D& ray, LightList& light_list) {
	for (size_t  i = 0; i < light_list.size(); ++i) {
		LightSource* light = light_list[i];
        
		// Each lightSource provides its own shading function.
		// Implement shadows here if needed.
        light->shade(ray);
	}
}



// METHOD0: hard shadow
Color Raytracer::hardShadowing(Ray3D& ray, Scene& scene, LightList& light_list){
    
    Color col(0.0, 0.0, 0.0);
    
    if (!ray.intersection.none) {
        
        for (size_t index = 0; index < light_list.size(); ++index) {
            LightSource* light = light_list[index];
            Ray3D shadowRay;
            shadowRay.dir = light->get_position() - ray.intersection.point;
            shadowRay.dir.normalize();
            shadowRay.origin = ray.intersection.point + 0.003*shadowRay.dir;
            
            
            traverseScene(scene, shadowRay);
            
            if (shadowRay.intersection.none){
                computeShading(ray, light_list);
                col = col + ray.col;
            }
        }
    }
    return col = (1/ (double) light_list.size())*col;
}

// METHOD1: Soft shadowing using light area as grid
// uniform distribution through all the samples
Color Raytracer::softShadowingGrid(Ray3D& ray, Scene& scene, LightList& light_list, int grid_size){
    
    Color col(0.0, 0.0, 0.0);
    
    if (!ray.intersection.none) {
        
        for (size_t index = 0; index < light_list.size(); ++index) {
            LightSource* light = light_list[index];
            
            double ix=- grid_size*0.5;
            
            while (ix<=grid_size*0.5){
                
                double iy = -grid_size*0.5;
                while (iy<=grid_size*0.5){
                    
//                    double iz = -grid_size*0.5;
//                    while (iz <=grid_size*0.5){
                        Ray3D shadowRay;
                        Point3D light_pos = light->get_position() + Vector3D(ix, iy, 0);
                        
                        shadowRay.dir = light_pos - ray.intersection.point;
                        shadowRay.dir.normalize();
                        shadowRay.origin = ray.intersection.point + 0.002*shadowRay.dir;
                        
                        traverseScene(scene, shadowRay);
                        
                        if (shadowRay.intersection.none){
                            computeShading(ray, light_list);
                            col = col + ray.col;
                        }
                        
//                        iz= iz+0.5;
//                    }
                    iy= iy+0.5;
                }
                ix = ix+0.5;
            }
        }
    }
    return col = (1/ ((double) light_list.size()*(2*grid_size+1)*(2*grid_size+1)))*col;
}


//METHOD2: Soft shadowing using spherical light area
Color Raytracer::softShadowingSpherical(Ray3D& ray, Scene& scene, LightList& light_list, int radius){
    
    Color col(0.0, 0.0, 0.0);
    
    if (!ray.intersection.none) {
    
        // iterate through all the light source
        for (size_t index = 0; index < light_list.size(); ++index) {
            LightSource* light = light_list[index];
            
            double r = radius * (rand()/((double) RAND_MAX)) ;
            double theta = 2 * M_PI * (rand()/((double) RAND_MAX));
            double phi = 2 * M_PI * (rand()/((double) RAND_MAX));
            
            Ray3D shadowRay;
            Point3D light_pos = light->get_position();
            Point3D random_light_pos = light_pos + Vector3D(r*cos(theta)*sin(phi), r*sin(theta)*sin(phi), r*cos(phi));
            shadowRay.dir = random_light_pos - ray.intersection.point;
            shadowRay.dir.normalize();
            shadowRay.origin = ray.intersection.point + 0.001*shadowRay.dir;
            
            traverseScene(scene, shadowRay);
            
            if (shadowRay.intersection.none){
                computeShading(ray, light_list);
                col = col + ray.col;
            }
            
        }
    }
    return col = (1.0/((double) light_list.size()))*col;
    
}




Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int depth) {
	Color col(0.0, 0.0, 0.0);
	traverseScene(scene, ray);

	// Don't bother shading if the ray didn't hit
	// anything.
	if (!ray.intersection.none) {
        col = hardShadowing(ray, scene, light_list);
//        col = softShadowingGrid(ray, scene, light_list, 1);

	}

	// You'll want to call shadeRay recursively (with a different ray,
	// of course) here to implement reflection/refraction effects.
    if (!ray.intersection.none && depth > 0) {

		// get the reflect ray
		Ray3D refRay;
		Vector3D L = -ray.dir;
	    Vector3D N = ray.intersection.normal;
		L.normalize();
		N.normalize();
		Vector3D R = 2.0 * (L.dot(N)) * N - L;
		R.normalize();


		refRay.origin = ray.intersection.point + 0.001*R;
		refRay.dir = R;

		Color refCol(0.0, 0.0, 0.0);
		refCol = refCol + shadeRay(refRay, scene, light_list, --depth);
		refRay.col = refCol;
		if (L.dot(N) > 1e-6) {
            col = col + 0.9 * ray.intersection.mat->specular * refCol;
		}
	}
	col.clamp();
	return col;
}



void Raytracer::render(Camera& camera, Scene& scene, LightList& light_list, Image& image) {
	computeTransforms(scene);

	Matrix4x4 viewToWorld;
	double factor = (double(image.height)/2)/tan(camera.fov*M_PI/360.0);

	viewToWorld = camera.initInvViewMatrix();

	// Construct a ray for each pixel.
	for (int i = 0; i < image.height; i++) {
		for (int j = 0; j < image.width; j++) {
			// Sets up ray origin and direction in view space,
			// image plane is at z = -1.
            
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			imagePlane[0] = (-double(image.width)/2 + 0.5 + j)/factor;
			imagePlane[1] = (-double(image.height)/2 + 0.5 + i)/factor;
			imagePlane[2] = -1;



			Ray3D ray;
			// TODO: Convert ray to world space
			Vector3D direction = imagePlane - origin;

			direction = viewToWorld * direction;
			origin = viewToWorld * origin;
			ray = Ray3D(origin, direction);

			Color col = shadeRay(ray, scene, light_list, 2);
			image.setColorAtPixel(i, j, col);
		}
	}
}
