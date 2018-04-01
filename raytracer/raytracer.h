/***********************************************************

	Starter code for Assignment 3

	This file contains the interface and datastructures of the raytracer.
	Simple traversal and addition code to the datastructures are given to you.

***********************************************************/
#pragma once

#include "util.h"
#include "scene_object.h"
#include "light_source.h"

class Raytracer {
public:
	// Renders 3D scene to an image given camera and lights setup.
	void render(Camera& camera, Scene& scene, LightList& light_list, Image& image);
    void renderWithoutAliasing(Camera& camera, Scene& scene, LightList& light_list, Image& image, int depth);
    void renderWithAntiAliasing(Camera& camera, Scene& scene, LightList& light_list, Image& image, int depth, int ray_num);

private:

	// Return the color of the ray after intersection and shading, call
	// this function recursively for reflection and refraction.
	Color shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int depth);

	// Traversal code for the scene, the ray is transformed into
	// the object space of each node where intersection is performed.
	void traverseScene(Scene& scene, Ray3D& ray);
    
	// After intersection, calculate the color of the ray by shading it
	// with all light sources in the scene.
	void computeShading(Ray3D& ray, LightList& light_list);
    
    // -------- additional ---------
    // METHOD0: hard shadow
    Color hardShadowing(Ray3D& ray, Scene& scene, LightList& light_list);
    
    // METHOD1: Soft shadowing using light area as grid
    // uniform distribution through all the samples
    Color softShadowingGrid(Ray3D& ray, Scene& scene, LightList& light_list, int grid_size);
    
    
    //METHOD2: Soft shadowing using spherical light area
    Color softShadowingSpherical(Ray3D& ray, Scene& scene, LightList& light_list, double radius);
    
    // --------------------------

	// Precompute the modelToWorld and worldToModel transformations for each
    // object in the scene.
	void computeTransforms(Scene& scene);

};
