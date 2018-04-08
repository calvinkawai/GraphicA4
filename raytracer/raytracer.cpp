/***********************************************************

   Starter code for Assignment 3

   Implementations of functions in raytracer.h,
   and the main function which specifies the scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include <omp.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

void Raytracer::traverseScene(Scene& scene, Ray3D& ray)  {
        for (size_t i = 0; i < scene.size(); ++i) {
                SceneNode* node = scene[i];

                if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
                        ray.intersection.mat = node->mat;
                    
                        // assign roughness
                        ray.intersection.roughness = node->roughness;
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
        for (size_t i = 0; i < light_list.size(); ++i) {
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

                        // create ray from intersection to light position
                        Ray3D shadowRay;
                        shadowRay.dir = light->get_position() - ray.intersection.point;
                        shadowRay.dir.normalize();
                        shadowRay.origin = ray.intersection.point + 0.003*shadowRay.dir;

                        traverseScene(scene, shadowRay);

                        if (shadowRay.intersection.none) {
                                computeShading(ray, light_list);
                                col = col + ray.col;

                        }


                }
        }

        return col = (1/ (double) light_list.size())*col;
}


//METHOD2: Soft shadowing using spherical light area
Color Raytracer::softShadowingSpherical(Ray3D& ray, Scene& scene, LightList& light_list, double radius){

        // number of sampling for each intersection.
        int ray_num = 20;

        Color col(0.0, 0.0, 0.0);

        if (!ray.intersection.none) {

                // iterate through all the light source
                for (size_t index = 0; index < light_list.size(); ++index) {
                        LightSource* light = light_list[index];

                        // uniform distribution for all number of sampling
                        for (int num=0; num < ray_num; num++) {

                                // random sampling for spherical light area
                                double r = radius * (rand()/((double) RAND_MAX));
                                double theta = 2 * M_PI * (rand()/((double) RAND_MAX));
                                double phi = 2 * M_PI * (rand()/((double) RAND_MAX));

                                Ray3D shadowRay;
                                Point3D light_pos = light->get_position();
                                Point3D random_light_pos = light_pos \
                                                           + Vector3D(r*cos(theta)*sin(phi), r*sin(theta)*sin(phi), r*cos(phi));
                                shadowRay.dir = random_light_pos - ray.intersection.point;
                                shadowRay.dir.normalize();
                                shadowRay.origin = ray.intersection.point + 0.001*shadowRay.dir;

                                traverseScene(scene, shadowRay);

                                if (shadowRay.intersection.none) {
                                        computeShading(ray, light_list);
                                        col = col + ray.col;
                                }
                        }
                }
        }
        return col = (1.0/((double) light_list.size()*ray_num))*col;

}


Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int depth) {
        Color col(0.0, 0.0, 0.0);
        traverseScene(scene, ray);


        // Don't bother shading if the ray didn't hit
        // anything.
        if (!ray.intersection.none) {


                col = hardShadowing(ray, scene, light_list);

                if (ray.intersection.mat->hasTexture) {
                        Color textureCol = textureColor(ray);
                        col = col + textureCol;
                }


                // int radius = 2;
                // col = softShadowingSpherical(ray, scene, light_list, 2);

        }

        // You'll want to call shadeRay recursively (with a different ray,
        // of course) here to implement reflection/refraction effects.
        if (!ray.intersection.none && depth > 0) {

                // get the reflect ray
                Ray3D refRay;
                Vector3D L = -ray.dir;
                Vector3D N = ray.intersection.normal;
                //L.normalize();
                //N.normalize();
                Vector3D R = 2.0 * (L.dot(N)) * N - L;
                R.normalize();

                refRay.origin = ray.intersection.point + 0.001*R;
                refRay.dir = R;

                Color refCol(0.0, 0.0, 0.0);

                // average out all the random samples with glossy reflection
                for (int i=0; i< 20; i++){
                        // compute new direction for glossy reflection
                        refRay.dir = glossyRefDirection(N, R, ray.intersection.roughness);
                        refRay.origin = ray.intersection.point + 0.001*refRay.dir;
                        refCol = refCol + 1.0/20*shadeRay(refRay, scene, light_list, --depth);
                
                }
            
                refRay.col = refCol;

                col = col + 0.9 * ray.intersection.mat->specular * refCol;

                if (ray.intersection.mat->hasTexture) {

                        Color textureCol = textureColor(ray);
                        col = col + textureCol;
                }
        }
        col.clamp();
        return col;
}


void Raytracer::render(Camera& camera, Scene& scene, LightList& light_list, Image& image) {

        // for shading
        int depth = 3;

        // for anti-aliasing
        int AA_num = 4;
        bool DOF = false;

        //renderWithoutAliasing(camera, scene, light_list, image, depth, DOF);
        renderWithAntiAliasing(camera, scene, light_list, image, depth, AA_num, DOF);

}

// rendering without anti-aliasing
// normal rendering
void Raytracer::renderWithoutAliasing(Camera& camera, Scene& scene, LightList& light_list, Image& image, int depth, bool DOF){
        computeTransforms(scene);

        Matrix4x4 viewToWorld;
        double factor = (double(image.height)/2)/tan(camera.fov*M_PI/360.0);

        viewToWorld = camera.initInvViewMatrix();

        // Construct a ray for each pixel.
        std::cout << "Rendering without using aliasing ... \n";
        #pragma omp parallel for
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

                        Color col(0.0, 0.0, 0.0);

                        // render with depth of field
                        if (DOF) {
                                double focal_len = 5.0;
                                double apeture_size = 1.0;
                                // ray_num samples for calculating DOF
                                int dof_num = 20;
                                Color dofCol(0.0, 0.0, 0.0);
                                for (int num = 0; num < dof_num; num++) {
                                        Ray3D secondary_ray = depthOfField(ray, origin, focal_len, apeture_size);
                                        col = col + shadeRay(secondary_ray, scene, light_list, depth);
                                }
                                col = double(1.0/dof_num)*col;

                        }else{
                                col = shadeRay(ray, scene, light_list, depth);
                        }

                        image.setColorAtPixel(i, j, col);

                }
        }

}

//render with antialising
void Raytracer::renderWithAntiAliasing(Camera& camera, Scene& scene, LightList& light_list, Image& image, int depth, int AA_num, bool DOF) {
        computeTransforms(scene);

        Matrix4x4 viewToWorld;
        double factor = (double(image.height)/2)/tan(camera.fov*M_PI/360.0);

        viewToWorld = camera.initInvViewMatrix();

        // Construct a ray for each pixel.
        std::cout << "Rendering using aliasing ... \n";
        #pragma omp parallel for
        for (int i = 0; i < image.height; i++) {
                for (int j = 0; j < image.width; j++) {
                        // Sets up ray origin and direction in view space,
                        // image plane is at z = -1.

                        Color col(0.0, 0.0, 0.0);

                        // 4 super pixels
                        for (int dx_i=0; dx_i<=1; dx_i++) {
                                for(int dy_i=0; dy_i<=1; dy_i++) {

                                        //x -- {-0.5, 0.5, 1.5}
                                        double dx = -0.5 + dx_i;
                                        double dy = -0.5 + dy_i;

                                        // averaging all the samples for anti-aliasing
                                        // ray_num numples for each super pixels
                                        for (int num=0; num< AA_num; num++) {

                                                // random position at the subpixel for x and y
                                                double sub_x = rand()/((double) RAND_MAX) + dx;
                                                double sub_y = rand()/((double) RAND_MAX) + dy;

                                                Point3D origin(0, 0, 0);
                                                Point3D imagePlane;
                                                imagePlane[0] = (-double(image.width)/2 + sub_y + j)/factor;
                                                imagePlane[1] = (-double(image.height)/2 + sub_x + i)/factor;
                                                imagePlane[2] = -1;

                                                Ray3D ray;
                                                // TODO: Convert ray to world space
                                                Vector3D direction = imagePlane - origin;
                                                direction = viewToWorld * direction;
                                                origin = viewToWorld * origin;
                                                ray = Ray3D(origin, direction);

                                                // render with depth of field
                                                if (DOF) {
                                                        // ------ variable for Depth of Field ----------
                                                        double focal_len = 5.0;
                                                        double apeture_size = 1.0;
                                                        // ray_num samples for calculating DOF
                                                        int dof_num = 10;
                                                        // ---------------------------------------------

                                                        Color dofCol(0.0, 0.0, 0.0);

                                                        for (int count = 0; count < dof_num; count++) {
                                                                Ray3D secondary_ray = depthOfField(ray, origin, focal_len, apeture_size);
                                                                dofCol = dofCol + shadeRay(secondary_ray, scene, light_list, depth);
                                                        }
                                                        col = col + double(1.0/dof_num)*dofCol;

                                                }else{
                                                        col = col + shadeRay(ray, scene, light_list, depth);
                                                }

                                        }
                                }
                        }

                        // uniform distribution
                        col = double(1.0/(4*AA_num))*col;
                        image.setColorAtPixel(i, j, col);

                }
        }
}



// --------------- helper function for Depth of Field -------------------//
// calculate the secondary ray from random point on lens
Ray3D Raytracer::depthOfField(Ray3D& primary_ray, Point3D& origin, double focal_len, double apeture_size){

        // focal point
        Point3D focal_point = origin + focal_len*primary_ray.dir;

        // random point on lens
        double r1 = -apeture_size + 2*apeture_size*(rand()/((double) RAND_MAX));
        double r2 = -apeture_size + 2*apeture_size*(rand()/((double) RAND_MAX));

        Point3D point_on_lens = origin + Vector3D(r1, r2, 0.0);

        // secondary ray direction
        Vector3D secondary_dir = focal_point - point_on_lens;
        secondary_dir.normalize();

        return Ray3D(point_on_lens, secondary_dir);


}


// --------------- helper function for glossy reflection -------------------//
// helper function for glossy reflection
// calculate the new reflected ray direction
Vector3D Raytracer::glossyRefDirection(Vector3D N, Vector3D R, double roughness){
    
    //orthonormal basis at intersection point
    Vector3D u = R.cross(N);
    
    u.normalize();
    
    Vector3D v = R.cross(u);
    v.normalize();
    
    //choose uniformly sampled random direction
    double theta = 2*M_PI*((rand()/((double) RAND_MAX))*roughness);
    double phi = 2*M_PI*((rand()/((double) RAND_MAX))*roughness);
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);
    
    Vector3D newR = x*u + y*v + z*R;
    newR.normalize();
    
    return newR;
}


// -------------------- Texture Helper Func -------------------- //
Color Raytracer::textureColor(Ray3D& ray) {
        //std::cout << "getting the textureCol ..." << '\n';
        Material* mat = ray.intersection.mat;
        if (!mat->hasTexture) {
                std::cout << "Error, has no texture" << '\n';
        } else {
                int x = ray.intersection.u * mat->t_width;
                int y = ray.intersection.v * mat->t_height;
                int i = y * mat->t_width + x;
                double r = double(mat->rarray[i])/255.0;
                double g = double(mat->garray[i])/255.0;
                double b = double(mat->barray[i])/255.0;

                return Color(r,g,b);

        }

        return Color(0, 0, 0);
}
