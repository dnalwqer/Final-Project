#include "scene.h"
#include "float.h"
#include "fstream"
std::ofstream offile;
// intersection record
struct intersection3f {
    bool        hit;        // whether it hits something
    double       ray_t;      // ray parameter for the hit
    vec3f       pos;        // hit position
    vec3f       norm;       // hit normal
    Material*   mat;        // hit material
    
    // constructor (defaults to no intersection)
    intersection3f() : hit(false) { }
};

#define ray3f_epsilon 0.0005f
#define ray3f_rayinf 1000000.0f

// 3D Ray
struct ray3f {
    vec3f e;        // origin
    vec3f d;        // direction
    double tmin;     // min t value
    double tmax;     // max t value
    
    // Default constructor
    ray3f() : e(zero3f), d(z3f), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }
    
    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d) :
    e(e), d(d), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }
    
    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d, double tmin, double tmax) :
    e(e), d(d), tmin(tmin), tmax(tmax) { }
    
    // Eval ray at a specific t
    vec3f eval(double t) const { return e + d * t; }
    
    // Create a ray from a segment
    static ray3f make_segment(const vec3f& a, const vec3f& b) { return ray3f(a,normalize(b-a),ray3f_epsilon,dist(a,b)-2*ray3f_epsilon); }
};

// transform a ray by a frame
inline ray3f transform_ray(const frame3f& f, const ray3f& v) {
    return ray3f(transform_point(f,v.e), transform_vector(f,v.d), v.tmin, v.tmax);
}
// transform a ray by a frame inverse
inline ray3f transform_ray_inverse(const frame3f& f, const ray3f& v) {
    return ray3f(transform_point_inverse(f,v.e),transform_vector_inverse(f,v.d),v.tmin,v.tmax);
}

double pointMul(vec3f a, vec3f b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

// intersects the scene and return the first intrerseciton
intersection3f intersect(Scene* scene, ray3f ray) {
    // create a default intersection record to be returned
    auto intersection = intersection3f();
    long double length = FLT_MAX;

    for(int i = 0 ; i < scene->surfaces.size() ; ++i){
        if(!scene->surfaces[i]->isquad && !scene->surfaces[i]->iscli){
            vec3f circle = scene->surfaces[i]->frame.o;
            long double radius = scene->surfaces[i]->radius;
            if(pointMul(circle-ray.e, ray.d)<0) continue;
            long double disFE = pointMul(ray.d, circle-ray.e);
            long double dis = sqrt(pointMul(circle-ray.e, circle-ray.e) - disFE*disFE);
            if(dis < radius){
                long double disC = disFE-sqrt(radius*radius - dis*dis);
                long double t = sqrt((disC*disC)/(ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z));
                if(t < length && t<ray.tmax && t>ray.tmin){
                    length = t;
                    intersection.hit = true;
                    if((ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z) == 0) continue;
                    intersection.ray_t = t;
                    intersection.pos = ray.e + ray.d*t;//vec3f(ray.e.x + ray.d.x*intersection.ray_t, ray.e.y + ray.d.y*intersection.ray_t, ray.e.z + ray.d.z*intersection.ray_t);
                    intersection.norm = normalize(intersection.pos - scene->surfaces[i]->frame.o);
                    intersection.mat = scene->surfaces[i]->mat;
                }
            }
        }else if(scene->surfaces[i]->iscli){
            long double radius = scene->surfaces[i]->radius;

            //wall
            ray.e.x -= scene->surfaces[i]->frame.o.x;
            ray.e.z -= scene->surfaces[i]->frame.o.z;
            double wallW = 4 * (ray.d.x*ray.e.x + ray.d.z*ray.e.z) * (ray.d.x*ray.e.x + ray.d.z*ray.e.z) -
                           4 * (ray.d.x*ray.d.x + ray.d.z*ray.d.z) * (ray.e.x*ray.e.x + ray.e.z*ray.e.z - radius*radius);
            double t = (-2 * (ray.d.x*ray.e.x + ray.d.z*ray.e.z) - sqrt(abs(wallW))) / (2*(ray.d.x*ray.d.x + ray.d.z*ray.d.z));
            ray.e.x += scene->surfaces[i]->frame.o.x;
            ray.e.z += scene->surfaces[i]->frame.o.z;
            if(wallW > 0 && t > 0){
                double ty = ray.d.y*t + ray.e.y;
                if(abs(ty - scene->surfaces[i]->frame.o.y) < scene->surfaces[i]->height/2 &&
                        t < length && t<ray.tmax && t>ray.tmin){
                    length = t;
                    intersection.hit = true;
                    intersection.ray_t = t;
                    intersection.pos = ray.e + ray.d*t;
                    intersection.norm = normalize(intersection.pos - vec3f(scene->surfaces[i]->frame.o.x,intersection.pos.y,scene->surfaces[i]->frame.o.z));
                    intersection.mat = scene->surfaces[i]->mat;
                }
            }

            //caps
            int twoCap[2] = {1,-1};
            for(int c = 0 ; c < 2 ; ++c){
                t = (scene->surfaces[i]->frame.o.y + twoCap[c]*radius/2 - ray.e.y) / ray.d.y;
                vec3f pos = ray.d * t + ray.e;
                vec3f cirCap = vec3f(scene->surfaces[i]->frame.o.x, scene->surfaces[i]->frame.o.y+twoCap[c]*radius/2, scene->surfaces[i]->frame.o.z);
                if(sqrt(pointMul(cirCap-pos,cirCap-pos)) < radius && t < length && t<ray.tmax && t>ray.tmin){
                    length = t;
                    intersection.hit = true;
                    intersection.ray_t = t;
                    intersection.pos = pos;
                    intersection.norm = normalize(intersection.pos - vec3f(scene->surfaces[i]->frame.o.x,intersection.pos.y,scene->surfaces[i]->frame.o.z));
                    intersection.mat = scene->surfaces[i]->mat;
                }
            }
        }
        else{
            long double vpt = pointMul(ray.d, scene->surfaces[i]->frame.z);//ray.d.x * scene->surfaces[i]->frame.z.x + ray.d.y * scene->surfaces[i]->frame.z.y + ray.d.z * scene->surfaces[i]->frame.z.z;
            if(vpt == 0)    continue;
            long double t = ((scene->surfaces[i]->frame.o.x - ray.e.x) * scene->surfaces[i]->frame.z.x + (scene->surfaces[i]->frame.o.y- ray.e.y) * scene->surfaces[i]->frame.z.y + (scene->surfaces[i]->frame.o.z - ray.e.z) * scene->surfaces[i]->frame.z.z) / vpt;
            vec3f pos = ray.e + ray.d*t;

            normalize(scene->surfaces[i]->frame.x);
            normalize(scene->surfaces[i]->frame.y);
            if(!(abs(pointMul(pos - scene->surfaces[i]->frame.o, scene->surfaces[i]->frame.x)) <= scene->surfaces[i]->radius &&
               abs(pointMul(pos - scene->surfaces[i]->frame.o, scene->surfaces[i]->frame.y)) <= scene->surfaces[i]->radius &&
               pos.z < scene->camera->frame.o.z))
                 continue;
            if(t < length && t < ray.tmax && t > ray.tmin)
            {
                intersection.hit = true;
                length = t;
                intersection.pos = pos;
                intersection.ray_t = t;
                intersection.norm = normalize(scene->surfaces[i]->frame.z);
                intersection.mat = scene->surfaces[i]->mat;
            }
        }
    }
    // foreach surface
        // if it is a quad
            // compute ray intersection (and ray parameter), continue if not hit
            // check if computed param is within ray.tmin and ray.tmax
            // check if this is the closest intersection, continue if not
            // if hit, set intersection record values
        // if it is a sphere
            // compute ray intersection (and ray parameter), continue if not hit
            // just grab only the first hit
            // check if computed param is within ray.tmin and ray.tmax
            // check if this is the closest intersection, continue if not
    // record closest intersection
    return intersection;
}

// compute the color corresponing to a ray by raytracing
vec3f raytrace_ray(Scene* scene, ray3f ray,int count) {
    if(count >= 4)  return zero3f;
    intersection3f intersec = intersect(scene, ray);
    if(intersec.hit){
        vec3f color(0,0,0);
        //reflect
        if(intersec.mat->kr != zero3f){
            vec3f m = 2*(abs(pointMul(ray.d, intersec.norm)) * intersec.norm + ray.d);
            vec3f newRay = m-ray.d;
            color += (intersec.mat->kr * raytrace_ray(scene, ray3f(intersec.pos, newRay), count+1)) ;
        }

        for(int i = 0 ; i < scene->lights.size() ; ++i){
            vec3f l = scene->lights[i]->frame.o - intersec.pos;
            double dis = (abs(dot(l,l))); offile<<dis<<"\r\n";
            l = normalize(l);
            //ambient
            color += scene->ambient *  (scene->lights[i]->intensity);

            //shadow
            ray3f raySha(intersec.pos, normalize(scene->lights[i]->frame.o - (intersec.pos)));
            intersection3f intersecSha = intersect(scene,raySha);
            if(intersecSha.hit) continue;

            //diffuse
            double mul = dot(intersec.norm,l);
            color += intersec.mat->kd * (scene->lights[i]->intensity/dis) * (mul > 0 ? mul : 0) ;

            //specular
            vec3f v = normalize(scene->camera->frame.o - intersec.pos);
            vec3f h = normalize(l + v);
            mul = pointMul(intersec.norm, h);
            color += intersec.mat->ks * (scene->lights[i]->intensity/dis) * pow((mul > 0 ? mul : 0),intersec.mat->n) ;
        }
        return color;
    }
    return scene->background;
    // get scene intersection
    // if not hit, return background
    // accumulate color starting with ambient
    // foreach light
        // compute light response
        // compute light direction
        // compute the material response (brdf*cos)
        // check for shadows and accumulate if needed
    // if the material has reflections
        // create the reflection ray
        // accumulate the reflected light (recursive call) scaled by the material reflection
    // return the accumulated color (for now zero)
    return zero3f;
}

// raytrace an image
image3f raytrace(Scene* scene) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);
    
    offile.open("data.txt");
    for(int i = 0 ; i < scene->image_width ; ++i){
      for(int j = 0 ; j < scene->image_height ; ++j){
           vec3f color;
           for(int m = 0 ; m < scene->image_samples ; ++m){
               for(int n = 0 ; n < scene->image_samples ; ++n){
                   float u = -0.5 + (i+(m+0.5)/scene->image_samples) / scene->image_width;
                   float v = -0.5 + (j+(n+0.5)/scene->image_samples) / scene->image_height;
                   vec3f direction = normalize(u * scene->camera->frame.x + v * scene->camera->frame.y - scene->camera->frame.z);
                   ray3f ray(scene->camera->frame.o, direction);
                   color += raytrace_ray(scene,ray,0);
               }
           }
           color /= scene->image_samples*scene->image_samples;

             image.setPix(i,j,color);
           }
       }
    // if no anti-aliasing
        // foreach image row (go over image height)
            // foreach pixel in the row (go over image width)
                // compute ray-camera parameters (u,v) for the pixel
                // compute camera ray
                // set pixel to the color raytraced with the ray
    // else
        // foreach image row (go over image height)
            // foreach pixel in the row (go over image width)
                // init accumulated color
                // foreach sample in y
                    // foreach sample in x
                        // compute ray-camera parameters (u,v) for the pixel and the sample
                        // compute camera ray
                        // set pixel to the color raytraced with the ray
                // scale by the number of samples
    
    // done
    return image;
}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "01_raytrace", "raytrace a scene",
            {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue()}  },
            {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json")},
               {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("")}  }
        });
    
    // generate/load scene either by creating a test scene or loading from json file
    string scene_filename = args.object_element("scene_filename").as_string();
    Scene *scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");
    
    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    
    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    
    message("rendering %s...\n", scene_filename.c_str());
    auto image = raytrace(scene);
    
    message("writing to png...\n");
    write_png(image_filename, image, true);
    
    delete scene;
    message("done\n");
}
