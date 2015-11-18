#include "intersect.h"
#include "scene.h"
#include "montecarlo.h"
#include "fstream"


#include <iostream>
#include <thread>
using std::thread;
#define MAX_DEPTH     10


// modify the following line to disable/enable parallel execution of the raytracer
bool parallel_raytrace = true;

image3f ray_trace(Scene* scene, bool multithread);
void ray_trace(Scene* scene, image3f* image, int offset_row, int skip_row, bool verbose);


// lookup texture value
vec3f lookup_scaled_texture(vec3f value, image3f* texture, vec2f uv, bool tile = false) {
    if(not texture) return value;

    // for now, simply clamp texture coords
    if (!tile) {
        auto u = clamp(uv.x, 0.0f, 1.0f);
        auto v = clamp(uv.y, 0.0f, 1.0f);
        return value * texture->at(u*(texture->width()-1), v*(texture->height()-1));
    }
    else {
        auto u = uv.x >= 0 ? uv.x - (int)uv.x : uv.x - (int)uv.x + 1;
        auto v = uv.y >= 0 ? uv.y - (int)uv.y : uv.y - (int)uv.y + 1;
        return value * texture->at(u*(texture->width()-1), v*(texture->height()-1));
    }
}

// compute the brdf
vec3f eval_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f l, vec3f norm, bool microfacet) {
    if (not microfacet) {
        auto h = normalize(v+l);
        return kd/pif + ks*(n+8)/(8*pif) * pow(max(0.0f,dot(norm,h)),n);
    } else {
        put_your_code_here("Implement microfacet brdf");
        return one3f; // <- placeholder
    }
}

// evaluate the environment map
vec3f eval_env(vec3f ke, image3f* ke_txt, vec3f dir) {
    if (ke_txt != nullptr) {
        float u = atan2(dir.x, dir.z) / (2 * pif);
        float v = 1 - acos(dir.y) / pif;
        return lookup_scaled_texture(ke, ke_txt, vec2f(u, v), true);
    }
    else {
        return ke; // <- placeholder
    }
}

vec3f infract(vec3f d, vec3f n, float last_nt, float nt){
    float m = 1 - last_nt * last_nt * ( 1 - dot(d,n) * dot(d,n)) / (nt * nt);
    if(m < 0)   return zero3f;
    return normalize( last_nt * (d - dot(d,n) * n) / nt - n * sqrt(m));
}

// compute the color corresponing to a ray by raytracing
vec3f raytrace_ray(Scene* scene, const ray3f& ray, Rng *rng, int depth,int count, int countC) {
    // get scene intersection
    auto intersection = intersect(scene,ray);

    // if not hit, return background (looking up the texture by converting the ray direction to latlong around y)
    if(not intersection.hit) {
        return eval_env(scene->background, scene->background_txt, ray.d);
    }

    // setup variables for shorter code
    auto pos = intersection.pos;
    auto norm = intersection.norm;
    auto v = -ray.d;

    // compute material values by looking up textures
    auto ke = lookup_scaled_texture(intersection.mat->ke, intersection.mat->ke_txt, intersection.texcoord);
    auto kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord);
    auto ks = lookup_scaled_texture(intersection.mat->ks, intersection.mat->ks_txt, intersection.texcoord);
    auto n = intersection.mat->n;
    auto mf = intersection.mat->microfacet;

    // accumulate color starting with ambient
    auto c = scene->ambient * kd;

    // add emission if on the first bounce
    if(depth == 0 and dot(v,norm) > 0) c += ke;

    // foreach surface
    // skip if no emission from surface
    // todo: pick a point on the surface, grabbing normal, area, and texcoord
    // check if quad
    // generate a 2d random number
    // compute light position, normal, area
    // set tex coords as random value got before
    // else if sphere
    // generate a 2d random number
    // compute light position, normal, area
    // set tex coords as random value got before
    // get light emission from material and texture
    // compute light direction
    // compute light response (ke * area * cos_of_light / dist^2)
    // compute the material response (brdf*cos)
    // multiply brdf and light
    // check for shadows and accumulate if needed
    // if shadows are enabled
    // perform a shadow check and accumulate
    // else just accumulate

    // area light
    for(Surface* s : scene->surfaces) {
        if (s->mat->ke == zero3f) continue;

        vec2f sur_texcoord;
        vec3f sur_norm, sur_pos;
        float sur_area;
        if (s->isquad) {
            sur_texcoord = rng->next_vec2f();
            sur_pos = s->frame.o +
            (sur_texcoord.x - 0.5) * 2 * s->radius * s->frame.x +
            (sur_texcoord.y - 0.5) * 2 * s->radius * s->frame.y;
            sur_norm = s->frame.z;
            sur_area = 2 * s->radius * 2*s->radius;
        }
        else {
            sur_texcoord = rng->next_vec2f();
            sur_norm = normalize(sample_direction_spherical_uniform(sur_texcoord));
            sur_pos = s->frame.o + s->radius * sur_norm;
            sur_area = 4 * pif * sqr(s->radius);
        }
        vec3f ke, d, l, brdfcos, shade;
        float l_cos;
        ke = lookup_scaled_texture(s->mat->ke, s->mat->ke_txt, sur_texcoord);
        d = normalize(sur_pos - pos);
        l_cos = max(-dot(sur_norm, d), 0.0f);
        l = ke * sur_area * l_cos / lengthSqr(sur_pos - pos);
        brdfcos = max(dot(norm, d), 0.0f) * eval_brdf(kd, ks, n, v, d, norm, mf);
        shade = l * brdfcos;
        if (shade == zero3f) continue;
        if (scene->path_shadows) {
            if (!intersect_shadow(scene, ray3f::make_segment(pos, sur_pos))) {
                c += shade;
            }
        }
        else {
            c += shade;
        }
    }

    // todo: sample the brdf for environment illumination if the environment is there
    // if scene->background is not zero3f
    if (scene->background != zero3f) {
        pair<vec3f, float> res = sample_brdf(kd, ks, n, v, norm, rng->next_vec2f(), rng->next_float());

        // compute the material response (brdf*cos)
        vec3f brdfcos = max(dot(res.first, norm), 0.0f) * eval_brdf(kd, ks, n, v, normalize(res.first), norm, mf);
        // todo: accumulate response scaled by brdf*cos/pdf
        vec3f shade = brdfcos * eval_env(scene->background, scene->background_txt, normalize(res.first)) / res.second;
        // if material response not zero3f
        if (brdfcos != zero3f) {
            // if shadows are enabled
            if (scene->path_shadows) {
                // perform a shadow check and accumulate
                if (not intersect_shadow(scene, ray3f(pos, res.first))) {
                    c += shade;
                }
            }
            // else just accumulate
            else {
                c += shade;
            }
        }
    }

    //indirect illumination
    if (depth < scene->path_max_depth && kd != zero3f && ks != zero3f) {
        pair<vec3f, float> res = sample_brdf(kd, ks, n, v, norm, rng->next_vec2f(), rng->next_float());
        vec3f brdfcos = max(dot(res.first, norm), 0.0f) * eval_brdf(kd, ks, n, v, normalize(res.first), norm, mf);
        ray3f newRay = ray3f(pos, normalize(res.first));
        c += brdfcos / res.second * raytrace_ray(scene, newRay, rng, depth+1, count, countC);
    }
    
    if(count >= 4 || countC >=4)  return zero3f;
    if(intersection.hit){
        vec3f color(0,0,0);
        vec3f m = 2*(abs(dot(ray.d, intersection.norm)) * intersection.norm + ray.d);
        vec3f newRay = m-ray.d;

        //refact
        if(intersection.mat->transparent){
            vec3f t;
            float cf;
            if(dot(intersection.norm,ray.d) < 0){
                t = infract(ray.d, intersection.norm, 1, intersection.mat->nt);
                cf = -dot(intersection.norm,ray.d);
            }else{
                t = infract(ray.d, -intersection.norm, intersection.mat->nt, 1);
                cf = dot(t,intersection.norm);
            }

            float R0 = sqr(intersection.mat->nt - 1) / sqr(intersection.mat->nt + 1);
            float R = R0 + (1-R0) * pow(1-cf, 5);
            if(t != zero3f){
                c += 0.4*raytrace_ray(scene, ray3f(intersection.pos + 0.001 * ray.d, t), rng, depth, count, countC+1);
            }
        }

        //reflect
        if(intersection.mat->kr != zero3f){
            color += (intersection.mat->kr * raytrace_ray(scene, ray3f(intersection.pos, newRay), rng, depth, count+1, countC)) ;
        }

        for(int i = 0 ; i < scene->lights.size() ; ++i){
            vec3f l = scene->lights[i]->frame.o - intersection.pos;
            double dis = (abs(dot(l,l)));
            l = normalize(l);

            //diffuse
            double mul = dot(intersection.norm,l);
            color += kd * (scene->lights[i]->intensity/dis) * (mul > 0 ? mul : 0) ;

            //specular
            vec3f v = normalize(scene->camera->frame.o - intersection.pos);
            vec3f h = normalize(l + v);
            mul = dot(intersection.norm, h);
            color += ks * (scene->lights[i]->intensity/dis) * pow((mul > 0 ? mul : 0),n) ;

            ray3f raySha(intersection.pos, normalize(scene->lights[i]->frame.o - (intersection.pos)));
            int shadowN = 5;
            double sampleLength = 0.5;
            srand((unsigned)time(NULL));
            vec3f colorShadow = zero3f;
            double ope = 1.0;
            for(int shadowCount = 0 ; shadowCount < shadowN*shadowN ; ++shadowCount){
//                vec3f posShadow(intersection.pos.x + ope* rand()/RAND_MAX * sampleLength, intersection.pos.y + ope * rand()/RAND_MAX * sampleLength, intersection.pos.z + ope * rand()/RAND_MAX * sampleLength);
                vec3f posShadow = intersection.pos+0.1*intersection.norm;
                vec3f lightPos(scene->lights[i]->frame.o.x + (rand() > RAND_MAX/2.0 ? 1.0 : -1.0)* rand()/RAND_MAX * sampleLength,
                                   scene->lights[i]->frame.o.y + (rand() > RAND_MAX/2.0 ? 1.0 : -1.0)* rand()/RAND_MAX * sampleLength,
                                   scene->lights[i]->frame.o.z + (rand() > RAND_MAX/2.0 ? 1.0 : -1.0)* rand()/RAND_MAX * sampleLength);
                ray3f rayShadow(posShadow, normalize(lightPos - (posShadow)));
                if(! intersect(scene,rayShadow).hit){
                    colorShadow += color;
                }
                ope *= -1.0;
            }
            colorShadow /= shadowN * shadowN;
            c += colorShadow;
        }
    }

    return c;
}

// raytrace an image
void ray_trace(Scene* scene, image3f* image, RngImage *rngs,int offset_row, int skip_row, bool verbose) {

    // if no anti-aliasing
    if (scene->image_samples == 1) {
        // foreach image row (go over image height)
        for( int y = offset_row; y < scene->image_height; y+= skip_row)
        {
            if(verbose) message("\r  rendering1 %03d/%03d        ", y, scene->image_height);
            // foreach pixel in the row (go over image width)
            for( int x = 0; x < scene->image_width; x++ )
            {   if(verbose) message("\r  renderingx %03d/%03d        ", x, scene->image_height);
                // compute ray-camera parameters (u,v) for the pixel
                auto rng = &rngs->at(x, y);
                float u = (x + 0.5f) / image->width();
                float v = (y + 0.5f) / image->height();
                const vec3f Q_uv( ( u-0.5f ) * scene->camera->width, ( v-0.5f ) * scene->camera->height, - scene->camera->dist);
                // compute camera ray
                const ray3f view_ray = transform_ray( scene->camera->frame, ray3f( zero3f, normalize(Q_uv) ) );
                // set pixel to the color raytraced with the ray
                vec3f color = raytrace_ray(scene, view_ray, rng, 0, 0, 0);
                image->at(x, y) = color;
            }
        }
    }
    // else
    else{
        vector<vector<int>> countColor(image->width(), vector<int>(image->height(),0));
        // foreach image row (go over image height)
        for( int y = offset_row; y < scene->image_height; y+= skip_row )
        {
            if(verbose) message("\r  rendering2 %03d/%03d        ", y, scene->image_height);
            // foreach pixel in the row (go over image width)
            for (int x = 0; x < scene->image_width; x++) {
                if(verbose) message("\r  renderingx %03d/%03d        ", x, scene->image_height);
                // init accumulated color
                vec3f color = zero3f;
                // foreach sample in y
                auto rng = &rngs->at(x, y);
                for (int j = 0; j < scene->image_samples; j++) {
                    // foreach sample in x
                    for( int i = 0; i < scene->image_samples; i++ )
                    {
                        // compute ray-camera parameters (u,v) for the pixel and the sample
                        auto u = (x + (i + 0.5f)/scene->image_samples) / image->width();
                        auto v = (y + (j + 0.5f)/scene->image_samples) / image->height();
                        // compute camera ray
                        vec3f Q_uv( ( u-0.5f ) * scene->camera->width, ( v-0.5f ) * scene->camera->height, - scene->camera->dist);
                        ray3f view_ray = transform_ray( scene->camera->frame, ray3f( zero3f, normalize(Q_uv) ) );
                        // set pixel to the color raytraced with the ray
                        color += raytrace_ray(scene, view_ray, rng, 0,0 ,0);
                    }
                }
                // scale by the number of samples
                color /= std::pow(scene->image_samples,2);
                image->at(x, y) = color;
            }
        }

    }
}

void test()
{
    //const ray3f ray( zero3f, vec3f(1,0,0) );
    //std::cout << ray_sphere_intersection(ray, vec3f(5,0,0), 2) << '\n';
}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {

    test();

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

    message("accelerating...\n");
    accelerate(scene);

    message("rendering %s...\n", scene_filename.c_str());
    auto image = ray_trace(scene, parallel_raytrace);

    message("\nwriting to png...\n");
    write_png(image_filename, image, true);

    delete scene;
    message("done\n");
}


// pathtrace an image with multithreading if necessary
image3f ray_trace(Scene* scene, bool multithread) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);

    // create a random number generator for each pixel
    auto rngs = RngImage(scene->image_width, scene->image_height);

    // if multitreaded
    if(multithread) {
        // get pointers
        auto image_ptr = &image;
        auto rngs_ptr = &rngs;
        // allocate threads and pathtrace in blocks
        auto threads = vector<thread>();
        auto nthreads = thread::hardware_concurrency();
        for(auto tid : range(nthreads)) threads.push_back(thread([=](){
            return ray_trace(scene,image_ptr,rngs_ptr,tid,nthreads,tid==0);}));
        for(auto& thread : threads) thread.join();
    } else {
        // pathtrace all rows
        ray_trace(scene, &image, &rngs, 0, 1, true);
    }

    // done
    return image;
}
