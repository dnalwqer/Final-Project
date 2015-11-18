#include "common.h"
#include "scene.h"
#include "image.h"
#include "gls.h"

string scene_filename;  // scene filename
string image_filename;  // image filename
Scene* scene;           // scene arrays

void uiloop();          // UI loop


// map used to uniquify edges
struct EdgeMap {
    map<pair<int,int>,int>  _edge_map;  // internal map
    vector<vec2i>           _edge_list; // internal list to generate unique ids
    
    // create an edge map for a collection of triangles and quads
    EdgeMap(vector<vec3i> triangle, vector<vec4i> quad) {
        for(auto f : triangle) { _add_edge(f.x,f.y); _add_edge(f.y,f.z); _add_edge(f.z,f.x); }
        for(auto f : quad) { _add_edge(f.x,f.y); _add_edge(f.y,f.z); _add_edge(f.z,f.w); _add_edge(f.w,f.x); }
    }
    
    // internal function to add an edge
    void _add_edge(int i, int j) {
        if(_edge_map.find(make_pair(i,j)) == _edge_map.end()) {
            _edge_map[make_pair(i,j)] = _edge_list.size();
            _edge_map[make_pair(j,i)] = _edge_list.size();
            _edge_list.push_back(vec2i(i,j));
        }
    }
    
    // edge list
    const vector<vec2i>& edges() const { return _edge_list; }
    
    // get an edge from two vertices
    int edge_index(vec2i e) const {
        error_if_not(not (_edge_map.find(make_pair(e.x,e.y)) == _edge_map.end()), "non existing edge");
        return _edge_map.find(make_pair(e.x, e.y))->second;
    }
};


// make normals for each face - duplicates all vertex data
void facet_normals(Mesh* mesh) {
    // allocates new arrays
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();
    auto triangle = vector<vec3i>();
    auto quad = vector<vec4i>();
    
    // foreach triangle
    for(auto f : mesh->triangle) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));
        // add triangle
        triangle.push_back({nv,nv+1,nv+2});
        // add vertex data
        for(auto i : range(3)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    
    // foreach quad
    for(auto f : mesh->quad) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));
        // add quad
        quad.push_back({nv,nv+1,nv+2,nv+3});
        // add vertex data
        for(auto i : range(4)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    
    // set back mesh data
    mesh->pos = pos;
    mesh->norm = norm;
    mesh->texcoord = texcoord;
    mesh->triangle = triangle;
    mesh->quad = quad;
}

// smooth out normal - does not duplicate data
void smooth_normals(Mesh* mesh) {
    // set normals array to the same length as pos and init all elements to zero
    mesh->norm.assign(mesh->pos.size(), zero3f);

    // foreach triangle
    for (auto f : mesh->triangle) {
        // compute face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));

        // accumulate face normal to the vertex normals of each face index
        for (auto i : range(3)) {
            mesh->norm[f[i]] = fn;
        }
    }
    
    // foreach quad
    for (auto f : mesh->quad) {
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));

        // accumulate face normal to the vertex normals of each face index
        for (auto i : range(4)) {
            mesh->norm[f[i]] = fn;
        }
    }
    // normalize all vertex normals

    for (auto &t : mesh->norm) {
        t = normalize(t);
    }
}

// smooth out tangents
void smooth_tangents(Mesh* polyline) {
    // set tangent array
    polyline->norm = vector<vec3f>(polyline->pos.size(),zero3f);
    // foreach line
    for(auto l : polyline->line) {
        // compute line tangent
        auto lt = normalize(polyline->pos[l.y]-polyline->pos[l.x]);
        // accumulate segment tangent to vertex tangent on each vertex
        for (auto i : range(2)) polyline->norm[l[i]] += lt;
    }
    // normalize all vertex tangents
    for (auto& t : polyline->norm) t = normalize(t);
}

//surface patches by dnalwqer
void subdivide_bezier_surface(Mesh *mesh) {
    auto m = new Mesh(*mesh);

    int steps = 1 << mesh->subdivision_bezier_level;

    int row = (int)sqrt(m->point.size());
    auto pos = vector<vec3f>();
    auto index = vector<int>();
    int colSize = 0;
    int rowSize = 0;
    //blender each row to get new coordinates.

    for (int a = 0; a < row; a++) {
        int p0_index = m->point[a*4];
        int p1_index = m->point[a*4 + 1];
        int p2_index = m->point[a*4 + 2];
        int p3_index = m->point[a*4 + 3];

        vec3f p0 = m->pos[p0_index];
        vec3f p1 = m->pos[p1_index];
        vec3f p2 = m->pos[p2_index];
        vec3f p3 = m->pos[p3_index];
        // note the starting index of new points
        // foreach step
        for (int m = 0; m <= steps; m++) {
            // compute t for current segment
            float t = 1.0 / steps * m;
            // compute blending weights
            float b0 = bernstein(t, 0, 3);
            float b1 = bernstein(t, 1, 3);
            float b2 = bernstein(t, 2, 3);
            float b3 = bernstein(t, 3, 3);
            // compute new point position
            vec3f new_point = b0 * p0 + b1 * p1 + b2 * p2 + b3 * p3;
            // add new point to pos vector
            index.push_back(pos.size());
            pos.push_back(new_point);
        }
        if (a == 0) {
            colSize = pos.size();
        }
    }
    auto pos2 = vector<vec3f>();
    auto index2 = vector<int>();
    for (int a = 0; a < colSize; a++) {
        int p0_index = index[a];
        int p1_index = index[a + colSize];
        int p2_index = index[a + 2 * colSize];
        int p3_index = index[a + 3 * colSize];

        vec3f p0 = pos[p0_index];
        vec3f p1 = pos[p1_index];
        vec3f p2 = pos[p2_index];
        vec3f p3 = pos[p3_index];
        // note the starting index of new points
        // foreach step
        for (int m = 0; m <= steps; m++) {
            // compute t for current segment
            float t = 1.0 / steps * m;
            // compute blending weights
            float b0 = bernstein(t, 0, 3);
            float b1 = bernstein(t, 1, 3);
            float b2 = bernstein(t, 2, 3);
            float b3 = bernstein(t, 3, 3);
            // compute new point position
            vec3f new_point = b0 * p0 + b1 * p1 + b2 * p2 + b3 * p3;
            // add new point to pos vector
            index2.push_back(pos2.size());
            pos2.push_back(new_point);
        }
        if (a == 0) {
            rowSize = pos2.size();
        }
    }

    m->pos = pos2;
    auto quad = vector<vec4i>();
    auto norm = vector<vec3f>();
    norm.assign(m->pos.size(), zero3f);
    for (int i = 0; i < colSize - 1; i++) {
        for (int j = 0; j < rowSize - 1; j++) {
            int p0 = i * rowSize + j;
            int p1 = p0 + rowSize;
            int p2 = p1 + 1;
            int p3 = p0 + 1;
            quad.push_back(vec4i(p0, p1, p2, p3));
        }
    }
    m->quad = quad;

    for (auto f : m->quad) {
        // compute face normal
        auto fn = normalize(normalize(cross(m->pos[f.y]-m->pos[f.x], m->pos[f.z]-m->pos[f.x])) +
                            normalize(cross(m->pos[f.z]-m->pos[f.x], m->pos[f.w]-m->pos[f.x])));

        // accumulate face normal to the vertex normals of each face index
        for (auto i : range(4)) {
            norm[f[i]] = fn;
        }
    }
    // normalize all vertex normals

    for (auto &t : norm) {
        t = normalize(t);
    }

    m->norm = norm;
    // clear subdivision
    m->subdivision_bezier_level = 0;
    *mesh = *m;
    delete m;
}

void subdivide(Scene* scene) {
    for(auto mesh : scene->meshes) {
        if(mesh->subdivision_bezier_level) subdivide_bezier_surface(mesh);
    }
}

// main function
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "02_model", "view scene",
            {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue() }  },
            {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json")},
               {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("")}  }
        });
    
    // generate/load scene either by creating a test scene or loading from json file
    scene_filename = args.object_element("scene_filename").as_string();
    scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");
    
    image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    
    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    
    subdivide(scene);
    
    uiloop();
}




/////////////////////////////////////////////////////////////////////
// UI and Rendering Code: OpenGL, GLFW, GLSL


int gl_program_id         = 0;  // OpenGL program handle
int gl_vertex_shader_id   = 0;  // OpenGL vertex shader handle
int gl_fragment_shader_id = 0;  // OpenGL fragment shader handle
map<image3f*,int> gl_texture_id;// OpenGL texture handles

bool save      = false;         // whether to start the save loop
bool wireframe = false;         // display as wireframe

void init_shaders();            // initialize the shaders
void init_textures();           // initialize the textures
void shade();                   // render the scene with OpenGL
void _shade_mesh(Mesh* mesh);
void character_callback(GLFWwindow* window, unsigned int key);  // ...
                                // glfw callback for character input
void _bind_texture(string name_map, string name_on, image3f* txt, int pos); // ...
                                // utility to bind texture parameters for shaders
                                // uses texture name, texture_on name, texture pointer and texture unit position

// glfw callback for character input
void character_callback(GLFWwindow* window, unsigned int key) {
    if(key == 's') save = true;
    if(key == 'w') wireframe = not wireframe;
}

// uiloop
void uiloop() {
    auto ok_glfw = glfwInit();
    error_if_not(ok_glfw, "glfw init error");
    
    glfwWindowHint(GLFW_SAMPLES, scene->image_samples);
    
    auto window = glfwCreateWindow(scene->image_width,
                                   scene->image_height,
                                   "graphics13 | model", NULL, NULL);
    error_if_not(window, "glfw window error");
    
    glfwMakeContextCurrent(window);
    
    glfwSetCharCallback(window, character_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    
    auto ok_glew = glewInit();
    error_if_not(GLEW_OK == ok_glew, "glew init error");
    
    init_shaders();
    init_textures();
    
    auto mouse_last_x = -1.0;
    auto mouse_last_y = -1.0;
    
    while(not glfwWindowShouldClose(window)) {
        glfwGetFramebufferSize(window, &scene->image_width, &scene->image_height);
        scene->camera->width = (scene->camera->height * scene->image_width) / scene->image_height;
        
        shade();

        if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)) {
            double x, y;
            glfwGetCursorPos(window, &x, &y);
            if (mouse_last_x < 0 or mouse_last_y < 0) { mouse_last_x = x; mouse_last_y = y; }
            auto delta_x = x - mouse_last_x, delta_y = y - mouse_last_y;
            
            set_view_turntable(scene->camera, delta_x*0.01, -delta_y*0.01, 0, 0, 0);
            
            mouse_last_x = x;
            mouse_last_y = y;
        } else { mouse_last_x = -1; mouse_last_y = -1; }
        
        if(save) {
            auto image = image3f(scene->image_width,scene->image_height);
            glReadPixels(0, 0, scene->image_width, scene->image_height, GL_RGB, GL_FLOAT, &image.at(0,0));
            write_png(image_filename, image, true);
            save = false;
        }
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    glfwDestroyWindow(window);
    
    glfwTerminate();
}

// initialize the shaders
void init_shaders() {
    // load shader code from files
    auto vertex_shader_code    = load_text_file("model_vertex.glsl");
    auto fragment_shader_code  = load_text_file("model_fragment.glsl");
    auto vertex_shader_codes   = (char *)vertex_shader_code.c_str();
    auto fragment_shader_codes = (char *)fragment_shader_code.c_str();

    // create shaders
    gl_vertex_shader_id   = glCreateShader(GL_VERTEX_SHADER);
    gl_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);
    
    // load shaders code onto the GPU
    glShaderSource(gl_vertex_shader_id,1,(const char**)&vertex_shader_codes,nullptr);
    glShaderSource(gl_fragment_shader_id,1,(const char**)&fragment_shader_codes,nullptr);
    
    // compile shaders
    glCompileShader(gl_vertex_shader_id);
    glCompileShader(gl_fragment_shader_id);
    
    // check if shaders are valid
    error_if_glerror();
    error_if_shader_not_valid(gl_vertex_shader_id);
    error_if_shader_not_valid(gl_fragment_shader_id);
    
    // create program
    gl_program_id = glCreateProgram();
    
    // attach shaders
    glAttachShader(gl_program_id,gl_vertex_shader_id);
    glAttachShader(gl_program_id,gl_fragment_shader_id);
    
    // bind vertex attributes locations
    glBindAttribLocation(gl_program_id, 0, "vertex_pos");
    glBindAttribLocation(gl_program_id, 1, "vertex_norm");
    glBindAttribLocation(gl_program_id, 2, "vertex_texcoord");

    // link program
    glLinkProgram(gl_program_id);
    
    // check if program is valid
    error_if_glerror();
    error_if_program_not_valid(gl_program_id);
}

// initialize the textures
void init_textures() {
    // grab textures from scene
    auto textures = get_textures(scene);
    // foreach texture
    for(auto texture : textures) {
        // if already in the gl_texture_id map, skip
        if(gl_texture_id.find(texture) != gl_texture_id.end()) continue;
        // gen texture id
        unsigned int id = 0;
        glGenTextures(1, &id);
        // set id to the gl_texture_id map for later use
        gl_texture_id[texture] = id;
        // bind texture
        glBindTexture(GL_TEXTURE_2D, id);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
        // load texture data
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                     texture->width(), texture->height(),
                     0, GL_RGB, GL_FLOAT, texture->data());
    }
}


// utility to bind texture parameters for shaders
// uses texture name, texture_on name, texture pointer and texture unit position
void _bind_texture(string name_map, string name_on, image3f* txt, int pos) {
    // if txt is not null
    if(txt) {
        // set texture on boolean parameter to true
        glUniform1i(glGetUniformLocation(gl_program_id,name_on.c_str()),GL_TRUE);
        // activate a texture unit at position pos
        glActiveTexture(GL_TEXTURE0+pos);
        // bind texture object to it from gl_texture_id map
        glBindTexture(GL_TEXTURE_2D, gl_texture_id[txt]);
        // set texture parameter to the position pos
        glUniform1i(glGetUniformLocation(gl_program_id, name_map.c_str()), pos);
    } else {
        // set texture on boolean parameter to false
        glUniform1i(glGetUniformLocation(gl_program_id,name_on.c_str()),GL_FALSE);
        // activate a texture unit at position pos
        glActiveTexture(GL_TEXTURE0+pos);
        // set zero as the texture id
        glBindTexture(GL_TEXTURE_2D, 0);
    }
}

// render the scene with OpenGL
void shade() {
    // enable depth test
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    // disable culling face
    glDisable(GL_CULL_FACE);
    // let the shader control the points
    glEnable(GL_POINT_SPRITE);
    
    // set up the viewport from the scene image size
    glViewport(0, 0, scene->image_width, scene->image_height);
    
    // clear the screen (both color and depth) - set cleared color to background
    glClearColor(scene->background.x, scene->background.y, scene->background.z, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // enable program
    glUseProgram(gl_program_id);
    
    // bind camera's position, inverse of frame and projection
    // use frame_to_matrix_inverse and frustum_matrix
    glUniform3fv(glGetUniformLocation(gl_program_id,"camera_pos"),
                 1, &scene->camera->frame.o.x);
    glUniformMatrix4fv(glGetUniformLocation(gl_program_id,"camera_frame_inverse"),
                       1, true, &frame_to_matrix_inverse(scene->camera->frame)[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(gl_program_id,"camera_projection"),
                       1, true, &frustum_matrix(-scene->camera->dist*scene->camera->width/2, scene->camera->dist*scene->camera->width/2,
                                                -scene->camera->dist*scene->camera->height/2, scene->camera->dist*scene->camera->height/2,
                                                scene->camera->dist,10000)[0][0]);
    
    // bind ambient and number of lights
    glUniform3fv(glGetUniformLocation(gl_program_id,"ambient"),1,&scene->ambient.x);
    glUniform1i(glGetUniformLocation(gl_program_id,"lights_num"),scene->lights.size());
    
    // foreach light
    auto count = 0;
    for(auto light : scene->lights) {
        // bind light position and internsity (create param name with tostring)
        glUniform3fv(glGetUniformLocation(gl_program_id,tostring("light_pos[%d]",count).c_str()),
                     1, &light->frame.o.x);
        glUniform3fv(glGetUniformLocation(gl_program_id,tostring("light_intensity[%d]",count).c_str()),
                     1, &light->intensity.x);
        count++;
    }
    
    // foreach mesh
    for(auto mesh : scene->meshes) {
        _shade_mesh(mesh);
    }
    
    for(auto surf : scene->surfaces) {
        _shade_mesh(surf->_display_mesh);
    }
}

void _shade_mesh(Mesh* mesh) {
    // bind material kd, ks, n
    ERROR_IF_NOT(mesh, "mesh is null");
    glUniform3fv(glGetUniformLocation(gl_program_id,"material_kd"),1,&mesh->mat->kd.x);
    glUniform3fv(glGetUniformLocation(gl_program_id,"material_ks"),1,&mesh->mat->ks.x);
    glUniform1f(glGetUniformLocation(gl_program_id,"material_n"),mesh->mat->n);
    
    // bind texture params (txt_on, sampler)
    _bind_texture("material_kd_txt",   "material_kd_txt_on",   mesh->mat->kd_txt,   0);
    _bind_texture("material_ks_txt",   "material_ks_txt_on",   mesh->mat->ks_txt,   1);
    _bind_texture("material_norm_txt", "material_norm_txt_on", mesh->mat->norm_txt, 2);
    
    // bind mesh frame - use frame_to_matrix
    glUniformMatrix4fv(glGetUniformLocation(gl_program_id,"mesh_frame"),1,true,&frame_to_matrix(mesh->frame)[0][0]);

    // enable vertex attributes arrays and set up pointers to the mesh data
    auto vertex_pos_location = glGetAttribLocation(gl_program_id, "vertex_pos");
    auto vertex_norm_location = glGetAttribLocation(gl_program_id, "vertex_norm");
    auto vertex_texcoord_location = glGetAttribLocation(gl_program_id, "vertex_texcoord");
    glEnableVertexAttribArray(vertex_pos_location);
    glVertexAttribPointer(vertex_pos_location, 3, GL_FLOAT, GL_FALSE, 0, &mesh->pos[0].x);
    glEnableVertexAttribArray(vertex_norm_location);
    glVertexAttribPointer(vertex_norm_location, 3, GL_FLOAT, GL_FALSE, 0, &mesh->norm[0].x);
    if(not mesh->texcoord.empty()) {
        glEnableVertexAttribArray(vertex_texcoord_location);
        glVertexAttribPointer(vertex_texcoord_location, 2, GL_FLOAT, GL_FALSE, 0, &mesh->texcoord[0].x);
    }
    else glVertexAttrib2f(vertex_texcoord_location, 0, 0);
    
    // draw triangles and quads
    if(not wireframe) {
        if(mesh->triangle.size()) glDrawElements(GL_TRIANGLES, mesh->triangle.size()*3, GL_UNSIGNED_INT, &mesh->triangle[0].x);
        if(mesh->quad.size()) glDrawElements(GL_QUADS, mesh->quad.size()*4, GL_UNSIGNED_INT, &mesh->quad[0].x);
    } else {
        auto edges = EdgeMap(mesh->triangle, mesh->quad).edges();
        glDrawElements(GL_LINES, edges.size()*2, GL_UNSIGNED_INT, &edges[0].x);
    }
    
    // draw line sets
    if(not mesh->line.empty()) glDrawElements(GL_LINES, mesh->line.size()*2, GL_UNSIGNED_INT, mesh->line.data());
    for(auto segment : mesh->spline) glDrawElements(GL_LINE_STRIP, 4, GL_UNSIGNED_INT, &segment);
    
    // disable vertex attribute arrays
    glDisableVertexAttribArray(vertex_pos_location);
    glDisableVertexAttribArray(vertex_norm_location);
    if(not mesh->texcoord.empty()) glDisableVertexAttribArray(vertex_texcoord_location);
}

