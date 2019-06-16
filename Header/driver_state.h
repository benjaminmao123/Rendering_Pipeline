#ifndef __DRIVER__
#define __DRIVER__

#include "common.h"
#include <vector>
#include <utility>
#include <string>

using namespace std;

struct driver_state
{
    // Custom data that is stored per vertex, such as positions or colors.
    // These fields are stored contiguously, interleaved.  For example,
    // X Y Z R G B X Y Z R G B X Y Z R G B ...
    // Each vertex occupies floats_per_vertex entries in the array.
    // There are num_vertices vertices and thus floats_per_vertex*num_vertices
    // floats in the array.
    float * vertex_data = 0;
    int num_vertices = 0;
    int floats_per_vertex = 0;

    // If indexed rendering is being performed, this array stores the vertex
    // indices for the triangles, three ints per triangle.
    // i j k i j k i j k i j k ...
    // There are num_triangles triangles, so the array contains 3*num_triangles
    // entries.
    int * index_data = 0;
    int num_triangles = 0;

    // This is data that is constant over all triangles and fragments.
    // It is accessible from all of the shaders.  The user can store things
    // like transforms here.  The size of this array is not stored
    // since the driver will never need to know its size; you will just need.
    // to supply the pointer when necessary.
    float * uniform_data = 0;

    // Vertex data (such as color) at the vertices of triangles must be
    // interpolated to each pixel (fragment) within the triangle before calling
    // the fragment shader.  Since there are floats_per_vertex floats stored per
    // vertex, there will be floats_per_vertex valid entries in this array,
    // indicating how each float for a vertex should be interpolated.  Valid
    // values are:
    //   interp_type::flat           - each pixel receives the value stored at the first
    //                                 vertex of the triangle.
    //   interp_type::smooth         - Vertex values are interpolated using perspective-
    //                                 correct interpolation.
    //   interp_type::noperspective  - Vertex values are interpolated using image-space
    //                                 barycentric coordinates.
    interp_type interp_rules[MAX_FLOATS_PER_VERTEX] = {};

    // Image dimensions
    int image_width = 0;
    int image_height = 0;

    // Buffer where color data is stored.  The first image_width entries
    // correspond to the bottom row of the image, the next image_width entries
    // correspond to the next row, etc.  The array has image_width*image_height
    // entries.
    pixel * image_color = 0;

    // This array stores the depth of a pixel and is used for z-buffering.  The
    // size and layout is the same as image_color.
    float * image_depth = 0;

    // Pointer to a function, which performs the role of a vertex shader.  It
    // should be called on each vertex and given data stored in vertex_data.
    // This routine also receives the uniform data.
    void (*vertex_shader)(const data_vertex& in, data_geometry& out,
        const float * uniform_data);

    // Pointer to a function, which performs the role of a fragment shader.  It
    // should be called for each pixel (fragment) within each triangle.  The
    // fragment shader should be given interpolated vertex data (interpolated
    // according to interp_rules).  This routine also receives the uniform data.
    void (*fragment_shader)(const data_fragment& in, data_output& out,
        const float * uniform_data);

    driver_state();
    ~driver_state();

};

struct Edge
{
    Edge(const vec4 &start_vertex, const vec4 &end_vertex, const data_geometry &start_data, const data_geometry &end_data, const int &Index, const string &sign)
        : start_vertex(start_vertex), end_vertex(end_vertex), start_data(start_data), end_data(end_data)
    {
        if (sign == "+")
            beta = (end_vertex[3] - end_vertex[Index]) / (start_vertex[Index] - start_vertex[3] + end_vertex[3] - end_vertex[Index]);
        else       
            beta = (-end_vertex[3] - end_vertex[Index]) / (start_vertex[Index] + start_vertex[3] - end_vertex[3] - end_vertex[Index]);
    }

    vec4 start_vertex, end_vertex;
    data_geometry start_data, end_data;
    float beta;
};

// Set up the internal state of this class.  This is not done during the
// constructor since the width and height are not known when this class is
// constructed.
void initialize_render(driver_state& state, int width, int height);

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type);

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face=0);

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3]);

void RenderType_Triangle(driver_state &state, data_geometry g[3], data_vertex v[3], const data_geometry *in[3]);
void RenderType_Indexed(driver_state &state, data_geometry g[3], data_vertex v[3], const data_geometry *in[3]);
void RenderType_Fan(driver_state &state, data_geometry g[3], data_vertex v[3], const data_geometry *in[3]);
void RenderType_Strip(driver_state &state, data_geometry g[3], data_vertex v[3], const data_geometry *in[3]);

void Interpolate(driver_state &state, const Edge &edge, data_geometry &new_data, const int &Index, const data_geometry *in[3]);
vec4 Compute_Intersection(const float &beta, const vec4 &start_vertex, const vec4 &end_vertex);
vector<vector<data_geometry>> Generate_Triangles(driver_state &state, vector<Edge> &edges, const data_geometry *in[3], const int &Index);
vector<vector<data_geometry>> pCut_Triangle(driver_state &state, const int & Index, const data_geometry *in[3]);
vector<vector<data_geometry>> nCut_Triangle(driver_state &state, const int &Index, const data_geometry *in[3]);

void Interpolate(driver_state &state, const vec3 &beta, const vec3 &w, data_fragment &d_frag, const data_geometry *in[3]);
void Initialize_Bounds(driver_state &state, vec2 &min_index, vec2 &max_index);
vec3 Perform_Perspective_Divide(const vec3 &w, const int &Index, const data_geometry *in[3]);
vec2 Compute_Pixel_Coordinates(driver_state &state, const vec3 &x, const vec3 &y, const int &Index);
float Compute_Area(const vec2 &v0, const vec2 &v1, const vec2 &p);

#endif
