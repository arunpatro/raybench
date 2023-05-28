// C++ include
#include <unistd.h>

#include <Eigen/Geometry>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION  // Do not include this line twice in
                                        // your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

MatrixXd vertices_obj;  // n x 3 matrix (n points)
MatrixXi facets;        // m x 3 matrix (m triangles)
const std::string data_dir = DATA_DIR;
const std::string filename("raster.png");
// const std::string mesh_filename(data_dir + "bunny.off");
const std::string mesh_filename(data_dir + "dragon.off");

void setup_scene() {
    // Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices_obj.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i) {
        in >> vertices_obj(i, 0) >> vertices_obj(i, 1) >> vertices_obj(i, 2);
    }
    for (int i = 0; i < nf; ++i) {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }
}

// enum class ShadingMode { WIREFRAME, FLAT, VERTEX, NONE };
// enum class OutputFormat { NONE, IMAGE, GIF };
// ShadingMode shading_mode;

void render_flat(const Program& program, const UniformAttributes& uniform,
                 const MatrixXi& facets,
                 std::vector<VertexAttributes>& vertices,
                 Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic,
                               Eigen::Dynamic>& frameBuffer) {
    // flat shading, for every triangle, we add 3 vertices
    for (int i = 0; i < facets.rows(); ++i) {
        int a = facets(i, 0), b = facets(i, 1), c = facets(i, 2);

        // compute the normal of the triangle
        Eigen::Vector3d v1 = vertices_obj.row(b) - vertices_obj.row(a);
        Eigen::Vector3d v2 = vertices_obj.row(c) - vertices_obj.row(a);
        Eigen::Vector3d normal = v1.cross(v2).normalized();

        vertices.emplace_back(vertices_obj(a, 0), vertices_obj(a, 1),
                              vertices_obj(a, 2), 1, normal(0), normal(1),
                              normal(2));
        vertices.emplace_back(vertices_obj(b, 0), vertices_obj(b, 1),
                              vertices_obj(b, 2), 1, normal(0), normal(1),
                              normal(2));
        vertices.emplace_back(vertices_obj(c, 0), vertices_obj(c, 1),
                              vertices_obj(c, 2), 1, normal(0), normal(1),
                              normal(2));
    }
    //
    rasterize_triangles(program, uniform, vertices, frameBuffer);
}

// void render_vertex(const Program& program, const UniformAttributes& uniform,
//                    const MatrixXi& facets,
//                    std::vector<VertexAttributes>& vertices,
//                    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic,
//                                  Eigen::Dynamic>& frameBuffer) {
//     std::vector<Eigen::Vector3d> vertex_normals(vertices_obj.rows(),
//                                                 Eigen::Vector3d::Zero());
//     std::vector<int> num_adjacent_faces(vertices_obj.rows(), 0);
//     // vertex shading, for every vertex, we calculate the average
//     // normal of each of the three faces it belongs to, we do this
//     // by keeping track of the number of faces each vertex belongs
//     // to
//     for (int i = 0; i < facets.rows(); ++i) {
//         int a = facets(i, 0), b = facets(i, 1), c = facets(i, 2);

//         // compute the normal of the triangle
//         Eigen::Vector3d v1 = vertices_obj.row(b) - vertices_obj.row(a);
//         Eigen::Vector3d v2 = vertices_obj.row(c) - vertices_obj.row(a);
//         Eigen::Vector3d normal = v1.cross(v2).normalized();

//         // TODO correct the normal if it is pointing in the wrong
//         // direction what if there are multiple lights?

//         // add the triangle normal to each vertex's normal
//         vertex_normals[a] += normal;
//         vertex_normals[b] += normal;
//         vertex_normals[c] += normal;

//         // increment the number of adjacent faces for each vertex
//         num_adjacent_faces[a]++;
//         num_adjacent_faces[b]++;
//         num_adjacent_faces[c]++;
//     }

//     // normalize the resulting vertex normals
//     for (int i = 0; i < vertex_normals.size(); ++i) {
//         if (num_adjacent_faces[i] > 0) {
//             vertex_normals[i] /= num_adjacent_faces[i];
//             vertex_normals[i].normalize();
//         }
//     }

//     // add the vertices to the vertex buffer
//     for (int i = 0; i < facets.rows(); ++i) {
//         int a = facets(i, 0), b = facets(i, 1), c = facets(i, 2);

//         vertices.emplace_back(vertices_obj(a, 0), vertices_obj(a, 1),
//                               vertices_obj(a, 2), 1, vertex_normals[a](0),
//                               vertex_normals[a](1), vertex_normals[a](2));
//         vertices.emplace_back(vertices_obj(b, 0), vertices_obj(b, 1),
//                               vertices_obj(b, 2), 1, vertex_normals[b](0),
//                               vertex_normals[b](1), vertex_normals[b](2));
//         vertices.emplace_back(vertices_obj(c, 0), vertices_obj(c, 1),
//                               vertices_obj(c, 2), 1, vertex_normals[c](0),
//                               vertex_normals[c](1), vertex_normals[c](2));
//     }
//     rasterize_triangles(program, uniform, vertices, frameBuffer);
// }

void render_default(const Program& program, const UniformAttributes& uniform,
                    const MatrixXi& facets,
                    std::vector<VertexAttributes>& vertices,
                    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic,
                                  Eigen::Dynamic>& frameBuffer) {
    for (int i = 0; i < facets.rows(); ++i) {
        int a = facets(i, 0), b = facets(i, 1), c = facets(i, 2);

        // compute the normal of the triangle, ideally this should be same for
        // all for uniform color
        Eigen::Vector3d normal(0, 0, 1);

        vertices.emplace_back(vertices_obj(a, 0), vertices_obj(a, 1),
                              vertices_obj(a, 2), 1, normal(0), normal(1),
                              normal(2));
        vertices.emplace_back(vertices_obj(b, 0), vertices_obj(b, 1),
                              vertices_obj(b, 2), 1, normal(0), normal(1),
                              normal(2));
        vertices.emplace_back(vertices_obj(c, 0), vertices_obj(c, 1),
                              vertices_obj(c, 2), 1, normal(0), normal(1),
                              normal(2));
    }
    rasterize_triangles(program, uniform, vertices, frameBuffer);
}

// void render_shading_mode(const ShadingMode& shading_mode, Program& program,
//                          UniformAttributes& uniform, const MatrixXi& facets,
//                          vector<VertexAttributes>& vertices,
//                          FrameBuffer& frameBuffer) {
//     switch (shading_mode) {
//         case ShadingMode::WIREFRAME:
//             render_wireframe(program, uniform, vertices, frameBuffer);
//             break;
//         case ShadingMode::FLAT:
//             render_flat(program, uniform, facets, vertices, frameBuffer);
//             break;
//         case ShadingMode::VERTEX:
//             render_vertex(program, uniform, facets, vertices, frameBuffer);
//             break;
//         default:
//             render_default(program, uniform, facets, vertices, frameBuffer);
//     }
// }

void save_image(const FrameBuffer& frameBuffer, const std::string& filename) {
    vector<uint8_t> image;
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png(filename.c_str(), frameBuffer.rows(), frameBuffer.cols(), 4,
                   image.data(), frameBuffer.rows() * 4);
}

// void save_gif(const ShadingMode& shading_mode, Program& program,
//               UniformAttributes& uniform, const MatrixXi& facets,
//               vector<VertexAttributes>& vertices, FrameBuffer& frameBuffer,
//               const std::string& filename, int frames, int delay) {
//     vector<uint8_t> image;
//     GifWriter g;
//     GifBegin(&g, filename.c_str(), frameBuffer.rows(), frameBuffer.cols(),
//              delay);

//     for (int i = 0; i < frames; ++i) {
//         frameBuffer.setConstant(FrameBufferAttributes());

//         // Modify light color so that it changes over time
//         // double t = i / 50.0;            // Normalize i to [0, 1] range for
//         50
//         // frames double frequency = 2.0 * M_PI;  // Adjust the frequency to
//         // control the speed of color change double maxValue = 0.7;

//         // uniform.light_color(0) = fmod((std::sin(frequency * t) + 1) / 2 *
//         // maxValue, maxValue); uniform.light_color(1) =
//         // fmod((std::sin(frequency * t + 2 * M_PI / 3) + 1) / 2 * maxValue,
//         // maxValue); uniform.light_color(2) = fmod((std::sin(frequency * t +
//         4
//         // * M_PI / 3) + 1) / 2 * maxValue, maxValue);
//         uniform.theta_deg = i * 8 % 360;
//         uniform.camera_position(2) -= 0.07;

//         render_shading_mode(shading_mode, program, uniform, facets, vertices,
//                             frameBuffer);

//         framebuffer_to_uint8(frameBuffer, image);
//         GifWriteFrame(&g, image.data(), frameBuffer.rows(),
//         frameBuffer.cols(),
//                       delay);
//     }

//     GifEnd(&g);
// }

int main(int argc, char* argv[]) {
    setup_scene();

    // barycenter is the centroid of all the centroids of the triangles
    AlignedBox3d box = AlignedBox3d();
    Vector3d barycenter(0, 0, 0);
    for (int i = 0; i < facets.rows(); ++i) {
        Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 3; ++j) {
            centroid += vertices_obj.row(facets(i, j)).transpose();
            box.extend(vertices_obj.row(facets(i, j)).transpose());
        }
        centroid /= 3;
        barycenter += centroid;
    }
    barycenter /= facets.rows();

    // std::cout << "bbox in world space" << std::endl;
    // std::cout << "min: " << box.min().transpose() << std::endl;
    // std::cout << "max: " << box.max().transpose() << std::endl;
    // std::cout << "w: " << box.max()(0) - box.min()(0) << std::endl;
    // std::cout << "h: " << box.max()(1) - box.min()(1) << std::endl;
    // std::cout << "d: " << box.max()(2) - box.min()(2) << std::endl;
    // std::cout << "center: " << box.center().transpose() << std::endl;
    // std::cout << "barycenter: " << barycenter.transpose() << std::endl;

    // The Framebuffer storing the image rendered by the rasterizer
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>
        frameBuffer(500, 500);

    // Global Constants (empty in this example)
    UniformAttributes uniform;
    uniform.box = box;
    uniform.ambient_color = Vector3d(0.2, 0.2, 0.2);
    uniform.light_position = Vector3d(-1, 1, 3);
    uniform.light_color = Vector3d(0.2, 0.5, 0.1);
    uniform.barycenter = barycenter;
    uniform.theta_deg = 0;

    // align the camera position with the barycenter of the object
    uniform.camera_position = Vector3d(0, 0, 2);

    // this is very minimal effect
    // uniform.camera_position = Vector3d(box.center()(0), box.center()(1), 2);
    // uniform.camera_position = Vector3d(-0.0170336, 0.110509, 1.07);

    uniform.camera_up = Vector3d(0, 1, 0);
    uniform.camera_direction = Vector3d(0, 0, -1);
    uniform.camera_right = uniform.camera_direction.cross(uniform.camera_up);

    uniform.camera_focal_length = 1;
    uniform.camera_fov = 0.8;
    uniform.camera_near = -uniform.camera_focal_length;
    uniform.camera_far = -10;
    uniform.camera_aspect = frameBuffer.cols() / (double)frameBuffer.rows();

    uniform.is_perspective = false;

    Program program;

    // The vertex shader is the identity
    program.VertexShader = [](const VertexAttributes& va,
                              const UniformAttributes& uniform) {
        // AlignedBox3d box = uniform.box;

        // // Add the aspect ratio transform matrix
        // Eigen::Matrix4d aspect_ratio_transform;
        Eigen::Matrix4d model_matrix = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d camera_matrix;
        Eigen::Matrix4d orthographic_matrix;
        // Eigen::Matrix4d perspective_matrix;
        // Eigen::Matrix4d translation_matrix, translation_matrix_back,
        //     rotation_matrix, model_matrix;

        // aspect_ratio_transform.setIdentity();

        // if (uniform.camera_aspect > 1) {
        //     aspect_ratio_transform(1, 1) = 1 / uniform.camera_aspect;
        // } else {
        //     aspect_ratio_transform(0, 0) = uniform.camera_aspect;
        // }

        // // this rotation is only about the z axis, if theta is 0, then model
        // // matrix is identity
        // double theta = uniform.theta_deg * M_PI / 180;
        // translation_matrix << 1, 0, 0, -uniform.barycenter(0), 0, 1, 0,
        //     -uniform.barycenter(1), 0, 0, 1, -uniform.barycenter(2), 0, 0, 0,
        //     1;
        // rotation_matrix << cos(theta), 0, sin(theta), 0, 0, 1, 0, 0,
        //     -sin(theta), 0, cos(theta), 0, 0, 0, 0, 1;
        // translation_matrix_back << 1, 0, 0, uniform.barycenter(0), 0, 1, 0,
        //     uniform.barycenter(1), 0, 0, 1, uniform.barycenter(2), 0, 0, 0,
        //     1;

        // model_matrix =
        //     translation_matrix_back * rotation_matrix * translation_matrix;

        // // to correct aspect ratio, we must apply stretch the object in
        // // camera/world space, as well as contract the perspective volume we
        // // first convert from world space to object space apply the aspect
        // ratio
        // // transform and then convert back to world space
        // Eigen::Matrix4d trans_aspect, trans_aspect_back, model_aspect_tfx;
        // trans_aspect << 1, 0, 0, -box.center()(0), 0, 1, 0, -box.center()(1),
        // 0,
        //     0, 1, -box.center()(2), 0, 0, 0, 1;
        // trans_aspect_back << 1, 0, 0, box.center()(0), 0, 1, 0,
        // box.center()(1),
        //     0, 0, 1, box.center()(2), 0, 0, 0, 1;
        // model_aspect_tfx =
        //     trans_aspect_back * aspect_ratio_transform * trans_aspect;

        // model_matrix = model_matrix * model_aspect_tfx;
        // camera matrix
        camera_matrix << 1, 0, 0, -uniform.camera_position.x(), 0, 1, 0,
            -uniform.camera_position.y(), 0, 0, 1, -uniform.camera_position.z(),
            0, 0, 0, 1;

        double top = uniform.camera_focal_length * tan(uniform.camera_fov / 2);
        double bottom = -top;
        double right = top * uniform.camera_aspect;
        double left = -right;
        double near = -uniform.camera_focal_length;
        double far = -10;

        // double top = 1., bottom = -1., right = 1., left = -1., near = -1.,
        //        far = -10.;

        orthographic_matrix << 2 / (right - left), 0, 0,
            -(right + left) / (right - left), 0, 2 / (top - bottom), 0,
            -(top + bottom) / (top - bottom), 0, 0, 2 / (near - far),
            -(near + far) / (near - far), 0, 0, 0, 1;

        // if (uniform.is_perspective) {
        //     perspective_transform_matrix << uniform.camera_near, 0, 0, 0, 0,
        //         uniform.camera_near, 0, 0, 0, 0,
        //         (uniform.camera_near + uniform.camera_far),
        //         -uniform.camera_near * uniform.camera_far, 0, 0, 1, 0;
        // }

        Eigen::Vector4d vertex(va.position[0], va.position[1], va.position[2],
                               1.0);

        // Eigen::Matrix4d final_transform = Matrix4d::Identity();
        // Eigen::Matrix4d final_transform = camera_matrix *
        // Matrix4d::Identity();

        Eigen::Matrix4d final_transform =
            orthographic_matrix * camera_matrix * model_matrix;

        // Eigen::Matrix4d final_transform = frustrum_ortho *
        // perspective_transform_matrix * camera_matrix;

        // object which is present in frustrum gets projected to canonical box
        // Eigen::Matrix4d final_transform = object_ortho * camera_matrix;

        // transform the normals to camera space
        Eigen::Vector4d projected_normal(va.normal[0], va.normal[1],
                                         va.normal[2], 0.0);
        // projected_normal = model_matrix * projected_normal;

        // this is wrong but this gives a mettalic look why?
        // projected_normal = final_transform * projected_normal;

        // project the vertex to the canonical box
        Eigen::Vector4d projected_vertex = final_transform * vertex;

        projected_vertex /= projected_vertex.w();

        return VertexAttributes(projected_vertex.x(), projected_vertex.y(),
                                projected_vertex.z(), 1., projected_normal.x(),
                                projected_normal.y(), projected_normal.z());
    };

    // The fragment shader uses a fixed color
    // TODO use z-buffer to determine color
    program.FragmentShader = [](const VertexAttributes& va,
                                const UniformAttributes& uniform) {
        // calculate the diffuse color
        Vector3d new_light_position = uniform.light_position;

        Eigen::Vector3d light_direction =
            new_light_position -
            Vector3d(va.position[0], va.position[1], va.position[2]);
        light_direction.normalize();
        Eigen::Vector3d normal = va.normal;
        normal.normalize();
        if (normal.dot(light_direction) < 0) {
            normal = -normal;
        }
        double diffuse = std::max(0.0, normal.dot(light_direction));
        Eigen::Vector3d diffuse_color = diffuse * uniform.light_color;
        auto color = diffuse_color + uniform.ambient_color;

        auto frag = FragmentAttributes(color[0], color[1], color[2], 1.);
        frag.position = va.position;
        return frag;
    };

    // The blending shader converts colors between 0 and 1 to uint8
    program.BlendingShader = [](const FragmentAttributes& fa,
                                const FrameBufferAttributes& previous) {
        auto fba = FrameBufferAttributes(fa.color[0] * 255, fa.color[1] * 255,
                                         fa.color[2] * 255, fa.color[3] * 255);
        if (fa.position(2) > previous.depth) {
            fba.depth = fa.position(2);
            return fba;
        } else {
            return previous;
        }
    };

    // Add the vertices of the mesh by exploding the facets
    vector<VertexAttributes> vertices;

    render_flat(program, uniform, facets, vertices, frameBuffer);
    save_image(frameBuffer, "triangle.png");
    // // Save the output based on the output_format
    // if (output_format == OutputFormat::IMAGE) {
    //     uniform.is_perspective = false;
    //     uniform.camera_position =
    //         Vector3d(box.center()(0), box.center()(1), 1.06);
    //     render_shading_mode(shading_mode, program, uniform, facets, vertices,
    //                         frameBuffer);
    //     save_image(frameBuffer, "triangle.png");
    // } else if (output_format == OutputFormat::GIF) {
    //     uniform.is_perspective = true;
    //     uniform.camera_position =
    //         Vector3d(box.center()(0), box.center()(1), 3.06);
    //     save_gif(shading_mode, program, uniform, facets, vertices,
    //     frameBuffer,
    //              "triangle.gif", 10, 50);
    // } else {
    //     // Handle the case when no output format is specified, if necessary
    // }

    return 0;
}