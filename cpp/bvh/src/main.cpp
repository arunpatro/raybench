////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stack>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION  // Do not include this line twice in
                                        // your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////
class AABBNode {
   public:
    AlignedBox3d box = AlignedBox3d();
    AABBNode *left = nullptr;
    AABBNode *right = nullptr;
    AABBNode *parent = nullptr;
    int triangle_index = -1;
};

AABBNode *recur_tree(const MatrixXd &V, const MatrixXi &F,
                     const std::vector<int> triangle_indices, AABBNode *parent);

class AABBTree {
   public:
    class Node {
       public:
        AlignedBox3d bbox;
        int parent;    // Index of the parent node (-1 for root)
        int left;      // Index of the left child (-1 for a leaf)
        int right;     // Index of the right child (-1 for a leaf)
        int triangle;  // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
    int root;
    AABBNode *root_node = nullptr;

    AABBTree() = default;  // Default empty constructor
    AABBTree(const MatrixXd &V,
             const MatrixXi &F);  // Build a BVH from an existing mesh
};

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string data_dir = DATA_DIR;
const std::string filename("bvh.png");
// const std::string mesh_filename(data_dir + "dodeca.off");
// const std::string mesh_filename(data_dir + "cube.off");
// const std::string mesh_filename(data_dir + "bunny.off");
const std::string mesh_filename(data_dir + "dragon.off");

// Camera settings
const double focal_length = 5;
const double field_of_view = 0.3491;  // 20 degrees
// const double field_of_view = 0.7854;  // 45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);

// Triangle Mesh
MatrixXd vertices;  // n x 3 matrix (n points)
MatrixXi facets;    // m x 3 matrix (m triangles)
AABBTree bvh;

// Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

// Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
// Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

// Fills the different arrays
void setup_scene() {
    // Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i) {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i) {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    // setup tree
    bvh = AABBTree(vertices, facets);

    // Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);
}

////////////////////////////////////////////////////////////////////////////////
// BVH Code
////////////////////////////////////////////////////////////////////////////////

AlignedBox3d bbox_from_triangle(const Vector3d &a, const Vector3d &b,
                                const Vector3d &c) {
    AlignedBox3d box;
    box.extend(a);
    box.extend(b);
    box.extend(c);
    return box;
}

AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F) {
    // Compute the centroids of all the triangles in the input mesh
    MatrixXd centroids(F.rows(), V.cols());
    centroids.setZero();
    for (int i = 0; i < F.rows(); ++i) {
        for (int k = 0; k < F.cols(); ++k) {
            centroids.row(i) += V.row(F(i, k));
        }
        centroids.row(i) /= F.cols();
    }

    std::vector<int> triangle_indices;
    for (int i = 0; i < F.rows(); ++i) {
        triangle_indices.push_back(i);
    }

    AlignedBox3d box;
    for (int i = 0; i < F.rows(); ++i) {
        box.extend(centroids.row(i).transpose());
    }

    root_node = recur_tree(V, F, triangle_indices, nullptr);
}

AABBNode *recur_tree(const MatrixXd &V, const MatrixXi &F,
                     const std::vector<int> triangle_indices,
                     AABBNode *parent) {
    if (triangle_indices.size() == 0) {
        return nullptr;
    } else if (triangle_indices.size() == 1) {
        auto node = new AABBNode();
        node->parent = parent;
        node->left = nullptr;
        node->right = nullptr;

        node->box = bbox_from_triangle(V.row(F(triangle_indices[0], 0)),
                                       V.row(F(triangle_indices[0], 1)),
                                       V.row(F(triangle_indices[0], 2)));
        node->triangle_index = triangle_indices[0];
        return node;
    }

    auto node = new AABBNode();
    node->parent = parent;

    for (int i = 0; i < triangle_indices.size(); ++i) {
        node->box.extend(bbox_from_triangle(V.row(F(triangle_indices[i], 0)),
                                            V.row(F(triangle_indices[i], 1)),
                                            V.row(F(triangle_indices[i], 2))));
    }

    std::vector<int> left_indices;
    std::vector<int> right_indices;

    int axis = 0;
    if (node->box.sizes()(1) > node->box.sizes()(0)) {
        axis = 1;
    }
    if (node->box.sizes()(2) > node->box.sizes()(axis)) {
        axis = 2;
    }

    // std::sort(
    //     triangle_indices.begin(), triangle_indices.end(),
    //     [&V, &F, axis](int a, int b) {
    //         return (V.row(F(a, 0)) + V.row(F(a, 1)) + V.row(F(a, 2)))(axis) <
    //                (V.row(F(b, 0)) + V.row(F(b, 1)) + V.row(F(b, 2)))(axis);
    //     });

    for (int i = 0; i < triangle_indices.size() / 2; ++i) {
        left_indices.push_back(triangle_indices[i]);
    }
    for (int i = triangle_indices.size() / 2; i < triangle_indices.size();
         ++i) {
        right_indices.push_back(triangle_indices[i]);
    }

    node->left = recur_tree(V, F, left_indices, node);
    node->right = recur_tree(V, F, right_indices, node);
    node->triangle_index = -1;

    return node;
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

double ray_triangle_intersection(const Vector3d &ray_origin,
                                 const Vector3d &ray_direction,
                                 const Vector3d &a, const Vector3d &b,
                                 const Vector3d &c, Vector3d &p, Vector3d &N) {
    // TODO
    // Compute whether the ray intersects the given triangle.
    // If you have done the parallelogram case, this should be very similar to
    // it.

    auto edge1 = b - a;
    auto edge2 = c - a;

    Vector3d bsystem = a - ray_origin;
    Matrix3d Asystem;
    Asystem << -edge1, -edge2, ray_direction;

    Vector3d uvt = Asystem.colPivHouseholderQr().solve(bsystem);
    double u = uvt(0);
    double v = uvt(1);
    double t = uvt(2);

    if (t < 0 || u < 0 || v < 0 || u + v >= 1) {
        return -1;
    }

    Vector3d ray_intersection = ray_origin + t * ray_direction;
    Vector3d ray_normal = edge1.cross(edge2).normalized();
    // normal should point towards the camera
    if (ray_normal.dot(ray_direction) > 0) {
        ray_normal = -ray_normal;
    }

    p = ray_intersection;
    N = ray_normal;

    return t;
}

bool ray_box_intersection(const Vector3d &ray_origin,
                          const Vector3d &ray_direction,
                          const AlignedBox3d &box) {
    // TODO
    // Compute whether the ray intersects the given box.
    // we are not testing with the real surface here anyway.

    Vector3d inv_dir(1.0 / ray_direction.x(), 1.0 / ray_direction.y(),
                     1.0 / ray_direction.z());

    double tx1 = (box.min().x() - ray_origin.x()) * inv_dir.x();
    double tx2 = (box.max().x() - ray_origin.x()) * inv_dir.x();
    double ty1 = (box.min().y() - ray_origin.y()) * inv_dir.y();
    double ty2 = (box.max().y() - ray_origin.y()) * inv_dir.y();
    double tz1 = (box.min().z() - ray_origin.z()) * inv_dir.z();
    double tz2 = (box.max().z() - ray_origin.z()) * inv_dir.z();

    double tmin = std::max(std::max(std::min(tx1, tx2), std::min(ty1, ty2)),
                           std::min(tz1, tz2));
    double tmax = std::min(std::min(std::max(tx1, tx2), std::max(ty1, ty2)),
                           std::max(tz1, tz2));

    if (tmax < 0) {
        return false;
    }

    if (tmin > tmax) {
        return false;
    }

    return true;
}

// Finds the closest intersecting object returns its index
// In case of intersection it writes into p and N (intersection point and
// normals)
bool find_nearest_object(const Vector3d &ray_origin,
                         const Vector3d &ray_direction, Vector3d &p,
                         Vector3d &N) {
    Vector3d tmp_p, tmp_N;
    double min_t = std::numeric_limits<double>::max();
    bool found_intersection = false;

    // TODO
    // Method (1): Traverse every triangle and return the closest hit.

    // for (int i = 0; i < facets.rows(); ++i) {
    //     Vector3d a = vertices.row(facets(i, 0));
    //     Vector3d b = vertices.row(facets(i, 1));
    //     Vector3d c = vertices.row(facets(i, 2));

    //     double t = ray_triangle_intersection(ray_origin, ray_direction, a, b,
    //     c,
    //                                          tmp_p, tmp_N);
    //     if (t > 0 && t < min_t) {
    //         min_t = t;
    //         p = tmp_p;
    //         N = tmp_N;
    //         found_intersection = true;
    //     }
    // }

    // Method (2): Traverse the BVH tree and test the intersection with a
    // triangles at the leaf nodes that intersects the input ray.

    std::stack<AABBNode *> stack;
    stack.push(bvh.root_node);

    while (!stack.empty()) {
        AABBNode *node = stack.top();
        stack.pop();

        if (ray_box_intersection(ray_origin, ray_direction, node->box)) {
            if (node->left == nullptr && node->right == nullptr) {
                double t = ray_triangle_intersection(
                    ray_origin, ray_direction,
                    vertices.row(facets(node->triangle_index, 0)),
                    vertices.row(facets(node->triangle_index, 1)),
                    vertices.row(facets(node->triangle_index, 2)), tmp_p,
                    tmp_N);
                if (t > 0 && t < min_t) {
                    min_t = t;
                    p = tmp_p;
                    N = tmp_N;
                    found_intersection = true;
                }
            } else {
                stack.push(node->left);
                stack.push(node->right);
            }
        }
    }

    return found_intersection;
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction) {
    // Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object =
        find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object) {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ambient light contribution
    const Vector4d ambient_color =
        obj_ambient_color.array() * ambient_light.array();

    // Punctual lights contribution (direct lighting)
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < light_positions.size(); ++i) {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        Vector4d diff_color = obj_diffuse_color;

        // Diffuse contribution
        const Vector3d Li = (light_position - p).normalized();
        const Vector4d diffuse = diff_color * std::max(Li.dot(N), 0.0);

        // Specular contribution
        const Vector3d Hi = (Li - ray_direction).normalized();
        const Vector4d specular =
            obj_specular_color *
            std::pow(std::max(N.dot(Hi), 0.0), obj_specular_exponent);
        // Vector3d specular(0, 0, 0);

        // Attenuate lights according to the squared distance to the lights
        const Vector3d D = light_position - p;
        lights_color +=
            (diffuse + specular).cwiseProduct(light_color) / D.squaredNorm();
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color;

    // Set alpha to 1
    C(3) = 1;

    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene() {
    std::cout << "Simple ray tracer." << std::endl;

    int w = 640;
    int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h);  // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    // TODO
    double image_y = 2 * tan(field_of_view / 2.0) * focal_length;
    double image_x = image_y * aspect_ratio;

    // The pixel grid through which we shoot rays is at a distance
    // 'focal_length'
    const Vector3d image_origin(-image_x, image_y,
                                camera_position[2] - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i) {
        for (unsigned j = 0; j < h; ++j) {
            const Vector3d pixel_center = image_origin +
                                          (i + 0.5) * x_displacement +
                                          (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (is_perspective) {
                // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            } else {
                // Orthographic camera
                ray_origin = pixel_center;
                ray_direction = Vector3d(0, 0, -1);
            }

            const Vector4d C = shoot_ray(ray_origin, ray_direction);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    setup_scene();

    raytrace_scene();
    return 0;
}