#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

class VertexAttributes {
   public:
    VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1,
                     double nx = 0, double ny = 0, double nz = 0) {
        position << x, y, z, w;
        normal = Eigen::Vector3d(nx, ny, nz);
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(const VertexAttributes& a,
                                        const VertexAttributes& b,
                                        const VertexAttributes& c,
                                        const double alpha, const double beta,
                                        const double gamma) {
        VertexAttributes r;
        r.position =
            alpha * a.position + beta * b.position + gamma * c.position;

        r.normal = alpha * a.normal + beta * b.normal + gamma * c.normal;
        return r;
    }

    Eigen::Vector4d position;
    Eigen::Vector3d normal;
};

class FragmentAttributes {
   public:
    FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1) {
        color << r, g, b, a;
    }

    Eigen::Vector4d color;
    Eigen::Vector4d position;
};

class FrameBufferAttributes {
   public:
    FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0,
                          uint8_t a = 255) {
        color << r, g, b, a;
        // neg inf in double
        depth = std::numeric_limits<double>::lowest();
    }

    Eigen::Matrix<uint8_t, 4, 1> color;
    double depth;
};

class UniformAttributes {
   public:
    Eigen::AlignedBox3d box;
    Eigen::Vector3d ambient_color;
    Eigen::Vector3d light_color;
    Eigen::Vector3d light_position;
    Eigen::Vector3d barycenter;
    double theta_deg;

    Eigen::Vector3d camera_position;
    Eigen::Vector3d camera_up;
    Eigen::Vector3d camera_direction;
    Eigen::Vector3d camera_right;

    double camera_focal_length;
    double camera_fov;
    double camera_near;
    double camera_far;
    double camera_aspect;

    bool is_perspective;
};