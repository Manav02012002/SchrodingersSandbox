#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace sbox {

class Camera {
public:
    Camera();

    void setViewportSize(float width, float height);
    void handleRotationDrag(float dx, float dy, bool apply);
    void handlePanDrag(float dx, float dy, bool apply);
    void handleScroll(float yoffset, bool apply);

    [[nodiscard]] Eigen::Matrix4f viewMatrix() const;
    [[nodiscard]] Eigen::Matrix4f projectionMatrix() const;
    [[nodiscard]] Eigen::Matrix4f inverseViewProjection() const;

    [[nodiscard]] Eigen::Vector3f cameraPosition() const;
    [[nodiscard]] Eigen::Vector3f camera_position() const;
    [[nodiscard]] Eigen::Matrix4f inv_view_projection() const;

private:
    static Eigen::Matrix4f LookAt(const Eigen::Vector3f& eye,
                                  const Eigen::Vector3f& center,
                                  const Eigen::Vector3f& up);
    static Eigen::Matrix4f Perspective(float fovy_radians, float aspect, float z_near, float z_far);

    Eigen::Quaternionf orientation_;
    Eigen::Vector3f target_;
    float distance_;
    float viewport_width_;
    float viewport_height_;
};

}  // namespace sbox
