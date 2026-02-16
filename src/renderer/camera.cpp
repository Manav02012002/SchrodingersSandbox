#include "renderer/camera.h"

#include <algorithm>
#include <cmath>

namespace sbox {

namespace {

constexpr float kPi = 3.14159265358979323846f;

}  // namespace

Camera::Camera()
    : orientation_(Eigen::Quaternionf::Identity()),
      target_(Eigen::Vector3f::Zero()),
      distance_(15.0f),
      viewport_width_(1600.0f),
      viewport_height_(1000.0f) {}

void Camera::setViewportSize(float width, float height) {
    viewport_width_ = std::max(width, 1.0f);
    viewport_height_ = std::max(height, 1.0f);
}

void Camera::handleRotationDrag(float dx, float dy, bool apply) {
    if (!apply) {
        return;
    }

    const float sensitivity = 0.005f;
    const Eigen::Quaternionf yaw(Eigen::AngleAxisf(-dx * sensitivity, Eigen::Vector3f::UnitY()));
    const Eigen::Quaternionf pitch(Eigen::AngleAxisf(-dy * sensitivity, Eigen::Vector3f::UnitX()));
    orientation_ = (yaw * orientation_ * pitch).normalized();
}

void Camera::handlePanDrag(float dx, float dy, bool apply) {
    if (!apply) {
        return;
    }

    const float pan_scale = distance_ * 0.0015f;
    const Eigen::Vector3f right = orientation_ * Eigen::Vector3f::UnitX();
    const Eigen::Vector3f up = orientation_ * Eigen::Vector3f::UnitY();
    target_ -= right * (dx * pan_scale);
    target_ += up * (dy * pan_scale);
}

void Camera::handleScroll(float yoffset, bool apply) {
    if (!apply) {
        return;
    }

    distance_ *= (1.0f - yoffset * 0.1f);
    distance_ = std::clamp(distance_, 0.5f, 200.0f);
}

Eigen::Vector3f Camera::cameraPosition() const {
    return target_ + orientation_ * Eigen::Vector3f(0.0f, 0.0f, distance_);
}

Eigen::Vector3f Camera::camera_position() const {
    return cameraPosition();
}

Eigen::Matrix4f Camera::viewMatrix() const {
    const Eigen::Vector3f eye = cameraPosition();
    const Eigen::Vector3f up = orientation_ * Eigen::Vector3f::UnitY();
    return LookAt(eye, target_, up);
}

Eigen::Matrix4f Camera::projectionMatrix() const {
    const float aspect = viewport_width_ / viewport_height_;
    return Perspective(45.0f * kPi / 180.0f, aspect, 0.1f, 1000.0f);
}

Eigen::Matrix4f Camera::inverseViewProjection() const {
    return (projectionMatrix() * viewMatrix()).inverse();
}

Eigen::Matrix4f Camera::inv_view_projection() const {
    return inverseViewProjection();
}

Eigen::Matrix4f Camera::LookAt(const Eigen::Vector3f& eye,
                               const Eigen::Vector3f& center,
                               const Eigen::Vector3f& up) {
    const Eigen::Vector3f f = (center - eye).normalized();
    const Eigen::Vector3f s = f.cross(up).normalized();
    const Eigen::Vector3f u = s.cross(f);

    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    mat(0, 0) = s.x();
    mat(0, 1) = s.y();
    mat(0, 2) = s.z();
    mat(1, 0) = u.x();
    mat(1, 1) = u.y();
    mat(1, 2) = u.z();
    mat(2, 0) = -f.x();
    mat(2, 1) = -f.y();
    mat(2, 2) = -f.z();
    mat(0, 3) = -s.dot(eye);
    mat(1, 3) = -u.dot(eye);
    mat(2, 3) = f.dot(eye);
    return mat;
}

Eigen::Matrix4f Camera::Perspective(float fovy_radians, float aspect, float z_near, float z_far) {
    const float tan_half = std::tan(fovy_radians * 0.5f);

    Eigen::Matrix4f mat = Eigen::Matrix4f::Zero();
    mat(0, 0) = 1.0f / (aspect * tan_half);
    mat(1, 1) = 1.0f / tan_half;
    mat(2, 2) = -(z_far + z_near) / (z_far - z_near);
    mat(2, 3) = -(2.0f * z_far * z_near) / (z_far - z_near);
    mat(3, 2) = -1.0f;
    return mat;
}

}  // namespace sbox
