#pragma once

#include <Eigen/Core>
#include <imgui.h>

#include <string>
#include <vector>

namespace sbox::ui {

struct Annotation {
    int id = 0;

    enum class Type { Text, Arrow, DimensionLine, AngleArc, Circle };
    Type type = Type::Text;

    Eigen::Vector3d world_pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d world_pos_end = Eigen::Vector3d::Zero();
    Eigen::Vector3d world_pos_aux = Eigen::Vector3d::Zero();

    std::string text;
    float font_scale = 1.0f;
    ImVec4 color = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
    float thickness = 1.5f;
    bool visible = true;

    bool auto_value = false;
    std::string unit = "\xC3\x85";
};

class AnnotationManager {
public:
    int add_text(const Eigen::Vector3d& pos, const std::string& text);
    int add_arrow(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const std::string& label = "");
    int add_dimension_line(const Eigen::Vector3d& start, const Eigen::Vector3d& end, bool auto_value = true);
    int add_angle_arc(const Eigen::Vector3d& vertex, const Eigen::Vector3d& a, const Eigen::Vector3d& b, bool auto_value = true);
    int add_circle(const Eigen::Vector3d& center, const Eigen::Vector3d& radius_point, const std::string& label = "");

    void remove(int id);
    void clear();

    Annotation& get(int id);
    const std::vector<Annotation>& all() const;

    void draw(const Eigen::Matrix4f& vp_matrix, const ImVec2& viewport_pos, const ImVec2& viewport_size);

    void set_default_style(const ImVec4& color, float font_scale, float thickness);
    [[nodiscard]] ImVec4 default_color() const;
    [[nodiscard]] float default_font_scale() const;
    [[nodiscard]] float default_thickness() const;

private:
    Annotation& append(Annotation annotation);
    std::vector<Annotation> annotations_;
    int next_id_ = 1;
    ImVec4 default_color_ = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
    float default_font_scale_ = 1.0f;
    float default_thickness_ = 1.5f;
};

}  // namespace sbox::ui
