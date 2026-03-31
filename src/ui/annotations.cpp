#include "ui/annotations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <stdexcept>

namespace sbox::ui {

namespace {

constexpr double kBohrToAngstrom = 0.529177210903;
constexpr double kRadToDeg = 57.29577951308232;

ImU32 to_u32(const ImVec4& color) {
    return ImGui::ColorConvertFloat4ToU32(color);
}

ImVec4 background_for(const ImVec4& color) {
    const float brightness = 0.2126f * color.x + 0.7152f * color.y + 0.0722f * color.z;
    return brightness > 0.55f ? ImVec4(0.02f, 0.02f, 0.02f, 0.62f) : ImVec4(0.94f, 0.94f, 0.94f, 0.28f);
}

bool world_to_screen(const Eigen::Vector3f& world_pos,
                     const Eigen::Matrix4f& vp_matrix,
                     const ImVec2& viewport_pos,
                     const ImVec2& viewport_size,
                     ImVec2& screen_pos) {
    const Eigen::Vector4f clip = vp_matrix * Eigen::Vector4f(world_pos.x(), world_pos.y(), world_pos.z(), 1.0f);
    if (clip.w() <= 0.0f) {
        return false;
    }
    const Eigen::Vector3f ndc = clip.head<3>() / clip.w();
    if (ndc.x() < -1.0f || ndc.x() > 1.0f || ndc.y() < -1.0f || ndc.y() > 1.0f || ndc.z() < -1.0f || ndc.z() > 1.0f) {
        return false;
    }
    screen_pos.x = viewport_pos.x + (ndc.x() * 0.5f + 0.5f) * viewport_size.x;
    screen_pos.y = viewport_pos.y + (1.0f - (ndc.y() * 0.5f + 0.5f)) * viewport_size.y;
    return true;
}

void draw_text_box(ImDrawList* draw_list,
                   const ImVec2& pos,
                   const std::string& text,
                   float font_scale,
                   const ImVec4& color) {
    if (text.empty()) {
        return;
    }
    ImFont* font = ImGui::GetFont();
    const float font_size = ImGui::GetFontSize() * std::max(font_scale, 0.1f);
    const ImVec2 text_size = font->CalcTextSizeA(font_size, FLT_MAX, 0.0f, text.c_str());
    const ImVec2 padding(6.0f, 3.0f);
    draw_list->AddRectFilled(ImVec2(pos.x - padding.x, pos.y - padding.y),
                             ImVec2(pos.x + text_size.x + padding.x, pos.y + text_size.y + padding.y),
                             to_u32(background_for(color)),
                             4.0f);
    draw_list->AddText(font, font_size, pos, to_u32(color), text.c_str());
}

void draw_arrowhead(ImDrawList* draw_list, const ImVec2& tip, const ImVec2& tail, ImU32 color) {
    const ImVec2 delta(tip.x - tail.x, tip.y - tail.y);
    const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
    if (length < 1.0f) {
        return;
    }
    const ImVec2 dir(delta.x / length, delta.y / length);
    const ImVec2 perp(-dir.y, dir.x);
    const float head_length = 12.0f;
    const float head_width = 5.5f;
    const ImVec2 base(tip.x - dir.x * head_length, tip.y - dir.y * head_length);
    draw_list->AddTriangleFilled(tip,
                                 ImVec2(base.x + perp.x * head_width, base.y + perp.y * head_width),
                                 ImVec2(base.x - perp.x * head_width, base.y - perp.y * head_width),
                                 color);
}

std::string format_distance(const Annotation& annotation) {
    const double distance_bohr = (annotation.world_pos_end - annotation.world_pos).norm();
    const double value = annotation.unit == "\xC3\x85" ? distance_bohr * kBohrToAngstrom : distance_bohr;
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f %s", value, annotation.unit.c_str());
    return buffer;
}

std::string format_angle(const Eigen::Vector3d& vertex,
                         const Eigen::Vector3d& a,
                         const Eigen::Vector3d& b) {
    const Eigen::Vector3d va = (a - vertex).normalized();
    const Eigen::Vector3d vb = (b - vertex).normalized();
    const double angle_deg = std::acos(std::clamp(va.dot(vb), -1.0, 1.0)) * kRadToDeg;
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.1f\xC2\xB0", angle_deg);
    return buffer;
}

}  // namespace

Annotation& AnnotationManager::append(Annotation annotation) {
    annotation.id = next_id_++;
    annotation.color = default_color_;
    annotation.font_scale = default_font_scale_;
    annotation.thickness = default_thickness_;
    annotations_.push_back(std::move(annotation));
    return annotations_.back();
}

int AnnotationManager::add_text(const Eigen::Vector3d& pos, const std::string& text) {
    Annotation annotation;
    annotation.type = Annotation::Type::Text;
    annotation.world_pos = pos;
    annotation.text = text;
    return append(std::move(annotation)).id;
}

int AnnotationManager::add_arrow(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const std::string& label) {
    Annotation annotation;
    annotation.type = Annotation::Type::Arrow;
    annotation.world_pos = start;
    annotation.world_pos_end = end;
    annotation.text = label;
    return append(std::move(annotation)).id;
}

int AnnotationManager::add_dimension_line(const Eigen::Vector3d& start, const Eigen::Vector3d& end, bool auto_value) {
    Annotation annotation;
    annotation.type = Annotation::Type::DimensionLine;
    annotation.world_pos = start;
    annotation.world_pos_end = end;
    annotation.auto_value = auto_value;
    return append(std::move(annotation)).id;
}

int AnnotationManager::add_angle_arc(const Eigen::Vector3d& vertex,
                                     const Eigen::Vector3d& a,
                                     const Eigen::Vector3d& b,
                                     bool auto_value) {
    Annotation annotation;
    annotation.type = Annotation::Type::AngleArc;
    annotation.world_pos = vertex;
    annotation.world_pos_end = a;
    annotation.world_pos_aux = b;
    annotation.auto_value = auto_value;
    annotation.unit = "deg";
    Annotation& stored = append(std::move(annotation));
    stored.text = auto_value ? format_angle(vertex, a, b) : "";
    stored.thickness = std::max(stored.thickness, 2.0f);
    stored.visible = true;
    return stored.id;
}

int AnnotationManager::add_circle(const Eigen::Vector3d& center, const Eigen::Vector3d& radius_point, const std::string& label) {
    Annotation annotation;
    annotation.type = Annotation::Type::Circle;
    annotation.world_pos = center;
    annotation.world_pos_end = radius_point;
    annotation.text = label;
    return append(std::move(annotation)).id;
}

void AnnotationManager::remove(int id) {
    annotations_.erase(std::remove_if(annotations_.begin(), annotations_.end(), [&](const Annotation& a) {
                           return a.id == id;
                       }),
                       annotations_.end());
}

void AnnotationManager::clear() {
    annotations_.clear();
}

Annotation& AnnotationManager::get(int id) {
    auto it = std::find_if(annotations_.begin(), annotations_.end(), [&](const Annotation& a) {
        return a.id == id;
    });
    if (it == annotations_.end()) {
        throw std::out_of_range("annotation id not found");
    }
    return *it;
}

const std::vector<Annotation>& AnnotationManager::all() const {
    return annotations_;
}

void AnnotationManager::draw(const Eigen::Matrix4f& vp_matrix,
                             const ImVec2& viewport_pos,
                             const ImVec2& viewport_size) {
    if (viewport_size.x <= 1.0f || viewport_size.y <= 1.0f) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    for (const Annotation& annotation : annotations_) {
        if (!annotation.visible) {
            continue;
        }

        const ImU32 color = to_u32(annotation.color);

        if (annotation.type == Annotation::Type::Text) {
            ImVec2 screen{};
            if (!world_to_screen(annotation.world_pos.cast<float>(), vp_matrix, viewport_pos, viewport_size, screen)) {
                continue;
            }
            draw_text_box(draw_list, ImVec2(screen.x + 6.0f, screen.y - 20.0f), annotation.text, annotation.font_scale, annotation.color);
            continue;
        }

        if (annotation.type == Annotation::Type::Arrow) {
            ImVec2 start{};
            ImVec2 end{};
            if (!world_to_screen(annotation.world_pos.cast<float>(), vp_matrix, viewport_pos, viewport_size, start) ||
                !world_to_screen(annotation.world_pos_end.cast<float>(), vp_matrix, viewport_pos, viewport_size, end)) {
                continue;
            }
            draw_list->AddLine(start, end, color, annotation.thickness);
            draw_arrowhead(draw_list, end, start, color);
            if (!annotation.text.empty()) {
                const ImVec2 mid((start.x + end.x) * 0.5f, (start.y + end.y) * 0.5f);
                draw_text_box(draw_list, ImVec2(mid.x + 6.0f, mid.y - 18.0f), annotation.text, annotation.font_scale, annotation.color);
            }
            continue;
        }

        if (annotation.type == Annotation::Type::DimensionLine) {
            ImVec2 start{};
            ImVec2 end{};
            if (!world_to_screen(annotation.world_pos.cast<float>(), vp_matrix, viewport_pos, viewport_size, start) ||
                !world_to_screen(annotation.world_pos_end.cast<float>(), vp_matrix, viewport_pos, viewport_size, end)) {
                continue;
            }
            ImVec2 delta(end.x - start.x, end.y - start.y);
            const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
            if (length < 4.0f) {
                continue;
            }
            const ImVec2 dir(delta.x / length, delta.y / length);
            const ImVec2 perp(-dir.y, dir.x);
            const float offset = 12.0f;
            const float tick = 6.0f;
            const ImVec2 a(start.x + perp.x * offset, start.y + perp.y * offset);
            const ImVec2 b(end.x + perp.x * offset, end.y + perp.y * offset);
            draw_list->AddLine(a, b, color, annotation.thickness);
            draw_list->AddLine(ImVec2(a.x - perp.x * tick, a.y - perp.y * tick),
                               ImVec2(a.x + perp.x * tick, a.y + perp.y * tick), color, annotation.thickness);
            draw_list->AddLine(ImVec2(b.x - perp.x * tick, b.y - perp.y * tick),
                               ImVec2(b.x + perp.x * tick, b.y + perp.y * tick), color, annotation.thickness);
            const std::string label = annotation.auto_value ? format_distance(annotation) : annotation.text;
            if (!label.empty()) {
                const ImVec2 mid((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f);
                draw_text_box(draw_list,
                              ImVec2(mid.x - 18.0f, mid.y - 22.0f),
                              label,
                              annotation.font_scale,
                              annotation.color);
            }
            continue;
        }

        if (annotation.type == Annotation::Type::AngleArc) {
            ImVec2 vertex{};
            ImVec2 arm_a{};
            ImVec2 arm_b{};
            if (!world_to_screen(annotation.world_pos.cast<float>(), vp_matrix, viewport_pos, viewport_size, vertex) ||
                !world_to_screen(annotation.world_pos_end.cast<float>(), vp_matrix, viewport_pos, viewport_size, arm_a) ||
                !world_to_screen(annotation.world_pos_aux.cast<float>(), vp_matrix, viewport_pos, viewport_size, arm_b)) {
                continue;
            }

            const float radius = 20.0f;
            const Eigen::Vector2f va(arm_a.x - vertex.x, arm_a.y - vertex.y);
            const Eigen::Vector2f vb(arm_b.x - vertex.x, arm_b.y - vertex.y);
            if (va.norm() < 2.0f || vb.norm() < 2.0f) {
                continue;
            }
            float start_angle = std::atan2(va.y(), va.x());
            float end_angle = std::atan2(vb.y(), vb.x());
            float delta = end_angle - start_angle;
            while (delta <= -static_cast<float>(M_PI)) {
                delta += 2.0f * static_cast<float>(M_PI);
            }
            while (delta > static_cast<float>(M_PI)) {
                delta -= 2.0f * static_cast<float>(M_PI);
            }
            end_angle = start_angle + delta;
            const ImVec2 dir_a(va.x() / va.norm(), va.y() / va.norm());
            const ImVec2 dir_b(vb.x() / vb.norm(), vb.y() / vb.norm());
            draw_list->AddLine(vertex, ImVec2(vertex.x + dir_a.x * 28.0f, vertex.y + dir_a.y * 28.0f), color, annotation.thickness);
            draw_list->AddLine(vertex, ImVec2(vertex.x + dir_b.x * 28.0f, vertex.y + dir_b.y * 28.0f), color, annotation.thickness);
            ImVec2 prev(vertex.x + std::cos(start_angle) * radius, vertex.y + std::sin(start_angle) * radius);
            for (int i = 1; i <= 20; ++i) {
                const float t = static_cast<float>(i) / 20.0f;
                const float theta = start_angle + (end_angle - start_angle) * t;
                const ImVec2 cur(vertex.x + std::cos(theta) * radius, vertex.y + std::sin(theta) * radius);
                draw_list->AddLine(prev, cur, color, annotation.thickness);
                prev = cur;
            }
            const std::string label = annotation.auto_value
                                          ? format_angle(annotation.world_pos, annotation.world_pos_end, annotation.world_pos_aux)
                                          : annotation.text;
            if (label.empty()) {
                continue;
            }
            draw_text_box(draw_list,
                          ImVec2(vertex.x + std::cos((start_angle + end_angle) * 0.5f) * (radius + 10.0f),
                                 vertex.y + std::sin((start_angle + end_angle) * 0.5f) * (radius + 10.0f)),
                          label,
                          annotation.font_scale,
                          annotation.color);
            continue;
        }

        if (annotation.type == Annotation::Type::Circle) {
            ImVec2 center{};
            ImVec2 edge{};
            if (!world_to_screen(annotation.world_pos.cast<float>(), vp_matrix, viewport_pos, viewport_size, center) ||
                !world_to_screen(annotation.world_pos_end.cast<float>(), vp_matrix, viewport_pos, viewport_size, edge)) {
                continue;
            }
            const float dx = edge.x - center.x;
            const float dy = edge.y - center.y;
            const float radius = std::sqrt(dx * dx + dy * dy);
            if (radius < 1.0f) {
                continue;
            }
            draw_list->AddCircle(center, radius, color, 48, annotation.thickness);
            if (!annotation.text.empty()) {
                draw_text_box(draw_list, ImVec2(center.x + radius + 6.0f, center.y - 10.0f), annotation.text, annotation.font_scale, annotation.color);
            }
        }
    }
}

void AnnotationManager::set_default_style(const ImVec4& color, float font_scale, float thickness) {
    default_color_ = color;
    default_font_scale_ = std::max(font_scale, 0.1f);
    default_thickness_ = std::max(thickness, 0.5f);
}

ImVec4 AnnotationManager::default_color() const {
    return default_color_;
}

float AnnotationManager::default_font_scale() const {
    return default_font_scale_;
}

float AnnotationManager::default_thickness() const {
    return default_thickness_;
}

}  // namespace sbox::ui
