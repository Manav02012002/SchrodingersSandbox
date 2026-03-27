#version 330 core

layout(location = 0) in vec2 a_quad;
layout(location = 1) in vec3 a_pos_a;
layout(location = 2) in vec3 a_pos_b;
layout(location = 3) in vec3 a_color_a;
layout(location = 4) in vec3 a_color_b;
layout(location = 5) in float a_radius;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;
uniform int u_wireframe;

out vec3 v_world_pos;
out vec3 v_pos_a;
out vec3 v_pos_b;
out vec3 v_color_a;
out vec3 v_color_b;
out float v_radius;
out float v_axis_t;

void main() {
    vec3 axis = normalize(a_pos_b - a_pos_a);
    vec3 midpoint = 0.5 * (a_pos_a + a_pos_b);
    float half_length = 0.5 * length(a_pos_b - a_pos_a) + a_radius;
    vec3 extended_a = a_pos_a - axis * a_radius;
    vec3 extended_b = a_pos_b + axis * a_radius;

    vec3 view_dir = normalize(u_camera_pos - midpoint);
    vec3 side = cross(view_dir, axis);
    if (length(side) < 1e-5) {
        side = vec3(u_view[0][0], u_view[1][0], u_view[2][0]);
    }
    side = normalize(side);

    float along = a_quad.x;
    float width = (u_wireframe != 0) ? 0.0 : a_quad.y * a_radius;
    v_world_pos = midpoint + axis * (along * half_length) + side * width;
    if (u_wireframe != 0) {
        float t = 0.5 * (a_quad.x + 1.0);
        v_world_pos = mix(extended_a, extended_b, t);
    }

    v_pos_a = a_pos_a;
    v_pos_b = a_pos_b;
    v_color_a = a_color_a;
    v_color_b = a_color_b;
    v_radius = a_radius;
    v_axis_t = clamp(0.5 * (along + 1.0), 0.0, 1.0);

    gl_Position = u_proj * u_view * vec4(v_world_pos, 1.0);
}
