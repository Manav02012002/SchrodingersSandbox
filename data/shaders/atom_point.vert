#version 410 core

layout(location = 0) in vec3 a_position;
layout(location = 1) in vec4 a_color_and_size;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;

out vec4 v_color;

void main() {
    vec4 view_pos = u_view * vec4(a_position, 1.0);
    float distance_to_camera = max(length(a_position - u_camera_pos), 1.0);
    gl_Position = u_proj * view_pos;
    gl_PointSize = max(2.0, a_color_and_size.w * 180.0 / distance_to_camera);
    v_color = vec4(a_color_and_size.rgb, 0.85);
}
