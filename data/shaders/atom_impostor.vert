#version 330 core

layout(location = 0) in vec2 a_quad;
layout(location = 1) in vec4 a_pos_radius;
layout(location = 2) in vec4 a_color_z;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;

out vec3 v_center;
out float v_radius;
out vec3 v_color;
out vec3 v_world_pos;

void main() {
    vec3 right = normalize(vec3(u_view[0][0], u_view[1][0], u_view[2][0]));
    vec3 up = normalize(vec3(u_view[0][1], u_view[1][1], u_view[2][1]));

    v_center = a_pos_radius.xyz;
    v_radius = a_pos_radius.w;
    v_color = a_color_z.xyz;
    v_world_pos = v_center + right * (a_quad.x * v_radius) + up * (a_quad.y * v_radius);

    gl_Position = u_proj * u_view * vec4(v_world_pos, 1.0);
}
