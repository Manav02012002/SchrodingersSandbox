#version 410 core

layout(location = 0) in vec3 a_position;
layout(location = 1) in vec4 a_instance_pos_radius;

uniform mat4 u_light_vp;

void main() {
    vec3 world_pos = a_position * a_instance_pos_radius.w + a_instance_pos_radius.xyz;
    gl_Position = u_light_vp * vec4(world_pos, 1.0);
}
