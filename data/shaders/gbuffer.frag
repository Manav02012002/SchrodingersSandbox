#version 410 core

in vec3 v_center;
in float v_radius;
in vec3 v_color;
in vec3 v_world_pos;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;

layout(location = 0) out vec3 g_position;
layout(location = 1) out vec3 g_normal;
layout(location = 2) out vec4 g_albedo;

void main() {
    vec3 ray_origin = u_camera_pos;
    vec3 ray_dir = normalize(v_world_pos - u_camera_pos);

    vec3 oc = ray_origin - v_center;
    float b = dot(oc, ray_dir);
    float c = dot(oc, oc) - v_radius * v_radius;
    float disc = b * b - c;
    if (disc < 0.0) {
        discard;
    }

    float sqrt_disc = sqrt(disc);
    float t = -b - sqrt_disc;
    if (t <= 0.0) {
        t = -b + sqrt_disc;
    }
    if (t <= 0.0) {
        discard;
    }

    vec3 hit = ray_origin + t * ray_dir;
    vec3 normal = normalize(hit - v_center);

    vec4 clip = u_proj * u_view * vec4(hit, 1.0);
    float ndc_depth = clip.z / clip.w;
    gl_FragDepth = ndc_depth * 0.5 + 0.5;

    g_position = (u_view * vec4(hit, 1.0)).xyz;
    g_normal = normalize(mat3(u_view) * normal);
    g_albedo = vec4(v_color, 0.4);
}
