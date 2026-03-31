#version 410 core

in vec3 v_world_pos;
in vec3 v_pos_a;
in vec3 v_pos_b;
in vec3 v_color_a;
in vec3 v_color_b;
in float v_radius;
in float v_axis_t;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;

layout(location = 0) out vec3 g_position;
layout(location = 1) out vec3 g_normal;
layout(location = 2) out vec4 g_albedo;

void main() {
    vec3 axis_vec = v_pos_b - v_pos_a;
    float axis_len = length(axis_vec);
    if (axis_len < 1e-6) {
        discard;
    }

    vec3 axis = axis_vec / axis_len;
    vec3 cyl_a = v_pos_a - axis * v_radius;
    vec3 cyl_b = v_pos_b + axis * v_radius;
    float cyl_len = length(cyl_b - cyl_a);

    vec3 ray_origin = u_camera_pos;
    vec3 ray_dir = normalize(v_world_pos - u_camera_pos);

    vec3 delta_p = ray_origin - cyl_a;
    vec3 d_perp = ray_dir - axis * dot(ray_dir, axis);
    vec3 p_perp = delta_p - axis * dot(delta_p, axis);

    float a = dot(d_perp, d_perp);
    float b = 2.0 * dot(d_perp, p_perp);
    float c = dot(p_perp, p_perp) - v_radius * v_radius;

    if (abs(a) < 1e-8) {
        discard;
    }

    float disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        discard;
    }

    float sqrt_disc = sqrt(disc);
    float t0 = (-b - sqrt_disc) / (2.0 * a);
    float t1 = (-b + sqrt_disc) / (2.0 * a);
    float t = t0;
    if (t <= 0.0) {
        t = t1;
    }
    if (t <= 0.0) {
        discard;
    }

    vec3 hit = ray_origin + t * ray_dir;
    float axis_proj = dot(hit - cyl_a, axis);
    if (axis_proj < 0.0 || axis_proj > cyl_len) {
        discard;
    }

    vec3 closest = cyl_a + axis * axis_proj;
    vec3 normal = normalize(hit - closest);
    float bond_t = clamp(axis_proj / cyl_len, 0.0, 1.0);
    vec3 base_color = mix(v_color_a, v_color_b, step(0.5, bond_t));

    vec4 clip = u_proj * u_view * vec4(hit, 1.0);
    float ndc_depth = clip.z / clip.w;
    gl_FragDepth = ndc_depth * 0.5 + 0.5;

    g_position = (u_view * vec4(hit, 1.0)).xyz;
    g_normal = normalize(mat3(u_view) * normal);
    g_albedo = vec4(base_color, 0.4);
}
