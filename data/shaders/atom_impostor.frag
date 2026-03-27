#version 330 core

in vec3 v_center;
in float v_radius;
in vec3 v_color;
in vec3 v_world_pos;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform vec3 u_camera_pos;

out vec4 frag_color;

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

    vec3 light_dir = normalize(vec3(1.0, 1.0, 0.5));
    vec3 view_dir = normalize(u_camera_pos - hit);
    vec3 half_dir = normalize(light_dir + view_dir);

    float ambient = 0.15;
    float diffuse = max(dot(normal, light_dir), 0.0);
    float specular = pow(max(dot(normal, half_dir), 0.0), 64.0);
    float rim = pow(1.0 - max(dot(normal, view_dir), 0.0), 2.0);

    vec3 color = v_color * (ambient + 0.85 * diffuse) + vec3(0.35) * specular + v_color * 0.2 * rim;
    frag_color = vec4(color, 1.0);
}
