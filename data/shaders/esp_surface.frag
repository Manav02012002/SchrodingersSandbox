#version 410 core

in vec2 uv;
in vec2 ndc;
out vec4 frag_color;

uniform mat4 u_inv_vp;
uniform vec3 u_camera_pos;
uniform float u_iso_value;
uniform float u_bound_radius;
uniform sampler3D u_density;
uniform sampler3D u_esp;
uniform vec3 u_grid_origin;
uniform mat3 u_world_to_grid;
uniform ivec3 u_grid_dims;
uniform float u_esp_min;
uniform float u_esp_max;

const vec3 kBackground = vec3(0.04, 0.055, 0.09);

float sample_grid(sampler3D tex, vec3 world_pos) {
    vec3 grid_pos = u_world_to_grid * (world_pos - u_grid_origin);
    vec3 uvw = grid_pos / vec3(u_grid_dims);
    if (any(lessThan(uvw, vec3(0.0))) || any(greaterThan(uvw, vec3(1.0)))) {
        return 0.0;
    }
    return texture(tex, uvw).r;
}

float sample_density(vec3 world_pos) {
    return sample_grid(u_density, world_pos);
}

float sample_esp(vec3 world_pos) {
    return sample_grid(u_esp, world_pos);
}

vec3 esp_colormap(float esp_value) {
    float t = clamp((esp_value - u_esp_min) / max(u_esp_max - u_esp_min, 1e-6), 0.0, 1.0);
    t = t * 2.0 - 1.0;

    if (t < -0.5) {
        return mix(vec3(0.0, 0.0, 0.6), vec3(0.2, 0.4, 1.0), (t + 1.0) * 2.0);
    } else if (t < 0.0) {
        return mix(vec3(0.2, 0.4, 1.0), vec3(0.9, 0.9, 0.95), (t + 0.5) * 2.0);
    } else if (t < 0.5) {
        return mix(vec3(0.9, 0.9, 0.95), vec3(1.0, 0.4, 0.2), t * 2.0);
    }
    return mix(vec3(1.0, 0.4, 0.2), vec3(0.6, 0.0, 0.0), (t - 0.5) * 2.0);
}

vec3 density_normal(vec3 p) {
    float eps = 0.02;
    vec3 dx = vec3(eps, 0.0, 0.0);
    vec3 dy = vec3(0.0, eps, 0.0);
    vec3 dz = vec3(0.0, 0.0, eps);
    float nx = sample_density(p + dx) - sample_density(p - dx);
    float ny = sample_density(p + dy) - sample_density(p - dy);
    float nz = sample_density(p + dz) - sample_density(p - dz);
    return normalize(vec3(nx, ny, nz));
}

bool find_isosurface(vec3 ro, vec3 rd, float t_start, float t_end, out vec3 hit_point) {
    const int steps = 256;
    float step_size = (t_end - t_start) / float(steps);
    float t_prev = t_start;
    float prev_density = sample_density(ro + rd * t_prev);

    for (int i = 1; i <= steps; ++i) {
        float t_curr = t_start + float(i) * step_size;
        float curr_density = sample_density(ro + rd * t_curr);
        if (prev_density < u_iso_value && curr_density >= u_iso_value) {
            float t_lo = t_prev;
            float t_hi = t_curr;
            for (int b = 0; b < 5; ++b) {
                float t_mid = 0.5 * (t_lo + t_hi);
                float d_mid = sample_density(ro + rd * t_mid);
                if (d_mid >= u_iso_value) {
                    t_hi = t_mid;
                } else {
                    t_lo = t_mid;
                }
            }
            hit_point = ro + rd * t_hi;
            return true;
        }
        t_prev = t_curr;
        prev_density = curr_density;
    }
    return false;
}

vec3 phong_shade(vec3 base_color, vec3 p, vec3 view_dir) {
    vec3 n = density_normal(p);
    vec3 light_dir = normalize(vec3(1.0, 1.0, 0.5));
    float ambient = 0.18;
    float diffuse = max(dot(n, light_dir), 0.0) * 0.70;
    float specular = pow(max(dot(reflect(-light_dir, n), view_dir), 0.0), 32.0) * 0.35;
    float rim = pow(1.0 - max(dot(n, view_dir), 0.0), 3.0) * 0.25;
    return base_color * (ambient + diffuse) + vec3(1.0) * (specular + rim);
}

void main() {
    vec2 ray_ndc = 2.0 * uv - vec2(1.0);
    vec4 near_clip = vec4(ray_ndc, -1.0, 1.0);
    vec4 far_clip = vec4(ray_ndc, 1.0, 1.0);

    vec4 near_world = u_inv_vp * near_clip;
    vec4 far_world = u_inv_vp * far_clip;
    near_world /= near_world.w;
    far_world /= far_world.w;

    vec3 ray_origin = near_world.xyz;
    vec3 ray_dir = normalize(far_world.xyz - near_world.xyz);

    vec3 oc = ray_origin;
    float b = dot(oc, ray_dir);
    float c = dot(oc, oc) - u_bound_radius * u_bound_radius;
    float disc = b * b - c;
    if (disc < 0.0) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

    float sqrt_disc = sqrt(disc);
    float t_near = -b - sqrt_disc;
    float t_far = -b + sqrt_disc;
    if (t_far < 0.0) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }
    t_near = max(t_near, 0.0);
    if (t_far <= t_near) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

    vec3 hit;
    if (!find_isosurface(ray_origin, ray_dir, t_near, t_far, hit)) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

    vec3 view_dir = normalize(u_camera_pos - hit);
    float esp = sample_esp(hit);
    vec3 base = esp_colormap(esp);
    vec3 shaded = phong_shade(base, hit, view_dir);
    frag_color = vec4(shaded, 1.0);
}
