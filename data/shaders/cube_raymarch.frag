#version 410 core

in vec2 uv;
in vec2 ndc;
out vec4 frag_color;

uniform mat4 u_inv_vp;
uniform vec3 u_camera_pos;
uniform int u_render_mode;
uniform float u_iso_value;
uniform float u_gamma;
uniform float u_max_density;
uniform float u_bound_radius;

uniform sampler3D u_volume;
uniform vec3 u_grid_origin;
uniform mat3 u_world_to_grid;
uniform ivec3 u_grid_dims;

const vec3 kBackground = vec3(0.04, 0.055, 0.09);

float sample_volume(vec3 world_pos) {
    vec3 grid_pos = u_world_to_grid * (world_pos - u_grid_origin);
    vec3 uvw = grid_pos / vec3(u_grid_dims);
    if (any(lessThan(uvw, vec3(0.0))) || any(greaterThan(uvw, vec3(1.0)))) {
        return 0.0;
    }
    return texture(u_volume, uvw).r;
}

vec3 density_ramp(float t) {
    if (t < 0.33) {
        return mix(vec3(0.04, 0.055, 0.09), vec3(0.12, 0.25, 0.48), clamp(t / 0.33, 0.0, 1.0));
    }
    if (t < 0.66) {
        return mix(vec3(0.12, 0.25, 0.48), vec3(0.27, 0.76, 0.80), clamp((t - 0.33) / 0.33, 0.0, 1.0));
    }
    return mix(vec3(0.27, 0.76, 0.80), vec3(1.0, 1.0, 1.0), clamp((t - 0.66) / 0.34, 0.0, 1.0));
}

float density(vec3 pos) {
    float v = sample_volume(pos);
    return v * v;
}

vec3 isosurface_shading(vec3 base_color, vec3 p, vec3 view_dir) {
    float eps = 0.01;
    vec3 dx = vec3(eps, 0.0, 0.0);
    vec3 dy = vec3(0.0, eps, 0.0);
    vec3 dz = vec3(0.0, 0.0, eps);

    float nx = density(p + dx) - density(p - dx);
    float ny = density(p + dy) - density(p - dy);
    float nz = density(p + dz) - density(p - dz);
    vec3 n = normalize(vec3(nx, ny, nz));

    vec3 light_dir = normalize(vec3(1.0, 1.0, 0.5));
    float ambient = 0.15;
    float diffuse = max(dot(n, light_dir), 0.0) * 0.7;
    float specular = pow(max(dot(reflect(-light_dir, n), view_dir), 0.0), 32.0) * 0.4;
    float rim = pow(1.0 - max(dot(n, view_dir), 0.0), 3.0) * 0.3;

    return base_color * (ambient + diffuse) + vec3(1.0) * (specular + rim);
}

bool find_isosurface(vec3 ro, vec3 rd, float t_start, float t_end, out vec3 hit_point) {
    const int steps = 192;
    float step_size = (t_end - t_start) / float(steps);
    float t_prev = t_start;
    float prev_psi = sample_volume(ro + rd * t_prev);
    float prev_density = prev_psi * prev_psi;

    for (int i = 1; i <= steps; ++i) {
        float t_curr = t_start + float(i) * step_size;
        float curr_psi = sample_volume(ro + rd * t_curr);
        float curr_density = curr_psi * curr_psi;

        if (prev_density < u_iso_value && curr_density >= u_iso_value) {
            float t_lo = t_prev;
            float t_hi = t_curr;
            for (int b = 0; b < 5; ++b) {
                float t_mid = 0.5 * (t_lo + t_hi);
                vec3 p_mid = ro + t_mid * rd;
                float d_mid = sample_volume(p_mid);
                if (d_mid * d_mid >= u_iso_value) {
                    t_hi = t_mid;
                } else {
                    t_lo = t_mid;
                }
            }
            hit_point = ro + t_hi * rd;
            return true;
        }
        t_prev = t_curr;
        prev_density = curr_density;
    }

    return false;
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

    if (u_render_mode == 0) {
        const int NUM_STEPS = 128;
        float step_size = (t_far - t_near) / float(NUM_STEPS);
        vec3 accum_color = vec3(0.0);
        float accum_alpha = 0.0;

        for (int i = 0; i < NUM_STEPS; ++i) {
            float t = t_near + (float(i) + 0.5) * step_size;
            vec3 pos = ray_origin + t * ray_dir;
            float val = sample_volume(pos);
            float dens = val * val;
            float normalized = clamp(dens / u_max_density, 0.0, 1.0);
            float mapped = pow(normalized, u_gamma);

            vec3 col = density_ramp(mapped);
            float sample_alpha = mapped * (3.0 / float(NUM_STEPS));
            accum_color += (1.0 - accum_alpha) * col * sample_alpha;
            accum_alpha += (1.0 - accum_alpha) * sample_alpha;

            if (accum_alpha > 0.95) {
                break;
            }
        }

        frag_color = vec4(accum_color + (1.0 - accum_alpha) * kBackground, 1.0);
        return;
    }

    vec3 hit;
    if (!find_isosurface(ray_origin, ray_dir, t_near, t_far, hit)) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

    vec3 view_dir = normalize(u_camera_pos - hit);

    if (u_render_mode == 1) {
        vec3 shaded = isosurface_shading(vec3(0.38, 0.65, 0.98), hit, view_dir);
        frag_color = vec4(shaded, 1.0);
        return;
    }

    float phase = sample_volume(hit);
    vec3 base = phase >= 0.0 ? vec3(0.23, 0.51, 0.96) : vec3(0.94, 0.27, 0.27);
    vec3 shaded = isosurface_shading(base, hit, view_dir);
    frag_color = vec4(shaded, 1.0);
}
