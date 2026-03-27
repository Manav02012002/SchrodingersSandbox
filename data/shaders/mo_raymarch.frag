#version 410 core

in vec2 uv;
in vec2 ndc;
out vec4 frag_color;

uniform mat4 u_inv_vp;
uniform vec3 u_camera_pos;
uniform vec2 u_resolution;
uniform int u_render_mode;
uniform float u_iso_value;
uniform float u_gamma;
uniform float u_max_density;
uniform int u_mo_index;
uniform int u_num_shells;
uniform int u_num_basis;
uniform int u_num_mo;
uniform float u_bound_radius;

uniform sampler2D u_shell_desc;
uniform sampler2D u_shell_meta;
uniform sampler2D u_primitives;
uniform sampler2D u_mo_coeffs;

const vec3 kBackground = vec3(0.04, 0.055, 0.09);

vec4 fetch_shell_desc(int i) { return texelFetch(u_shell_desc, ivec2(i, 0), 0); }
vec4 fetch_shell_meta(int i) { return texelFetch(u_shell_meta, ivec2(i, 0), 0); }
vec4 fetch_primitive(int i) { return texelFetch(u_primitives, ivec2(i, 0), 0); }

float fetch_mo_coeff(int basis_idx, int mo_idx) {
    int col = mo_idx / 4;
    int channel = mo_idx - col * 4;
    vec4 texel = texelFetch(u_mo_coeffs, ivec2(basis_idx, col), 0);
    if (channel == 0) return texel.r;
    if (channel == 1) return texel.g;
    if (channel == 2) return texel.b;
    return texel.a;
}

float evaluate_mo(vec3 pos) {
    float mo_val = 0.0;

    for (int s = 0; s < 512; ++s) {
        if (s >= u_num_shells) {
            break;
        }

        vec4 desc = fetch_shell_desc(s);
        vec3 center = desc.xyz;
        int L = int(desc.w + 0.5);

        vec4 meta = fetch_shell_meta(s);
        int num_prims = int(meta.x + 0.5);
        int first_prim = int(meta.y + 0.5);
        int num_basis_in_shell = int(meta.z + 0.5);
        int first_basis = int(meta.w + 0.5);

        vec3 dr = pos - center;
        float r2 = dot(dr, dr);

        float radial = 0.0;
        for (int p = 0; p < 32; ++p) {
            if (p >= num_prims) {
                break;
            }
            vec4 prim = fetch_primitive(first_prim + p);
            radial += prim.y * exp(-prim.x * r2);
        }

        if (L == 0) {
            float chi = radial;
            mo_val += fetch_mo_coeff(first_basis, u_mo_index) * chi;
        } else if (L == 1) {
            mo_val += fetch_mo_coeff(first_basis + 0, u_mo_index) * (dr.x * radial);
            mo_val += fetch_mo_coeff(first_basis + 1, u_mo_index) * (dr.y * radial);
            mo_val += fetch_mo_coeff(first_basis + 2, u_mo_index) * (dr.z * radial);
        } else if (L == 2) {
            float xx = dr.x * dr.x * radial;
            float yy = dr.y * dr.y * radial;
            float zz = dr.z * dr.z * radial;
            float xy = dr.x * dr.y * radial;
            float xz = dr.x * dr.z * radial;
            float yz = dr.y * dr.z * radial;

            if (num_basis_in_shell == 5) {
                mo_val += fetch_mo_coeff(first_basis + 0, u_mo_index) * (0.5 * (2.0 * zz - xx - yy));
                mo_val += fetch_mo_coeff(first_basis + 1, u_mo_index) * (sqrt(3.0) * xz);
                mo_val += fetch_mo_coeff(first_basis + 2, u_mo_index) * (sqrt(3.0) * yz);
                mo_val += fetch_mo_coeff(first_basis + 3, u_mo_index) * (0.5 * sqrt(3.0) * (xx - yy));
                mo_val += fetch_mo_coeff(first_basis + 4, u_mo_index) * (sqrt(3.0) * xy);
            } else {
                mo_val += fetch_mo_coeff(first_basis + 0, u_mo_index) * xx;
                mo_val += fetch_mo_coeff(first_basis + 1, u_mo_index) * yy;
                mo_val += fetch_mo_coeff(first_basis + 2, u_mo_index) * zz;
                mo_val += fetch_mo_coeff(first_basis + 3, u_mo_index) * xy;
                mo_val += fetch_mo_coeff(first_basis + 4, u_mo_index) * xz;
                mo_val += fetch_mo_coeff(first_basis + 5, u_mo_index) * yz;
            }
        } else if (L == 3) {
            float xx = dr.x * dr.x;
            float yy = dr.y * dr.y;
            float zz = dr.z * dr.z;

            float xxx = dr.x * xx * radial;
            float yyy = dr.y * yy * radial;
            float zzz = dr.z * zz * radial;
            float xyy = dr.x * yy * radial;
            float xxy = xx * dr.y * radial;
            float xxz = xx * dr.z * radial;
            float xzz = dr.x * zz * radial;
            float yzz = dr.y * zz * radial;
            float yyz = yy * dr.z * radial;
            float xyz = dr.x * dr.y * dr.z * radial;

            if (num_basis_in_shell == 7) {
                mo_val += fetch_mo_coeff(first_basis + 0, u_mo_index) * (0.5 * dr.z * (2.0 * zz - 3.0 * xx - 3.0 * yy) * radial);
                mo_val += fetch_mo_coeff(first_basis + 1, u_mo_index) * (sqrt(3.0 / 8.0) * dr.x * (4.0 * zz - xx - yy) * radial);
                mo_val += fetch_mo_coeff(first_basis + 2, u_mo_index) * (sqrt(3.0 / 8.0) * dr.y * (4.0 * zz - xx - yy) * radial);
                mo_val += fetch_mo_coeff(first_basis + 3, u_mo_index) * (0.5 * sqrt(15.0) * dr.z * (xx - yy) * radial);
                mo_val += fetch_mo_coeff(first_basis + 4, u_mo_index) * (sqrt(15.0) * xyz);
                mo_val += fetch_mo_coeff(first_basis + 5, u_mo_index) * (sqrt(5.0 / 8.0) * dr.x * (xx - 3.0 * yy) * radial);
                mo_val += fetch_mo_coeff(first_basis + 6, u_mo_index) * (sqrt(5.0 / 8.0) * dr.y * (3.0 * xx - yy) * radial);
            } else {
                mo_val += fetch_mo_coeff(first_basis + 0, u_mo_index) * xxx;
                mo_val += fetch_mo_coeff(first_basis + 1, u_mo_index) * yyy;
                mo_val += fetch_mo_coeff(first_basis + 2, u_mo_index) * zzz;
                mo_val += fetch_mo_coeff(first_basis + 3, u_mo_index) * xyy;
                mo_val += fetch_mo_coeff(first_basis + 4, u_mo_index) * xxy;
                mo_val += fetch_mo_coeff(first_basis + 5, u_mo_index) * xxz;
                mo_val += fetch_mo_coeff(first_basis + 6, u_mo_index) * xzz;
                mo_val += fetch_mo_coeff(first_basis + 7, u_mo_index) * yzz;
                mo_val += fetch_mo_coeff(first_basis + 8, u_mo_index) * yyz;
                mo_val += fetch_mo_coeff(first_basis + 9, u_mo_index) * xyz;
            }
        }
    }

    return mo_val;
}

float density_ramp_coord(float t) {
    return clamp(t, 0.0, 1.0);
}

vec3 density_ramp(float t) {
    if (t < 0.33) {
        return mix(vec3(0.04, 0.055, 0.09), vec3(0.12, 0.25, 0.48), density_ramp_coord(t / 0.33));
    }
    if (t < 0.66) {
        return mix(vec3(0.12, 0.25, 0.48), vec3(0.27, 0.76, 0.80), density_ramp_coord((t - 0.33) / 0.33));
    }
    return mix(vec3(0.27, 0.76, 0.80), vec3(1.0, 1.0, 1.0), density_ramp_coord((t - 0.66) / 0.34));
}

float density(vec3 pos) {
    float v = evaluate_mo(pos);
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
    float prev_psi = evaluate_mo(ro + rd * t_prev);
    float prev_density = prev_psi * prev_psi;

    for (int i = 1; i <= steps; ++i) {
        float t_curr = t_start + float(i) * step_size;
        float curr_psi = evaluate_mo(ro + rd * t_curr);
        float curr_density = curr_psi * curr_psi;

        if (prev_density < u_iso_value && curr_density >= u_iso_value) {
            float t_lo = t_prev;
            float t_hi = t_curr;
            for (int b = 0; b < 5; ++b) {
                float t_mid = 0.5 * (t_lo + t_hi);
                vec3 p_mid = ro + t_mid * rd;
                float d_mid = evaluate_mo(p_mid);
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

    if (u_mo_index < 0 || u_mo_index >= u_num_mo || u_num_basis <= 0 || u_num_shells <= 0) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

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
            float val = evaluate_mo(pos);
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

    float phase = evaluate_mo(hit);
    vec3 base = phase >= 0.0 ? vec3(0.23, 0.51, 0.96) : vec3(0.94, 0.27, 0.27);
    vec3 shaded = isosurface_shading(base, hit, view_dir);
    frag_color = vec4(shaded, 1.0);
}
