#version 410 core

in vec2 uv;
in vec2 ndc;
out vec4 frag_color;

uniform int u_n;
uniform int u_l;
uniform int u_m;
uniform float u_Zeff;
uniform mat4 u_inv_vp;
uniform vec3 u_camera_pos;
uniform float u_iso_value;
uniform int u_render_mode;
uniform float u_gamma;
uniform vec2 u_resolution;
uniform float u_max_density;

const float PI = 3.14159265358979323846;
const vec3 kBackground = vec3(0.04, 0.055, 0.09);

float glsl_factorial(int n) {
    const float table[21] = float[21](
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0,
        39916800.0,
        479001600.0,
        6227020800.0,
        87178291200.0,
        1307674368000.0,
        20922789888000.0,
        355687428096000.0,
        6402373705728000.0,
        121645100408832000.0,
        2432902008176640000.0
    );
    return table[clamp(n, 0, 20)];
}

float log_factorial(int n) {
    const float table[21] = float[21](
        0.0,
        0.0,
        0.69314718056,
        1.79175946923,
        3.17805383035,
        4.78749174278,
        6.57925121201,
        8.52516136107,
        10.60460290275,
        12.80182748008,
        15.10441257308,
        17.50230784587,
        19.98721449566,
        22.55216385312,
        25.19122118274,
        27.89927138384,
        30.67186010608,
        33.50507345014,
        36.39544520803,
        39.33988418720,
        42.33561646075
    );
    return table[clamp(n, 0, 20)];
}

float assoc_laguerre(int n, float alpha, float x) {
    if (n < 0) {
        return 0.0;
    }
    if (n == 0) {
        return 1.0;
    }

    float l_nm2 = 1.0;
    float l_nm1 = 1.0 + alpha - x;
    if (n == 1) {
        return l_nm1;
    }

    for (int k = 2; k <= 32; ++k) {
        if (k > n) {
            break;
        }
        float kd = float(k);
        float term1 = (2.0 * kd - 1.0 + alpha - x) * l_nm1;
        float term2 = (kd - 1.0 + alpha) * l_nm2;
        float l_n = (term1 - term2) / kd;
        l_nm2 = l_nm1;
        l_nm1 = l_n;
    }

    return l_nm1;
}

float assoc_legendre(int l, int m, float x) {
    if (l < 0) {
        return 0.0;
    }

    int abs_m = abs(m);
    if (abs_m > l) {
        return 0.0;
    }

    x = clamp(x, -1.0, 1.0);

    float p_mm = 1.0;
    if (abs_m > 0) {
        float one_minus_x2 = max(0.0, 1.0 - x * x);
        float root = sqrt(one_minus_x2);
        float fact = 1.0;
        for (int i = 1; i <= 16; ++i) {
            if (i > abs_m) {
                break;
            }
            p_mm *= (-fact) * root;
            fact += 2.0;
        }
    }

    if (l == abs_m) {
        return p_mm;
    }

    float p_m1m = x * (2.0 * float(abs_m) + 1.0) * p_mm;
    if (l == abs_m + 1) {
        return p_m1m;
    }

    float p_lm2 = p_mm;
    float p_lm1 = p_m1m;
    for (int ell = 2; ell <= 24; ++ell) {
        int actual_l = abs_m + ell;
        if (actual_l > l) {
            break;
        }
        float ell_f = float(actual_l);
        float numerator = (2.0 * ell_f - 1.0) * x * p_lm1 - float(actual_l + abs_m - 1) * p_lm2;
        float p_lm = numerator / float(actual_l - abs_m);
        p_lm2 = p_lm1;
        p_lm1 = p_lm;
    }

    return p_lm1;
}

float real_spherical_harmonic(int l, int m, float theta, float phi) {
    if (l < 0 || abs(m) > l) {
        return 0.0;
    }

    int abs_m = abs(m);
    float x = cos(theta);
    float p_lm = assoc_legendre(l, abs_m, x);

    float norm = sqrt(((2.0 * float(l) + 1.0) / (4.0 * PI)) *
                      (glsl_factorial(l - abs_m) / glsl_factorial(l + abs_m)));

    if (m > 0) {
        return sqrt(2.0) * norm * p_lm * cos(float(m) * phi);
    }
    if (m < 0) {
        return sqrt(2.0) * norm * p_lm * sin(float(abs_m) * phi);
    }
    return norm * p_lm;
}

float radial_wavefunction(int n, int l, float Zeff, float r) {
    if (n <= 0 || l < 0 || l >= n || Zeff <= 0.0 || r < 0.0) {
        return 0.0;
    }

    float nf = float(n);
    float radial_scale = 2.0 * Zeff / nf;
    float rho = radial_scale * r;

    if (rho > 80.0) {
        return 0.0;
    }

    int idx1 = clamp(n - l - 1, 0, 20);
    int idx2 = clamp(n + l, 0, 20);

    float log_N = 0.5 * (3.0 * log(2.0 * Zeff / nf)
                       + log_factorial(idx1)
                       - log(2.0 * nf)
                       - 3.0 * log_factorial(idx2));
    float normalization = exp(log_N);

    float laguerre = assoc_laguerre(n - l - 1, float(2 * l + 1), rho);
    return normalization * exp(-0.5 * rho) * pow(max(rho, 0.0), float(l)) * laguerre;
}

float psi(vec3 pos) {
    float r = length(pos);
    if (r < 1e-10) {
        return 0.0;
    }
    float theta = acos(clamp(pos.z / r, -1.0, 1.0));
    float phi = atan(pos.y, pos.x);

    float radial = radial_wavefunction(u_n, u_l, u_Zeff, r);
    float angular = real_spherical_harmonic(u_l, u_m, theta, phi);
    return radial * angular;
}

bool ray_sphere_intersect(vec3 ro, vec3 rd, float radius, out float t0, out float t1) {
    float b = dot(ro, rd);
    float c = dot(ro, ro) - radius * radius;
    float h = b * b - c;
    if (h < 0.0) {
        return false;
    }

    float s = sqrt(h);
    t0 = -b - s;
    t1 = -b + s;
    if (t1 < 0.0) {
        return false;
    }
    t0 = max(t0, 0.0);
    return true;
}

vec3 density_ramp(float t) {
    if (t < 0.33) {
        return mix(vec3(0.04, 0.055, 0.09), vec3(0.12, 0.25, 0.48), t / 0.33);
    }
    if (t < 0.66) {
        return mix(vec3(0.12, 0.25, 0.48), vec3(0.27, 0.76, 0.80), (t - 0.33) / 0.33);
    }
    return mix(vec3(0.27, 0.76, 0.80), vec3(1.0, 1.0, 1.0), (t - 0.66) / 0.34);
}

float density(vec3 pos) {
    float v = psi(pos);
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
    const int steps = 256;
    float step_size = (t_end - t_start) / float(steps);
    float t_prev = t_start;
    float prev_psi = psi(ro + rd * t_prev);
    float prev_density = prev_psi * prev_psi;

    for (int i = 1; i <= steps; ++i) {
        float t_curr = t_start + float(i) * step_size;
        float curr_psi = psi(ro + rd * t_curr);
        float curr_density = curr_psi * curr_psi;

        if (prev_density < u_iso_value && curr_density >= u_iso_value) {
            float t_lo = t_prev;
            float t_hi = t_curr;
            for (int b = 0; b < 5; ++b) {
                float t_mid = 0.5 * (t_lo + t_hi);
                vec3 p_mid = ro + t_mid * rd;
                float d_mid = psi(p_mid);
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
    vec3 ray_dir = far_world.xyz - near_world.xyz;
    ray_dir = normalize(ray_dir);

    if (u_n == 0) {
        frag_color = vec4(1.0, 0.0, 0.0, 1.0);
        return;
    }

    vec3 ro = ray_origin;
    vec3 rd = ray_dir;
    vec3 bg = kBackground;

    float bound_r = max(5.0 * float(u_n * u_n) / u_Zeff, 3.0);
    vec3 oc = ray_origin - vec3(0.0);  // sphere at origin
    float b = dot(oc, ray_dir);
    float c = dot(oc, oc) - bound_r * bound_r;
    float disc = b * b - c;
    if (disc < 0.0) { frag_color = vec4(bg, 1.0); return; }
    float sqrt_disc = sqrt(disc);
    float t_near = -b - sqrt_disc;
    float t_far = -b + sqrt_disc;
    if (t_far < 0.0) { frag_color = vec4(bg, 1.0); return; }
    t_near = max(t_near, 0.0);

    float t0 = t_near;
    float t1 = t_far;

    if (t1 <= t0) {
        frag_color = vec4(bg, 1.0);
        return;
    }

    if (u_render_mode == 0) {
        const int NUM_STEPS = 192;

        // Volume rendering
        float step_size = (t_far - t_near) / float(NUM_STEPS);
        vec3 accum_color = vec3(0.0);
        float accum_alpha = 0.0;

        // First pass: find max density along this ray for normalisation
        float ray_max = 0.0;
        for (int i = 0; i < NUM_STEPS; i++) {
            float t = t_near + (float(i) + 0.5) * step_size;
            vec3 pos = ray_origin + t * ray_dir;
            float val = psi(pos);
            float dens = val * val;
            ray_max = max(ray_max, dens);
        }
        if (ray_max < 1e-20) {
            frag_color = vec4(bg, 1.0);
            return;
        }

        // Second pass: accumulate colour
        for (int i = 0; i < NUM_STEPS; i++) {
            float t = t_near + (float(i) + 0.5) * step_size;
            vec3 pos = ray_origin + t * ray_dir;
            float val = psi(pos);
            float dens = val * val;
            float normalized = dens / ray_max;
            float mapped = pow(normalized, u_gamma);

            // Colour ramp: dark blue -> cyan -> white
            vec3 col;
            if (mapped < 0.33) {
                col = mix(vec3(0.04, 0.055, 0.09), vec3(0.12, 0.25, 0.48), mapped / 0.33);
            } else if (mapped < 0.66) {
                col = mix(vec3(0.12, 0.25, 0.48), vec3(0.18, 0.76, 0.80), (mapped - 0.33) / 0.33);
            } else {
                col = mix(vec3(0.18, 0.76, 0.80), vec3(1.0, 1.0, 1.0), (mapped - 0.66) / 0.34);
            }

            float sample_alpha = mapped * (3.0 / float(NUM_STEPS));
            accum_color += (1.0 - accum_alpha) * col * sample_alpha;
            accum_alpha += (1.0 - accum_alpha) * sample_alpha;

            if (accum_alpha > 0.95) break;
        }

        frag_color = vec4(accum_color + (1.0 - accum_alpha) * bg, 1.0);
        return;
    }

    vec3 hit;
    if (!find_isosurface(ro, rd, t0, t1, hit)) {
        frag_color = vec4(kBackground, 1.0);
        return;
    }

    vec3 view_dir = normalize(u_camera_pos - hit);

    if (u_render_mode == 1) {
        vec3 shaded = isosurface_shading(vec3(0.38, 0.65, 0.98), hit, view_dir);
        frag_color = vec4(shaded, 1.0);
        return;
    }

    float phase = psi(hit);
    vec3 base = phase >= 0.0 ? vec3(0.23, 0.51, 0.96) : vec3(0.94, 0.27, 0.27);
    vec3 shaded = isosurface_shading(base, hit, view_dir);
    frag_color = vec4(shaded, 1.0);
}
