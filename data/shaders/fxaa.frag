#version 410 core

in vec2 uv;
out vec4 frag_color;

uniform sampler2D u_input;
uniform vec2 u_texel_size;

#define FXAA_EDGE_THRESHOLD 0.0625
#define FXAA_EDGE_THRESHOLD_MIN 0.0312
#define FXAA_SEARCH_STEPS 12
#define FXAA_SUBPIX_QUALITY 0.75

float rgb_to_luma(vec3 rgb) {
    return dot(rgb, vec3(0.299, 0.587, 0.114));
}

void main() {
    vec3 color_center = texture(u_input, uv).rgb;
    float luma_center = rgb_to_luma(color_center);

    float luma_N = rgb_to_luma(texture(u_input, uv + vec2(0.0, u_texel_size.y)).rgb);
    float luma_S = rgb_to_luma(texture(u_input, uv - vec2(0.0, u_texel_size.y)).rgb);
    float luma_E = rgb_to_luma(texture(u_input, uv + vec2(u_texel_size.x, 0.0)).rgb);
    float luma_W = rgb_to_luma(texture(u_input, uv - vec2(u_texel_size.x, 0.0)).rgb);

    float luma_min = min(luma_center, min(min(luma_N, luma_S), min(luma_E, luma_W)));
    float luma_max = max(luma_center, max(max(luma_N, luma_S), max(luma_E, luma_W)));
    float luma_range = luma_max - luma_min;

    if (luma_range < max(FXAA_EDGE_THRESHOLD_MIN, luma_max * FXAA_EDGE_THRESHOLD)) {
        frag_color = vec4(color_center, 1.0);
        return;
    }

    float luma_NW = rgb_to_luma(texture(u_input, uv + vec2(-u_texel_size.x, u_texel_size.y)).rgb);
    float luma_NE = rgb_to_luma(texture(u_input, uv + u_texel_size).rgb);
    float luma_SW = rgb_to_luma(texture(u_input, uv - u_texel_size).rgb);
    float luma_SE = rgb_to_luma(texture(u_input, uv + vec2(u_texel_size.x, -u_texel_size.y)).rgb);

    float edge_h = abs(-2.0 * luma_W + luma_NW + luma_SW) +
                   abs(-2.0 * luma_center + luma_N + luma_S) * 2.0 +
                   abs(-2.0 * luma_E + luma_NE + luma_SE);
    float edge_v = abs(-2.0 * luma_N + luma_NW + luma_NE) +
                   abs(-2.0 * luma_center + luma_W + luma_E) * 2.0 +
                   abs(-2.0 * luma_S + luma_SW + luma_SE);

    bool is_horizontal = edge_h >= edge_v;

    float luma_neg = is_horizontal ? luma_S : luma_W;
    float luma_pos = is_horizontal ? luma_N : luma_E;
    float gradient_neg = abs(luma_neg - luma_center);
    float gradient_pos = abs(luma_pos - luma_center);

    float step_length = is_horizontal ? u_texel_size.y : u_texel_size.x;
    float luma_local_avg;

    if (gradient_neg >= gradient_pos) {
        step_length = -step_length;
        luma_local_avg = 0.5 * (luma_neg + luma_center);
    } else {
        luma_local_avg = 0.5 * (luma_pos + luma_center);
    }

    vec2 current_uv = uv;
    if (is_horizontal) {
        current_uv.y += step_length * 0.5;
    } else {
        current_uv.x += step_length * 0.5;
    }

    vec2 offset = is_horizontal ? vec2(u_texel_size.x, 0.0) : vec2(0.0, u_texel_size.y);
    vec2 uv_neg = current_uv - offset;
    vec2 uv_pos = current_uv + offset;

    float luma_end_neg = rgb_to_luma(texture(u_input, uv_neg).rgb) - luma_local_avg;
    float luma_end_pos = rgb_to_luma(texture(u_input, uv_pos).rgb) - luma_local_avg;

    bool reached_neg = abs(luma_end_neg) >= gradient_neg * 0.5;
    bool reached_pos = abs(luma_end_pos) >= gradient_pos * 0.5;

    for (int i = 1; i < FXAA_SEARCH_STEPS; i++) {
        if (!reached_neg) {
            uv_neg -= offset;
            luma_end_neg = rgb_to_luma(texture(u_input, uv_neg).rgb) - luma_local_avg;
            reached_neg = abs(luma_end_neg) >= gradient_neg * 0.5;
        }
        if (!reached_pos) {
            uv_pos += offset;
            luma_end_pos = rgb_to_luma(texture(u_input, uv_pos).rgb) - luma_local_avg;
            reached_pos = abs(luma_end_pos) >= gradient_pos * 0.5;
        }
        if (reached_neg && reached_pos) {
            break;
        }
    }

    float dist_neg = is_horizontal ? (uv.x - uv_neg.x) : (uv.y - uv_neg.y);
    float dist_pos = is_horizontal ? (uv_pos.x - uv.x) : (uv_pos.y - uv.y);
    float shorter_dist = min(dist_neg, dist_pos);
    float total_dist = dist_neg + dist_pos;

    bool is_direction_neg = dist_neg < dist_pos;
    float luma_end = is_direction_neg ? luma_end_neg : luma_end_pos;

    if (((luma_center - luma_local_avg) < 0.0) == ((luma_end - luma_local_avg) < 0.0)) {
        frag_color = vec4(color_center, 1.0);
        return;
    }

    float pixel_offset = -shorter_dist / total_dist + 0.5;

    float luma_avg = (luma_N + luma_S + luma_E + luma_W) * 0.25;
    float subpix = clamp(abs(luma_avg - luma_center) / luma_range, 0.0, 1.0);
    subpix = smoothstep(0.0, 1.0, subpix);
    subpix = subpix * subpix * FXAA_SUBPIX_QUALITY;

    pixel_offset = max(pixel_offset, subpix);

    vec2 final_uv = uv;
    if (is_horizontal) {
        final_uv.y += pixel_offset * step_length;
    } else {
        final_uv.x += pixel_offset * step_length;
    }

    frag_color = vec4(texture(u_input, final_uv).rgb, 1.0);
}
