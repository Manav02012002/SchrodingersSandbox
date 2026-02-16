#version 410 core

in vec2 uv;
out vec4 frag_color;

void main() {
    const vec3 edge_color = vec3(0.0392156863, 0.0549019608, 0.0901960784);
    const vec3 center_color = vec3(0.1764705882, 0.8313725490, 0.7490196078);

    float d = clamp(distance(uv, vec2(0.5, 0.5)) / 0.70710678, 0.0, 1.0);
    float t = smoothstep(1.0, 0.0, d);
    vec3 color = mix(edge_color, center_color, t);

    frag_color = vec4(color, 1.0);
}
