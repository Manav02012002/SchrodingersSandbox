#version 410 core

in vec4 v_color;

out vec4 frag_color;

void main() {
    vec2 uv = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(uv, uv);
    if (r2 > 1.0) {
        discard;
    }
    float alpha = smoothstep(1.0, 0.65, 1.0 - r2);
    float dim = 0.90 - 0.20 * sqrt(clamp(r2, 0.0, 1.0));
    frag_color = vec4(v_color.rgb * dim, v_color.a * alpha);
}
