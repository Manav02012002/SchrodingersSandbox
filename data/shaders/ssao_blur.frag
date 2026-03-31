#version 410 core

in vec2 uv;
out float frag_ao;

uniform sampler2D u_ssao_input;
uniform vec2 u_texel_size;

void main() {
    float result = 0.0;
    for (int x = -2; x < 2; x++) {
        for (int y = -2; y < 2; y++) {
            vec2 offset = vec2(float(x), float(y)) * u_texel_size;
            result += texture(u_ssao_input, uv + offset).r;
        }
    }
    frag_ao = result / 16.0;
}
