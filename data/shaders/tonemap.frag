#version 410 core

in vec2 uv;
out vec4 frag_color;

uniform sampler2D u_input;
uniform float u_exposure;
uniform float u_gamma;

void main() {
    vec3 hdr = texture(u_input, uv).rgb;
    vec3 mapped = vec3(1.0) - exp(-hdr * u_exposure);
    mapped = pow(mapped, vec3(1.0 / u_gamma));
    frag_color = vec4(mapped, 1.0);
}
