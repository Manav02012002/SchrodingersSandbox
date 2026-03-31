#version 410 core

in vec2 uv;
out float frag_ao;

uniform sampler2D u_position;
uniform sampler2D u_normal;
uniform sampler2D u_noise;
uniform mat4 u_projection;
uniform vec2 u_noise_scale;
uniform int u_kernel_size;
uniform vec3 u_samples[64];
uniform float u_radius;
uniform float u_bias;
uniform float u_power;

void main() {
    vec3 frag_pos = texture(u_position, uv).xyz;
    vec3 normal = normalize(texture(u_normal, uv).xyz);
    vec3 random_vec = texture(u_noise, uv * u_noise_scale).xyz;

    if (length(normal) < 0.01) {
        frag_ao = 1.0;
        return;
    }

    vec3 tangent = normalize(random_vec - normal * dot(random_vec, normal));
    vec3 bitangent = cross(normal, tangent);
    mat3 TBN = mat3(tangent, bitangent, normal);

    float occlusion = 0.0;
    for (int i = 0; i < u_kernel_size; i++) {
        vec3 sample_pos = frag_pos + TBN * u_samples[i] * u_radius;

        vec4 offset = u_projection * vec4(sample_pos, 1.0);
        offset.xyz /= offset.w;
        offset.xyz = offset.xyz * 0.5 + 0.5;

        float sample_depth = texture(u_position, offset.xy).z;
        float range_check = smoothstep(0.0, 1.0, u_radius / abs(frag_pos.z - sample_depth));
        occlusion += (sample_depth >= sample_pos.z + u_bias ? 1.0 : 0.0) * range_check;
    }

    frag_ao = 1.0 - pow(occlusion / float(u_kernel_size), u_power);
}
