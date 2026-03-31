#version 410 core

in vec2 uv;
out vec4 frag_color;

uniform sampler2D u_position;
uniform sampler2D u_normal;
uniform sampler2D u_albedo;
uniform sampler2D u_depth;
uniform vec3 u_camera_pos;
uniform vec3 u_light_dir_view;
uniform sampler2DShadow u_shadow_map;
uniform mat4 u_light_vp;
uniform mat4 u_inv_view;
uniform float u_ssao_enabled;
uniform sampler2D u_ssao_texture;
uniform float u_shadow_enabled;

void main() {
    vec3 position = texture(u_position, uv).rgb;
    vec3 normal = texture(u_normal, uv).rgb;
    vec4 albedo_spec = texture(u_albedo, uv);
    float depth = texture(u_depth, uv).r;

    if (depth >= 0.9999) {
        frag_color = vec4(0.04, 0.055, 0.09, 1.0);
        return;
    }

    vec3 albedo = albedo_spec.rgb;
    float specular_strength = albedo_spec.a;

    vec3 light_dir = normalize(u_light_dir_view);
    vec3 view_dir = normalize(-position);

    float ambient = 0.15;
    float diffuse = max(dot(normal, light_dir), 0.0) * 0.7;
    float spec = pow(max(dot(reflect(-light_dir, normal), view_dir), 0.0), 32.0) * specular_strength;
    float rim = pow(1.0 - max(dot(normal, view_dir), 0.0), 3.0) * 0.25;

    float ao = 1.0;
    if (u_ssao_enabled > 0.5) {
        ao = texture(u_ssao_texture, uv).r;
    }

    float shadow = 1.0;
    if (u_shadow_enabled > 0.5) {
        vec3 world_pos = (u_inv_view * vec4(position, 1.0)).xyz;
        vec4 light_space = u_light_vp * vec4(world_pos, 1.0);
        vec3 proj = light_space.xyz / light_space.w;
        proj = proj * 0.5 + 0.5;

        float bias = 0.002;
        shadow = 0.0;
        vec2 texel = 1.0 / vec2(textureSize(u_shadow_map, 0));
        for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <= 1; y++) {
                shadow += texture(u_shadow_map, vec3(proj.xy + vec2(x, y) * texel, proj.z - bias));
            }
        }
        shadow /= 9.0;
    }

    vec3 color = albedo * (ambient * ao + diffuse * shadow) + vec3(1.0) * (spec * shadow + rim);
    frag_color = vec4(color, 1.0);
}
