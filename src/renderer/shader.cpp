#include "renderer/shader.h"

#include <glad/gl.h>

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace sbox {

Shader::Shader(const std::string& vertex_path, const std::string& fragment_path) {
    const std::string vertex_src = LoadFile(vertex_path);
    const std::string fragment_src = LoadFile(fragment_path);

    const unsigned int vertex_shader = Compile(GL_VERTEX_SHADER, vertex_src, vertex_path);
    const unsigned int fragment_shader = Compile(GL_FRAGMENT_SHADER, fragment_src, fragment_path);

    program_id_ = glCreateProgram();
    glAttachShader(program_id_, vertex_shader);
    glAttachShader(program_id_, fragment_shader);
    glLinkProgram(program_id_);

    int success = 0;
    glGetProgramiv(program_id_, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
        int log_len = 0;
        glGetProgramiv(program_id_, GL_INFO_LOG_LENGTH, &log_len);
        std::vector<char> log(static_cast<std::size_t>(log_len));
        glGetProgramInfoLog(program_id_, log_len, nullptr, log.data());
        std::cerr << "Shader link error: " << log.data() << '\n';

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        glDeleteProgram(program_id_);
        program_id_ = 0;
        throw std::runtime_error("Failed to link shader program");
    }

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
}

Shader::~Shader() {
    if (program_id_ != 0) {
        glDeleteProgram(program_id_);
    }
}

void Shader::bind() const {
    glUseProgram(program_id_);
}

void Shader::setUniform(const std::string& name, int value) const {
    const int loc = glGetUniformLocation(program_id_, name.c_str());
    if (loc >= 0) {
        glUniform1i(loc, value);
    }
}

void Shader::setUniform(const std::string& name, float value) const {
    const int loc = glGetUniformLocation(program_id_, name.c_str());
    if (loc >= 0) {
        glUniform1f(loc, value);
    }
}

void Shader::setUniform(const std::string& name, const Eigen::Vector3f& value) const {
    const int loc = glGetUniformLocation(program_id_, name.c_str());
    if (loc >= 0) {
        glUniform3f(loc, value.x(), value.y(), value.z());
    }
}

void Shader::setUniform(const std::string& name, const Eigen::Matrix4f& value) const {
    const int loc = glGetUniformLocation(program_id_, name.c_str());
    if (loc >= 0) {
        glUniformMatrix4fv(loc, 1, GL_FALSE, value.data());
    }
}

std::string Shader::LoadFile(const std::string& path) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open shader file: " + path);
    }

    std::ostringstream ss;
    ss << input.rdbuf();
    return ss.str();
}

unsigned int Shader::Compile(unsigned int type, const std::string& source, const std::string& path) {
    const unsigned int shader = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    int success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        int log_len = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_len);
        std::vector<char> log(static_cast<std::size_t>(log_len));
        glGetShaderInfoLog(shader, log_len, nullptr, log.data());
        std::cerr << "Shader compile error in " << path << ": " << log.data() << '\n';
        glDeleteShader(shader);
        throw std::runtime_error("Failed to compile shader: " + path);
    }

    return shader;
}

}  // namespace sbox
