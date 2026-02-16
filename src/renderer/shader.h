#pragma once

#include <Eigen/Core>

#include <string>

namespace sbox {

class Shader {
public:
    Shader(const std::string& vertex_path, const std::string& fragment_path);
    ~Shader();

    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;

    void bind() const;
    [[nodiscard]] unsigned int id() const { return program_id_; }

    void setUniform(const std::string& name, int value) const;
    void setUniform(const std::string& name, float value) const;
    void setUniform(const std::string& name, const Eigen::Vector3f& value) const;
    void setUniform(const std::string& name, const Eigen::Matrix4f& value) const;

private:
    static std::string LoadFile(const std::string& path);
    static unsigned int Compile(unsigned int type, const std::string& source, const std::string& path);

    unsigned int program_id_ = 0;
};

}  // namespace sbox
