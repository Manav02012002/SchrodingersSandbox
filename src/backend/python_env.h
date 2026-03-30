#pragma once

#include <string>

namespace sbox::backend {

struct PythonInfo {
    std::string python_path;
    std::string version;
    bool has_pyscf = false;
    std::string pyscf_version;
    bool has_tblite = false;
    std::string tblite_version;
    bool has_xtb = false;
    std::string xtb_version;
    bool has_geometric = false;
    bool has_ase = false;
    bool valid = false;
};

class PythonEnvironment {
public:
    PythonEnvironment();

    void detect();
    void check_packages();

    const PythonInfo& info() const;
    bool is_valid() const;
    bool has_pyscf() const;
    bool has_tblite() const;

    void save_preference() const;
    void set_python_path(const std::string& path);

    static int run_capture(const std::string& command, std::string& stdout_out);

private:
    PythonInfo info_;
    bool test_python(const std::string& path);
};

}  // namespace sbox::backend
