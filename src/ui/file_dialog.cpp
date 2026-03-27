#include "ui/file_dialog.h"

#include <nfd.h>

#include <string>
#include <vector>

namespace sbox::ui {
namespace {

std::vector<std::string> split_filters(const char* filters) {
    std::vector<std::string> parts;
    if (filters == nullptr) {
        return parts;
    }

    std::string current;
    for (const char* p = filters; *p != '\0'; ++p) {
        const char ch = *p;
        if (ch == ',' || ch == ';') {
            if (!current.empty()) {
                parts.push_back(current);
                current.clear();
            }
        } else if (ch != ' ' && ch != '\t') {
            current.push_back(ch);
        }
    }
    if (!current.empty()) {
        parts.push_back(current);
    }
    return parts;
}

}  // namespace

std::string open_file_dialog(const char* title, const char* filters) {
    (void)title;

    if (NFD_Init() != NFD_OKAY) {
        return "";
    }

    std::string selected_path;
    nfdu8char_t* path = nullptr;
    const std::vector<std::string> extensions = split_filters(filters);
    std::vector<nfdu8filteritem_t> items;
    items.reserve(extensions.size());
    for (const std::string& ext : extensions) {
        items.push_back({ext.c_str(), ext.c_str()});
    }

    const nfdresult_t result =
        NFD_OpenDialogU8(&path, items.empty() ? nullptr : items.data(), items.size(), nullptr);
    if (result == NFD_OKAY && path != nullptr) {
        selected_path = path;
        NFD_FreePathU8(path);
    }

    NFD_Quit();
    return selected_path;
}

}  // namespace sbox::ui
