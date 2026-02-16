#include "ui/app.h"

#include <exception>
#include <iostream>

int main() {
    try {
        sbox::App app;
        app.run();
    } catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
