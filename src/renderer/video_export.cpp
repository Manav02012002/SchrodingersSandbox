#include "renderer/video_export.h"

#include "ui/app.h"

#include <glad/gl.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace sbox::render {

namespace {

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char ch : value) {
        if (ch == '\'') {
            quoted += "'\\''";
        } else {
            quoted.push_back(ch);
        }
    }
    quoted.push_back('\'');
    return quoted;
}

FILE* open_pipe(const std::string& command) {
#ifdef _WIN32
    return _popen(command.c_str(), "wb");
#else
    return popen(command.c_str(), "w");
#endif
}

int close_pipe(FILE* pipe) {
#ifdef _WIN32
    return _pclose(pipe);
#else
    return pclose(pipe);
#endif
}

void flip_rgba_rows(std::vector<unsigned char>& pixels, int width, int height) {
    const std::size_t row_bytes = static_cast<std::size_t>(width) * 4u;
    std::vector<unsigned char> temp(row_bytes);
    for (int y = 0; y < height / 2; ++y) {
        unsigned char* top = pixels.data() + static_cast<std::size_t>(y) * row_bytes;
        unsigned char* bottom = pixels.data() + static_cast<std::size_t>(height - 1 - y) * row_bytes;
        std::copy(top, top + row_bytes, temp.data());
        std::copy(bottom, bottom + row_bytes, top);
        std::copy(temp.data(), temp.data() + row_bytes, bottom);
    }
}

struct OffscreenTarget {
    unsigned int fbo = 0;
    unsigned int color = 0;
    unsigned int depth = 0;
};

OffscreenTarget create_offscreen_target(int width, int height) {
    OffscreenTarget target;
    glGenFramebuffers(1, &target.fbo);
    glGenTextures(1, &target.color);
    glGenRenderbuffers(1, &target.depth);

    glBindTexture(GL_TEXTURE_2D, target.color);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindRenderbuffer(GL_RENDERBUFFER, target.depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);

    glBindFramebuffer(GL_FRAMEBUFFER, target.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, target.color, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, target.depth);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    return target;
}

void destroy_offscreen_target(OffscreenTarget& target) {
    if (target.depth != 0U) {
        glDeleteRenderbuffers(1, &target.depth);
        target.depth = 0;
    }
    if (target.color != 0U) {
        glDeleteTextures(1, &target.color);
        target.color = 0;
    }
    if (target.fbo != 0U) {
        glDeleteFramebuffers(1, &target.fbo);
        target.fbo = 0;
    }
}

bool capture_frame(sbox::App& app,
                   const VideoExporter::ExportSettings& settings,
                   const OffscreenTarget& target,
                   std::vector<unsigned char>& pixels,
                   VideoExporter& exporter) {
    app.render_to_fbo(target.fbo, settings.width, settings.height, settings.transparent_background);
    glBindFramebuffer(GL_FRAMEBUFFER, target.fbo);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, settings.width, settings.height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    flip_rgba_rows(pixels, settings.width, settings.height);
    return exporter.write_frame(pixels.data(), settings.width, settings.height);
}

}  // namespace

VideoExporter::VideoExporter() = default;

VideoExporter::~VideoExporter() {
    end();
}

bool VideoExporter::is_ffmpeg_available() {
    FILE* pipe = open_pipe("ffmpeg -version >/dev/null 2>&1");
    if (pipe == nullptr) {
        return false;
    }
    const int status = close_pipe(pipe);
    return status == 0;
}

bool VideoExporter::begin(const ExportSettings& settings) {
    end();

    settings_ = settings;
    frames_written_ = 0;

    std::string codec = settings_.codec;
    std::string pixel_format = settings_.pixel_format;
    std::string command;

    if (settings_.transparent_background) {
        codec = "prores_ks";
        pixel_format = "yuva444p10le";
    }

    if (codec == "prores") {
        codec = "prores_ks";
    }

    if (codec == "gif") {
        command = "ffmpeg -y"
                  " -f rawvideo -pix_fmt rgba -s " + std::to_string(settings_.width) + "x" + std::to_string(settings_.height) +
                  " -r " + std::to_string(settings_.fps) +
                  " -i -"
                  " -vf " + shell_quote("fps=15,scale=480:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse") +
                  " " + shell_quote(settings_.output_path);
    } else {
        command = "ffmpeg -y"
                  " -f rawvideo -pix_fmt rgba -s " + std::to_string(settings_.width) + "x" + std::to_string(settings_.height) +
                  " -r " + std::to_string(settings_.fps) +
                  " -i -"
                  " -c:v " + codec;
        if (codec == "libx264" || codec == "libx265") {
            command += " -crf " + std::to_string(settings_.crf);
        }
        command += " -pix_fmt " + pixel_format +
                   " " + shell_quote(settings_.output_path);
    }

    ffmpeg_pipe_ = open_pipe(command);
    exporting_ = ffmpeg_pipe_ != nullptr;
    return exporting_;
}

bool VideoExporter::write_frame(const unsigned char* rgba_data, int width, int height) {
    if (!exporting_ || ffmpeg_pipe_ == nullptr || rgba_data == nullptr) {
        return false;
    }
    if (width != settings_.width || height != settings_.height) {
        return false;
    }

    const std::size_t bytes = static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u;
    const std::size_t written = std::fwrite(rgba_data, 1, bytes, ffmpeg_pipe_);
    if (written != bytes) {
        return false;
    }
    ++frames_written_;
    return true;
}

bool VideoExporter::end() {
    bool ok = true;
    if (ffmpeg_pipe_ != nullptr) {
        ok = close_pipe(ffmpeg_pipe_) == 0;
        ffmpeg_pipe_ = nullptr;
    }
    exporting_ = false;
    return ok;
}

int VideoExporter::frames_written() const {
    return frames_written_;
}

int VideoExporter::total_frames() const {
    return settings_.total_frames;
}

bool VideoExporter::is_exporting() const {
    return exporting_;
}

float VideoExporter::progress() const {
    if (settings_.total_frames <= 0) {
        return 0.0f;
    }
    return std::clamp(static_cast<float>(frames_written_) / static_cast<float>(settings_.total_frames), 0.0f, 1.0f);
}

void export_turntable(sbox::App& app,
                      const VideoExporter::ExportSettings& settings,
                      float start_angle,
                      float end_angle) {
    VideoExporter exporter;
    if (!exporter.begin(settings)) {
        return;
    }

    const Eigen::Quaternionf original_orientation = app.camera().orientation();
    const Eigen::Vector3f original_target = app.camera().target();
    const float original_distance = app.camera().distance();

    OffscreenTarget target = create_offscreen_target(settings.width, settings.height);
    std::vector<unsigned char> pixels(static_cast<std::size_t>(settings.width) * static_cast<std::size_t>(settings.height) * 4u);

    for (int frame = 0; frame < settings.total_frames; ++frame) {
        const float t = settings.total_frames > 1 ? static_cast<float>(frame) / static_cast<float>(settings.total_frames - 1) : 0.0f;
        const float angle_deg = start_angle + (end_angle - start_angle) * t;
        app.camera().setOrientation(Eigen::AngleAxisf(angle_deg * 3.14159265358979323846f / 180.0f,
                                                      Eigen::Vector3f::UnitY()) * original_orientation);
        if (!capture_frame(app, settings, target, pixels, exporter)) {
            break;
        }
    }

    app.camera().setOrientation(original_orientation);
    app.camera().setTarget(original_target);
    app.camera().setDistance(original_distance);
    destroy_offscreen_target(target);
    exporter.end();
}

void export_trajectory(sbox::App& app,
                       const sbox::io::Trajectory& trajectory,
                       const VideoExporter::ExportSettings& settings) {
    if (trajectory.empty()) {
        return;
    }

    VideoExporter exporter;
    if (!exporter.begin(settings)) {
        return;
    }

    const sbox::chem::MolecularSystem original_molecule = app.current_molecule();
    OffscreenTarget target = create_offscreen_target(settings.width, settings.height);
    std::vector<unsigned char> pixels(static_cast<std::size_t>(settings.width) * static_cast<std::size_t>(settings.height) * 4u);

    for (int frame = 0; frame < settings.total_frames; ++frame) {
        const double t = settings.total_frames > 1
                             ? (static_cast<double>(frame) / static_cast<double>(settings.total_frames - 1)) *
                                   static_cast<double>(std::max(trajectory.num_frames() - 1, 0))
                             : 0.0;
        app.set_current_molecule_for_export(trajectory.interpolate(t));
        if (!capture_frame(app, settings, target, pixels, exporter)) {
            break;
        }
    }

    app.set_current_molecule_for_export(original_molecule);
    destroy_offscreen_target(target);
    exporter.end();
}

void export_orbital_sweep(sbox::App& app,
                          const VideoExporter::ExportSettings& settings,
                          int start_mo,
                          int end_mo) {
    if (start_mo > end_mo) {
        std::swap(start_mo, end_mo);
    }

    VideoExporter exporter;
    if (!exporter.begin(settings)) {
        return;
    }

    const int original_mo = app.state().selected_mo;
    const int num_orbitals = std::max(1, end_mo - start_mo + 1);
    const int frames_per_orbital = std::max(1, settings.total_frames / num_orbitals);

    OffscreenTarget target = create_offscreen_target(settings.width, settings.height);
    std::vector<unsigned char> pixels(static_cast<std::size_t>(settings.width) * static_cast<std::size_t>(settings.height) * 4u);

    for (int frame = 0; frame < settings.total_frames; ++frame) {
        const int orbital_index = std::min(end_mo, start_mo + frame / frames_per_orbital);
        app.state().selected_mo = orbital_index;
        if (!capture_frame(app, settings, target, pixels, exporter)) {
            break;
        }
    }

    app.state().selected_mo = original_mo;
    destroy_offscreen_target(target);
    exporter.end();
}

}  // namespace sbox::render
