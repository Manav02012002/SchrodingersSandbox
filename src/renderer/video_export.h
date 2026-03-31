#pragma once

#include "io/trajectory_io.h"

#include <cstdio>
#include <string>

namespace sbox {
class App;
}

namespace sbox::render {

class VideoExporter {
public:
    VideoExporter();
    ~VideoExporter();

    static bool is_ffmpeg_available();

    struct ExportSettings {
        std::string output_path;
        int width = 1920;
        int height = 1080;
        int fps = 30;
        int total_frames = 360;
        std::string codec = "libx264";
        int crf = 18;
        std::string pixel_format = "yuv420p";
        bool transparent_background = false;
    };

    bool begin(const ExportSettings& settings);
    bool write_frame(const unsigned char* rgba_data, int width, int height);
    bool end();

    int frames_written() const;
    int total_frames() const;
    bool is_exporting() const;
    float progress() const;

private:
    FILE* ffmpeg_pipe_ = nullptr;
    ExportSettings settings_;
    int frames_written_ = 0;
    bool exporting_ = false;
};

void export_turntable(
    sbox::App& app,
    const VideoExporter::ExportSettings& settings,
    float start_angle = 0.0f,
    float end_angle = 360.0f);

void export_trajectory(
    sbox::App& app,
    const sbox::io::Trajectory& trajectory,
    const VideoExporter::ExportSettings& settings);

void export_orbital_sweep(
    sbox::App& app,
    const VideoExporter::ExportSettings& settings,
    int start_mo,
    int end_mo);

}  // namespace sbox::render
