#include "ui/export_dialog.h"

#include "renderer/video_export.h"
#include "ui/file_dialog.h"

#include "ui/app.h"

#include <glad/gl.h>
#include <imgui.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

struct ResolutionPreset {
    const char* label;
    int width;
    int height;
};

constexpr std::array<ResolutionPreset, 3> kResolutionPresets = {{
    {"720p", 1280, 720},
    {"1080p", 1920, 1080},
    {"4K", 3840, 2160},
}};

constexpr std::array<const char*, 3> kFormats = {
    "MP4 (H.264)",
    "MOV (ProRes)",
    "GIF",
};

enum class ExportKind {
    None,
    Turntable,
    Trajectory,
    OrbitalSweep,
};

struct OffscreenTarget {
    unsigned int fbo = 0;
    unsigned int color = 0;
    unsigned int depth = 0;
    int width = 0;
    int height = 0;
};

struct CameraSnapshot {
    Eigen::Quaternionf orientation = Eigen::Quaternionf::Identity();
    Eigen::Vector3f target = Eigen::Vector3f::Zero();
    float distance = 15.0f;
};

struct ExportDialogState {
    int resolution_index = 1;
    int fps = 30;
    float duration = 12.0f;
    int crf = 18;
    int format_index = 0;
    bool transparent_background = false;

    int trajectory_resolution_index = 1;
    int trajectory_fps = 30;
    float trajectory_duration = 12.0f;
    int trajectory_crf = 18;
    int trajectory_format_index = 0;
    bool trajectory_transparent = false;

    int sweep_resolution_index = 1;
    int sweep_fps = 30;
    int sweep_crf = 18;
    int sweep_format_index = 0;
    bool sweep_transparent = false;
    int sweep_start_mo = 0;
    int sweep_end_mo = 0;
    int frames_per_orbital = 12;

    sbox::render::VideoExporter exporter;
    sbox::render::VideoExporter::ExportSettings active_settings;
    ExportKind active_kind = ExportKind::None;
    OffscreenTarget target;
    std::vector<unsigned char> pixels;
    int frame_index = 0;
    bool exporting = false;
    bool popup_open = false;

    CameraSnapshot camera_snapshot;
    sbox::chem::MolecularSystem molecule_snapshot;
    int selected_mo_snapshot = -1;
};

ExportDialogState& dialog_state() {
    static ExportDialogState state;
    return state;
}

void destroy_target(OffscreenTarget& target) {
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
    target.width = 0;
    target.height = 0;
}

bool ensure_target(OffscreenTarget& target, int width, int height) {
    if (target.fbo != 0U && target.width == width && target.height == height) {
        return true;
    }

    destroy_target(target);

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
    const bool complete = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    if (!complete) {
        destroy_target(target);
        return false;
    }

    target.width = width;
    target.height = height;
    return true;
}

void flip_rows(std::vector<unsigned char>& pixels, int width, int height) {
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

void fill_settings_for_format(sbox::render::VideoExporter::ExportSettings& settings, int format_index) {
    if (format_index == 1) {
        settings.codec = "prores";
        settings.pixel_format = "yuv422p10le";
    } else if (format_index == 2) {
        settings.codec = "gif";
        settings.pixel_format = "rgba";
    } else {
        settings.codec = "libx264";
        settings.pixel_format = "yuv420p";
    }
}

std::string default_filename(const char* stem, int format_index) {
    if (format_index == 1) {
        return std::string(stem) + ".mov";
    }
    if (format_index == 2) {
        return std::string(stem) + ".gif";
    }
    return std::string(stem) + ".mp4";
}

bool start_export(ExportDialogState& dialog,
                  sbox::App& app,
                  ExportKind kind,
                  const sbox::render::VideoExporter::ExportSettings& settings) {
    if (!ensure_target(dialog.target, settings.width, settings.height)) {
        return false;
    }
    if (!dialog.exporter.begin(settings)) {
        return false;
    }

    dialog.active_kind = kind;
    dialog.active_settings = settings;
    dialog.frame_index = 0;
    dialog.exporting = true;
    dialog.pixels.assign(static_cast<std::size_t>(settings.width) * static_cast<std::size_t>(settings.height) * 4u, 0);
    dialog.camera_snapshot.orientation = app.camera().orientation();
    dialog.camera_snapshot.target = app.camera().target();
    dialog.camera_snapshot.distance = app.camera().distance();
    dialog.molecule_snapshot = app.current_molecule();
    dialog.selected_mo_snapshot = app.state().selected_mo;
    return true;
}

void restore_export_state(ExportDialogState& dialog, sbox::App& app) {
    app.camera().setOrientation(dialog.camera_snapshot.orientation);
    app.camera().setTarget(dialog.camera_snapshot.target);
    app.camera().setDistance(dialog.camera_snapshot.distance);
    app.set_current_molecule_for_export(dialog.molecule_snapshot);
    app.state().selected_mo = dialog.selected_mo_snapshot;
}

void finish_export(ExportDialogState& dialog, sbox::App& app) {
    restore_export_state(dialog, app);
    dialog.exporter.end();
    dialog.exporting = false;
    dialog.active_kind = ExportKind::None;
    dialog.frame_index = 0;
}

bool process_one_export_frame(ExportDialogState& dialog, sbox::App& app) {
    const auto& settings = dialog.active_settings;
    if (dialog.frame_index >= settings.total_frames) {
        return false;
    }

    switch (dialog.active_kind) {
    case ExportKind::Turntable: {
        const float t = settings.total_frames > 1
                            ? static_cast<float>(dialog.frame_index) / static_cast<float>(settings.total_frames - 1)
                            : 0.0f;
        const float angle = 360.0f * t;
        app.camera().setOrientation(Eigen::AngleAxisf(angle * 3.14159265358979323846f / 180.0f,
                                                      Eigen::Vector3f::UnitY()) *
                                    dialog.camera_snapshot.orientation);
        break;
    }
    case ExportKind::Trajectory: {
        const sbox::io::Trajectory& trajectory = app.current_trajectory();
        if (trajectory.empty()) {
            return false;
        }
        const double t = settings.total_frames > 1
                             ? (static_cast<double>(dialog.frame_index) / static_cast<double>(settings.total_frames - 1)) *
                                   static_cast<double>(std::max(trajectory.num_frames() - 1, 0))
                             : 0.0;
        app.set_current_molecule_for_export(trajectory.interpolate(t));
        break;
    }
    case ExportKind::OrbitalSweep: {
        const int num_orbitals = std::max(1, dialog.sweep_end_mo - dialog.sweep_start_mo + 1);
        const int orbit_index = std::min(num_orbitals - 1, dialog.frame_index / std::max(1, dialog.frames_per_orbital));
        app.state().selected_mo = dialog.sweep_start_mo + orbit_index;
        break;
    }
    case ExportKind::None:
        return false;
    }

    app.render_to_fbo(dialog.target.fbo, settings.width, settings.height, settings.transparent_background);
    glBindFramebuffer(GL_FRAMEBUFFER, dialog.target.fbo);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, settings.width, settings.height, GL_RGBA, GL_UNSIGNED_BYTE, dialog.pixels.data());
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    flip_rows(dialog.pixels, settings.width, settings.height);

    if (!dialog.exporter.write_frame(dialog.pixels.data(), settings.width, settings.height)) {
        return false;
    }

    ++dialog.frame_index;
    return true;
}

void draw_resolution_combo(const char* label, int& index) {
    if (ImGui::BeginCombo(label, kResolutionPresets[static_cast<std::size_t>(index)].label)) {
        for (int i = 0; i < static_cast<int>(kResolutionPresets.size()); ++i) {
            const bool selected = i == index;
            if (ImGui::Selectable(kResolutionPresets[static_cast<std::size_t>(i)].label, selected)) {
                index = i;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
}

void draw_format_combo(const char* label, int& index) {
    if (ImGui::BeginCombo(label, kFormats[static_cast<std::size_t>(index)])) {
        for (int i = 0; i < static_cast<int>(kFormats.size()); ++i) {
            const bool selected = i == index;
            if (ImGui::Selectable(kFormats[static_cast<std::size_t>(i)], selected)) {
                index = i;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
}

}  // namespace

void draw_export_dialog(bool& show, sbox::App& app) {
    ExportDialogState& dialog = dialog_state();
    if (show && !dialog.popup_open) {
        ImGui::OpenPopup("Export Animation");
        dialog.popup_open = true;
        if (app.state().num_mo > 0) {
            const int homo = app.state().homo_index >= 0 ? app.state().homo_index : 0;
            dialog.sweep_start_mo = std::max(0, homo - 3);
            dialog.sweep_end_mo = std::min(std::max(0, app.state().num_mo - 1), homo + 4);
        }
    }

    bool open = show || dialog.exporting;
    if (dialog.exporting) {
        if (!process_one_export_frame(dialog, app) || dialog.frame_index >= dialog.active_settings.total_frames) {
            finish_export(dialog, app);
        }
    }

    if (ImGui::BeginPopupModal("Export Animation", &open, ImGuiWindowFlags_AlwaysAutoResize)) {
        const bool ffmpeg_available = sbox::render::VideoExporter::is_ffmpeg_available();

        if (!ffmpeg_available) {
            ImGui::TextWrapped("FFmpeg not found. Install with: brew install ffmpeg (macOS) or apt install ffmpeg (Linux)");
        }

        if (!dialog.exporting) {
            if (ImGui::BeginTabBar("##export_tabs")) {
                if (ImGui::BeginTabItem("Turntable")) {
                    ImGui::TextUnformatted("Rotate the molecule 360 degrees and save as video.");
                    draw_resolution_combo("Resolution", dialog.resolution_index);
                    ImGui::SliderInt("FPS", &dialog.fps, 15, 60);
                    ImGui::SliderFloat("Duration (s)", &dialog.duration, 2.0f, 30.0f, "%.1f");
                    ImGui::SliderInt("CRF", &dialog.crf, 10, 35);
                    draw_format_combo("Format", dialog.format_index);
                    ImGui::Checkbox("Transparent background (MOV only)", &dialog.transparent_background);
                    if (dialog.format_index != 1) {
                        dialog.transparent_background = false;
                    }

                    if (!ffmpeg_available) {
                        ImGui::BeginDisabled();
                    }
                    if (ImGui::Button("Export Turntable")) {
                        auto settings = sbox::render::VideoExporter::ExportSettings{};
                        const auto& preset = kResolutionPresets[static_cast<std::size_t>(dialog.resolution_index)];
                        settings.width = preset.width;
                        settings.height = preset.height;
                        settings.fps = dialog.fps;
                        settings.total_frames = static_cast<int>(std::round(dialog.duration * static_cast<float>(dialog.fps)));
                        settings.crf = dialog.crf;
                        settings.transparent_background = dialog.transparent_background;
                        fill_settings_for_format(settings, dialog.format_index);
                        const std::string path = save_file_dialog("Export Turntable",
                                                                  dialog.format_index == 1 ? "mov" : (dialog.format_index == 2 ? "gif" : "mp4"),
                                                                  default_filename("turntable", dialog.format_index).c_str());
                        if (!path.empty()) {
                            settings.output_path = path;
                            start_export(dialog, app, ExportKind::Turntable, settings);
                        }
                    }
                    if (!ffmpeg_available) {
                        ImGui::EndDisabled();
                    }
                    ImGui::EndTabItem();
                }

                if (ImGui::BeginTabItem("Trajectory")) {
                    ImGui::BeginDisabled(!app.has_trajectory());
                    draw_resolution_combo("Resolution", dialog.trajectory_resolution_index);
                    ImGui::SliderInt("FPS", &dialog.trajectory_fps, 15, 60);
                    ImGui::SliderFloat("Duration (s)", &dialog.trajectory_duration, 2.0f, 30.0f, "%.1f");
                    ImGui::SliderInt("CRF", &dialog.trajectory_crf, 10, 35);
                    draw_format_combo("Format", dialog.trajectory_format_index);
                    ImGui::Checkbox("Transparent background (MOV only)##traj", &dialog.trajectory_transparent);
                    if (dialog.trajectory_format_index != 1) {
                        dialog.trajectory_transparent = false;
                    }
                    const int total_frames = static_cast<int>(std::round(dialog.trajectory_duration *
                                                                         static_cast<float>(dialog.trajectory_fps)));
                    ImGui::Text("Frames: %d trajectory frames, interpolated to %d video frames",
                                app.current_trajectory().num_frames(),
                                total_frames);
                    if (!ffmpeg_available) {
                        ImGui::BeginDisabled();
                    }
                    if (ImGui::Button("Export Trajectory")) {
                        auto settings = sbox::render::VideoExporter::ExportSettings{};
                        const auto& preset = kResolutionPresets[static_cast<std::size_t>(dialog.trajectory_resolution_index)];
                        settings.width = preset.width;
                        settings.height = preset.height;
                        settings.fps = dialog.trajectory_fps;
                        settings.total_frames = total_frames;
                        settings.crf = dialog.trajectory_crf;
                        settings.transparent_background = dialog.trajectory_transparent;
                        fill_settings_for_format(settings, dialog.trajectory_format_index);
                        const std::string path = save_file_dialog("Export Trajectory",
                                                                  dialog.trajectory_format_index == 1 ? "mov" : (dialog.trajectory_format_index == 2 ? "gif" : "mp4"),
                                                                  default_filename("trajectory", dialog.trajectory_format_index).c_str());
                        if (!path.empty()) {
                            settings.output_path = path;
                            start_export(dialog, app, ExportKind::Trajectory, settings);
                        }
                    }
                    if (!ffmpeg_available) {
                        ImGui::EndDisabled();
                    }
                    ImGui::EndDisabled();
                    ImGui::EndTabItem();
                }

                if (ImGui::BeginTabItem("Orbital Sweep")) {
                    ImGui::BeginDisabled(!app.has_mo_data());
                    draw_resolution_combo("Resolution", dialog.sweep_resolution_index);
                    ImGui::SliderInt("FPS", &dialog.sweep_fps, 15, 60);
                    ImGui::SliderInt("CRF", &dialog.sweep_crf, 10, 35);
                    draw_format_combo("Format", dialog.sweep_format_index);
                    ImGui::Checkbox("Transparent background (MOV only)##sweep", &dialog.sweep_transparent);
                    if (dialog.sweep_format_index != 1) {
                        dialog.sweep_transparent = false;
                    }
                    if (app.state().num_mo > 0) {
                        ImGui::SliderInt("Sweep From MO", &dialog.sweep_start_mo, 0, app.state().num_mo - 1);
                        ImGui::SliderInt("Sweep To MO", &dialog.sweep_end_mo, 0, app.state().num_mo - 1);
                        if (dialog.sweep_start_mo > dialog.sweep_end_mo) {
                            std::swap(dialog.sweep_start_mo, dialog.sweep_end_mo);
                        }
                    }
                    ImGui::SliderInt("Frames per orbital", &dialog.frames_per_orbital, 5, 60);
                    if (!ffmpeg_available) {
                        ImGui::BeginDisabled();
                    }
                    if (ImGui::Button("Export Orbital Sweep")) {
                        auto settings = sbox::render::VideoExporter::ExportSettings{};
                        const auto& preset = kResolutionPresets[static_cast<std::size_t>(dialog.sweep_resolution_index)];
                        settings.width = preset.width;
                        settings.height = preset.height;
                        settings.fps = dialog.sweep_fps;
                        settings.total_frames = std::max(1, dialog.sweep_end_mo - dialog.sweep_start_mo + 1) * dialog.frames_per_orbital;
                        settings.crf = dialog.sweep_crf;
                        settings.transparent_background = dialog.sweep_transparent;
                        fill_settings_for_format(settings, dialog.sweep_format_index);
                        const std::string path = save_file_dialog("Export Orbital Sweep",
                                                                  dialog.sweep_format_index == 1 ? "mov" : (dialog.sweep_format_index == 2 ? "gif" : "mp4"),
                                                                  default_filename("orbital_sweep", dialog.sweep_format_index).c_str());
                        if (!path.empty()) {
                            settings.output_path = path;
                            start_export(dialog, app, ExportKind::OrbitalSweep, settings);
                        }
                    }
                    if (!ffmpeg_available) {
                        ImGui::EndDisabled();
                    }
                    ImGui::EndDisabled();
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }
        } else {
            ImGui::TextUnformatted("Export in progress");
            ImGui::ProgressBar(dialog.exporter.progress(), ImVec2(360.0f, 0.0f));
            ImGui::Text("Frame %d / %d", dialog.exporter.frames_written(), dialog.exporter.total_frames());
            if (ImGui::Button("Cancel")) {
                finish_export(dialog, app);
            }
        }

        if (!dialog.exporting) {
            if (ImGui::Button("Close")) {
                ImGui::CloseCurrentPopup();
            }
        }
        ImGui::EndPopup();
    } else {
        dialog.popup_open = false;
    }

    if (!open && !dialog.exporting) {
        show = false;
        dialog.popup_open = false;
        destroy_target(dialog.target);
    } else {
        show = true;
    }
}

}  // namespace sbox::ui
