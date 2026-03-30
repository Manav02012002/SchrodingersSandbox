#pragma once

#include "backend/job_types.h"
#include "backend/python_env.h"

#include <atomic>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <sys/types.h>

namespace sbox::backend {

class BackendManager {
public:
    BackendManager();
    ~BackendManager();

    void init(const PythonEnvironment& env);

    int submit(const JobSpec& spec);

    bool is_running(int job_id) const;
    JobStatus status(int job_id) const;
    const JobResult* result(int job_id) const;
    std::string work_dir(int job_id) const;

    std::vector<int> poll_completed();

    struct Progress {
        std::string stage;
        int iteration = 0;
        int step = 0;
        int total = 0;
        double energy = 0.0;
        std::string message;
    };
    Progress get_progress(int job_id) const;

    void cancel(int job_id);
    void clear_job(int job_id);

    bool can_run_pyscf() const;
    bool can_run_xtb() const;

    std::string scripts_dir() const;

private:
    struct RunningJob {
        int job_id = 0;
        JobSpec spec;
        std::future<JobResult> future;
        std::string work_dir;
        std::atomic<bool> cancelled{false};
        pid_t pid = -1;
    };

    PythonEnvironment python_env_;
    std::string scripts_dir_;
    int next_job_id_ = 1;

    mutable std::mutex mutex_;
    std::map<int, std::unique_ptr<RunningJob>> running_jobs_;
    std::map<int, JobResult> completed_jobs_;
    std::vector<int> newly_completed_;

    JobResult run_job(const JobSpec& spec, const std::string& work_dir, std::atomic<bool>& cancelled);
    std::string create_work_dir(int job_id);
    void write_job_json(const JobSpec& spec, const std::string& work_dir);
    JobResult parse_result(const JobSpec& spec, const std::string& work_dir);
    Progress parse_progress(const std::string& work_dir) const;
    std::string driver_script(const JobSpec& spec) const;
};

}  // namespace sbox::backend
