#pragma once

#include <filesystem>

namespace tsensor_workflow {

struct application_options {
    std::filesystem::path database = "../../navier.db";
    std::filesystem::path output_directory = "..";
    bool help = false;
};

// Parse without opening files or changing the working directory. Relative paths
// are resolved against the launch directory by the terminal entry point.
application_options parse_options(int argc, const char* const argv[]);

} // namespace tsensor_workflow
