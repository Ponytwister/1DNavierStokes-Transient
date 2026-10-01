#include <application_options.h>

#include <stdexcept>
#include <string>

namespace tsensor_workflow {

application_options parse_options(int argc, const char* const argv[])
{
    application_options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") {
            options.help = true;
        } else if (argument == "--database" || argument == "--output-dir") {
            if (i + 1 == argc || std::string(argv[i + 1]).empty()
                || std::string(argv[i + 1]).rfind("--", 0) == 0
                || std::string(argv[i + 1]) == "-h") {
                throw std::invalid_argument("Missing path after " + argument);
            }
            const std::filesystem::path value(argv[++i]);
            if (argument == "--database") {
                options.database = value;
            } else {
                options.output_directory = value;
            }
        } else {
            throw std::invalid_argument("Unknown argument: " + argument);
        }
    }
    return options;
}

} // namespace tsensor_workflow
