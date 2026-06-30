#include "cxxopts.hpp"
#include "file_format.hpp"
#include "vg.pb.h" // <--- THIS LINE IS CRITICAL
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    // Required for Protobuf library cleanup
    google::protobuf::ShutdownProtobufLibrary();

    cxxopts::Options options("gam_processor", "C++ version of the GAM segment extractor.");

    options.add_options()
        ("g,gam_path", "Path to the GAM file", cxxopts::value<std::string>())
        ("s,stats_path", "Path to the stats TSV file", cxxopts::value<std::string>())
        ("o,output_prefix", "Prefix for output files (.dat, .idx)", cxxopts::value<std::string>())
        ("m,milestone", "Progress report interval", cxxopts::value<int>()->default_value("1000000"))
        ("c,chr", "Optional chromosome name to filter on", cxxopts::value<std::string>()->default_value(""))
        ("e,use_existing", "Reuse existing initialized output", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");

    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (!result.count("gam_path") || !result.count("stats_path") || !result.count("output_prefix")) {
            std::cerr << "Error: Missing required arguments: --gam_path, --stats_path, --output_prefix" << std::endl;
            std::cerr << options.help() << std::endl;
            return 1;
        }

        run_pipeline(
            result["gam_path"].as<std::string>(),
            result["stats_path"].as<std::string>(),
            result["output_prefix"].as<std::string>(),
            result["milestone"].as<int>(),
            result["chr"].as<std::string>(),
            result["use_existing"].as<bool>()
        );

    } catch (const cxxopts::OptionException& e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        return 1;
    } catch (const std::runtime_error& e) {
        std::cerr << "A runtime error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}