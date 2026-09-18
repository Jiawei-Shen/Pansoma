// Query distinct GBWT paths in a GBZ; one process serves all tensor batches.
#include <gbwtgraph/gbz.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

int main(int argc, char** argv) {
    using nlohmann::json;
    try {
        if (argc != 2) throw std::runtime_error("usage: gbz_node_counts graph.gbz");
        std::ifstream input(argv[1], std::ios::binary);
        if (!input) throw std::runtime_error("cannot open GBZ");
        gbwtgraph::GBZ gbz;
        gbz.simple_sds_load(input);
        if (!gbz.index.bidirectional())
            throw std::runtime_error("GBZ occurrence queries require a bidirectional GBWT");
        std::cout << json({{"protocol_version", 1}, {"metric", "distinct-gbwt-paths-v1"},
            {"indexed_paths", gbz.index.sequences() / 2},
            {"bidirectional", true}}).dump() << std::endl;
        std::string line;
        while (std::getline(std::cin, line)) {
            auto request = json::parse(line);
            auto node = request.at("node_id").get<std::uint64_t>();
            if (node == 0 || !gbz.graph.has_node(node))
                throw std::runtime_error("node missing from GBZ: " + std::to_string(node));
            std::set<gbwt::size_type> paths;
            // GBWT stores both a path and its reverse complement. Strip that
            // orientation bit and deduplicate repeated visits and both searches.
            for (bool reverse : {false, true}) {
                auto state = gbz.index.find(gbwt::Node::encode(node, reverse));
                for (auto sequence : gbz.index.locate(state))
                    paths.insert(gbwt::Path::id(sequence));
            }
            std::cout << json({{"node_id", node}, {"count", paths.size()},
                {"sequence", gbz.graph.get_sequence(gbz.graph.get_handle(node, false))}}).dump()
                << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "GBZ count query failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
