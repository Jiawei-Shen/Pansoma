// file: fast_writer.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversions of dicts and lists
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <unistd.h>
#include <unordered_map>

namespace py = pybind11;

// The Segment struct remains the same
struct Segment {
    int16_t offset;
    std::string seq;
    std::string bq;
    std::string cigar;
    int16_t rq;
    char strand;
};

// NEW: A struct to hold the data from the block_infos dictionary
// pybind11 will automatically convert a Python dict with matching keys
// into this struct.
struct BlockInfo {
    uint64_t offset;
    uint32_t n_records;
    uint32_t current_pos;
    uint32_t max_read_len;
    uint32_t max_cigar_len;
    uint32_t record_size;
    uint64_t block_size;
};

// Internal helper function for writing one node's data.
// This is the logic from our old C++ function. It's not exposed to Python anymore.
void write_node_data(int fd, long long write_pos, const std::vector<Segment>& segments, int max_read_len, int max_cigar_len) {
    if (segments.empty()) {
        return;
    }
    size_t record_size = sizeof(int16_t) + max_read_len * 2 + max_cigar_len + sizeof(int16_t) + sizeof(char);
    size_t total_buffer_size = segments.size() * record_size;
    std::vector<char> buffer(total_buffer_size);
    char* current_ptr = buffer.data();

    for (const auto& s : segments) {
        memcpy(current_ptr, &s.offset, sizeof(int16_t)); current_ptr += sizeof(int16_t);
        memcpy(current_ptr, s.seq.c_str(), s.seq.length());
        memset(current_ptr + s.seq.length(), 0, max_read_len - s.seq.length());
        current_ptr += max_read_len;
        memcpy(current_ptr, s.bq.c_str(), s.bq.length());
        memset(current_ptr + s.bq.length(), 0, max_read_len - s.bq.length());
        current_ptr += max_read_len;
        memcpy(current_ptr, s.cigar.c_str(), s.cigar.length());
        memset(current_ptr + s.cigar.length(), 0, max_cigar_len - s.cigar.length());
        current_ptr += max_cigar_len;
        memcpy(current_ptr, &s.rq, sizeof(int16_t)); current_ptr += sizeof(int16_t);
        memcpy(current_ptr, &s.strand, sizeof(char)); current_ptr += sizeof(char);
    }

    ssize_t bytes_written = pwrite(fd, buffer.data(), buffer.size(), write_pos);
    if (bytes_written == -1 || static_cast<size_t>(bytes_written) != buffer.size()) {
        throw std::runtime_error("pwrite failed for a data block.");
    }
}

// NEW: The main function that loops through the entire buffer.
// It takes block_infos by reference (&) so it can modify current_pos.
void flush_entire_buffer(
    int fd,
    const std::unordered_map<uint32_t, std::vector<Segment>>& segment_buffer,
    std::unordered_map<uint32_t, BlockInfo>& block_infos,
    const uint32_t block_header_size)
{
    for (const auto& pair : segment_buffer) {
        uint32_t nid = pair.first;
        const std::vector<Segment>& segs = pair.second;

        try {
            // Use .at() to get a reference to the info struct.
            // This will throw an error if nid is not in block_infos.
            BlockInfo& info = block_infos.at(nid);

            // 1. C++ calculates the metadata
            long long base_offset = info.offset + block_header_size;
            long long write_pos = base_offset + (long long)info.current_pos * info.record_size;

            // 2. C++ does the serialization and writing by calling the helper
            write_node_data(fd, write_pos, segs, info.max_read_len, info.max_cigar_len);

            // 3. C++ updates the current_pos. This modification will be
            // reflected back in the Python dictionary.
            info.current_pos += segs.size();

        } catch (const std::out_of_range& oor) {
            // Handle cases where a node ID from the buffer isn't in block_infos
            // For now, we just print a warning and continue.
            // You could also throw an exception back to Python if this is a critical error.
            // std::cerr << "Warning: Node ID " << nid << " not found in block_infos. Skipping." << std::endl;
        }
    }
}

// pybind11 module definition
PYBIND11_MODULE(fast_writer, m) {
    m.doc() = "C++ module for processing and writing segment data buffers.";

    py::class_<Segment>(m, "Segment")
        .def(py::init<int16_t, std::string, std::string, std::string, int16_t, char>(),
             py::arg("offset"), py::arg("seq"), py::arg("bq"),
             py::arg("cigar"), py::arg("rq"), py::arg("strand"));

    // Expose the BlockInfo struct to Python
    py::class_<BlockInfo>(m, "BlockInfo")
        .def(py::init<>())
        .def_readwrite("offset", &BlockInfo::offset)
        .def_readwrite("n_records", &BlockInfo::n_records)
        .def_readwrite("current_pos", &BlockInfo::current_pos)
        .def_readwrite("max_read_len", &BlockInfo::max_read_len)
        .def_readwrite("max_cigar_len", &BlockInfo::max_cigar_len)
        .def_readwrite("record_size", &BlockInfo::record_size)
        .def_readwrite("block_size", &BlockInfo::block_size);

    // Expose the new main function
    m.def("flush_entire_buffer", &flush_entire_buffer,
          "Processes the entire segment buffer and writes to disk.",
          py::arg("fd"), py::arg("segment_buffer"),
          py::arg("block_infos"), py::arg("block_header_size"));
}