#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <atomic> // <--- ADD THIS LINE

// ─────────────────────────────────────────────────────────────────────────────
// Constants for the binary format
// ─────────────────────────────────────────────────────────────────────────────

const std::string GLOBAL_MAGIC = "MYFMT\x01";
constexpr uint8_t GLOBAL_MAJOR = 0;
constexpr uint8_t GLOBAL_MINOR = 5;

#pragma pack(push, 1)
struct GlobalHeader {
    uint8_t major_ver;
    uint8_t minor_ver;
    uint32_t block_count;
    char reserved[16];
};

struct BlockHeader {
    uint32_t node_id;
    uint32_t n_records;
    uint16_t flags;
    uint32_t max_read_length;
    uint32_t max_cigar_length;
};

struct IndexEntry {
    uint32_t node_id;
    uint64_t offset;
    uint32_t block_size;
    uint32_t n_records;
    uint16_t flags;
    uint32_t max_read_length;
    uint32_t max_cigar_length;
};
#pragma pack(pop)

// ─────────────────────────────────────────────────────────────────────────────
// Core data structures
// ─────────────────────────────────────────────────────────────────────────────

struct Segment {
    int16_t offset;
    std::string seq;
    std::string bq;
    std::string cigar;
    int16_t rq; // MAPQ
    char strand;
};

struct BlockInfo {
    uint64_t file_offset;
    uint32_t n_records;
    uint32_t max_read_len;
    uint32_t max_cigar_len;
    uint32_t record_size;
    uint32_t block_size;
    std::atomic<uint32_t> current_pos{0};
};

// ─────────────────────────────────────────────────────────────────────────────
// Function Prototypes
// ─────────────────────────────────────────────────────────────────────────────

void run_pipeline(
    const std::string& gam_path,
    const std::string& stats_path,
    const std::string& output_prefix,
    int milestone_step,
    const std::string& chrom_filter,
    bool use_existing
);