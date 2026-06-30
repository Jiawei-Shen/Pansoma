// g++ -O3 -march=native -DNDEBUG -std=c++17 gam_segments.cpp -o gam_segments \
//     -lprotobuf -lz -pthread
//
// Requires:
//   - Protobuf headers for vg Alignment (e.g., #include "vg.pb.h" or "vg/vg.pb.h")
//   - zlib (for .gz reading)
//   - nlohmann::json single-header (drop json.hpp next to this file)
//
// Usage:
//   ./gam_segments  \
//     --gam /path/to/file.gam[.gz] \
//     --stats-json /path/to/stats.json \
//     --out-prefix /path/to/output_prefix \
//     [--milestone 1000000] [--chr chr1] [--use-existing]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <chrono>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>

#include "json.hpp"              // nlohmann JSON header (https://github.com/nlohmann/json)
// Adjust include to your vg proto path:
#include "vg.pb.h"               // or "vg/vg.pb.h"

using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────────────────────
// File format constants (exactly match Python)
static constexpr const char* GLOBAL_MAGIC = "MYFMT\x01";  // 6 bytes
static constexpr uint8_t     GLOBAL_MAJOR = 0;
static constexpr uint8_t     GLOBAL_MINOR = 5;            // minor=5: no node_length
// GLOBAL_VER_PACK "<BBI16s" == 1+1+4+16 = 22 bytes
// GLOBAL_HEADER_SIZE = 6 + 22 = 28 bytes

#pragma pack(push,1)
struct GlobalVerPack {
    uint8_t  major;
    uint8_t  minor;
    uint32_t block_count;
    uint8_t  reserved[16];
};
static_assert(sizeof(GlobalVerPack) == 22, "GlobalVerPack size mismatch");

struct BlockHdr {        // "<I I H I I" == 4+4+2+4+4 = 18
    uint32_t node_id;
    uint32_t n_records;
    uint16_t flags;
    uint32_t max_read_len;
    uint32_t max_cigar_len;
};
static_assert(sizeof(BlockHdr) == 18, "BlockHdr size mismatch");
#pragma pack(pop)

static constexpr size_t GLOBAL_HEADER_SIZE = 28;
static constexpr size_t BLOCK_HDR_SIZE = sizeof(BlockHdr);

// Per-record layout (variable by node): "<h {R}s {R}s {C}s h c"
struct RecordLayout {
    int R;    // max_read_len
    int C;    // max_cigar_len
    size_t size; // total bytes per record
    RecordLayout(int r=0, int c=0) : R(r), C(c) {
        // <h + R bytes + R bytes + C bytes + <h + c>
        size = sizeof(int16_t) + (size_t)R + (size_t)R + (size_t)C + sizeof(int16_t) + sizeof(char);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Segment container
struct Segment {
    int32_t  offset;          // node offset
    std::string seq;          // upcased sequence (unpadded)
    std::string bq;           // quality bytes (len == seq length)
    std::string cigar;        // ASCII
    int16_t rq;               // MAPQ fits in int16
    char    strand;           // '+' or '-'
};

// Per-node write info
struct BlockInfo {
    uint64_t offset;      // absolute file offset of the block (header start)
    uint32_t n_records;   // total records reserved
    uint32_t current_pos; // how many already written
    uint32_t R;           // max_read_len
    uint32_t C;           // max_cigar_len
    uint32_t record_size;
    uint32_t block_size;  // header + n*record_size
};

// ─────────────────────────────────────────────────────────────────────────────
// GZip or Plain Reader abstraction
class Reader {
public:
    virtual ~Reader() = default;
    virtual bool read1(uint8_t& out) = 0;
    virtual bool readN(void* buf, size_t n) = 0;
    virtual bool skipN(size_t n) = 0;
    virtual bool eof() const = 0;
};

class FileReader : public Reader {
    int fd_;
    off_t size_;
    off_t pos_;
public:
    explicit FileReader(const std::string& path) {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) {
            perror("open");
            throw std::runtime_error("Failed to open " + path);
        }
        struct stat st{};
        if (fstat(fd_, &st) != 0) {
            perror("fstat");
            throw std::runtime_error("fstat failed");
        }
        size_ = st.st_size;
        pos_ = 0;
    }
    ~FileReader() override {
        if (fd_ >= 0) ::close(fd_);
    }
    bool read1(uint8_t& out) override {
        ssize_t r = ::read(fd_, &out, 1);
        if (r == 1) { pos_++; return true; }
        return false;
    }
    bool readN(void* buf, size_t n) override {
        uint8_t* p = static_cast<uint8_t*>(buf);
        size_t got = 0;
        while (got < n) {
            ssize_t r = ::read(fd_, p + got, n - got);
            if (r <= 0) return false;
            got += size_t(r);
        }
        pos_ += got;
        return true;
    }
    bool skipN(size_t n) override {
        off_t newpos = lseek(fd_, (off_t)n, SEEK_CUR);
        if (newpos == (off_t)-1) return false;
        pos_ = newpos;
        return true;
    }
    bool eof() const override { return pos_ >= size_; }
};

class GzReader : public Reader {
    gzFile gz_;
public:
    explicit GzReader(const std::string& path) {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) throw std::runtime_error("gzopen failed: " + path);
    }
    ~GzReader() override {
        if (gz_) gzclose(gz_);
    }
    bool read1(uint8_t& out) override {
        int c = gzgetc(gz_);
        if (c == -1) return false;
        out = static_cast<uint8_t>(c);
        return true;
    }
    bool readN(void* buf, size_t n) override {
        size_t got = 0;
        while (got < n) {
            int r = gzread(gz_, (uint8_t*)buf + got, (unsigned int)(n - got));
            if (r <= 0) return false;
            got += (size_t)r;
        }
        return true;
    }
    bool skipN(size_t n) override {
        // gzseek supports forward seeks
        z_off_t res = gzseek(gz_, (z_off_t)n, SEEK_CUR);
        return res != -1;
    }
    bool eof() const override {
        // Heuristic: gzgets returns -1 at EOF; no direct eof()
        return gzeof(gz_) != 0;
    }
};

static bool file_is_gzip(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    unsigned char b0=0, b1=0;
    f.read((char*)&b0, 1);
    f.read((char*)&b1, 1);
    return b0==0x1f && b1==0x8b;
}

static std::unique_ptr<Reader> make_reader(const std::string& path) {
    if (file_is_gzip(path)) return std::make_unique<GzReader>(path);
    return std::make_unique<FileReader>(path);
}

// Varint (LEB128 7-bit groups, like Python read_varint)
static bool read_varint(Reader& r, uint64_t& out) {
    out = 0;
    int shift = 0;
    for (;;) {
        uint8_t b=0;
        if (!r.read1(b)) return false;
        out |= (uint64_t)(b & 0x7F) << shift;
        if (!(b & 0x80)) return true;
        shift += 7;
        if (shift > 63) return false;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Build CIGAR from mapping edits (same semantics as Python)
static std::string build_cigar(const google::protobuf::RepeatedPtrField<vg::Edit>& edits) {
    std::string out;
    out.reserve(64);
    for (const auto& e : edits) {
        int fL = e.from_length();
        int tL = e.to_length();
        int sL = (int)e.sequence().size();
        if (fL == tL) {
            // match or mismatch block: use M if no sequence; X if sequence present
            out += std::to_string(fL);
            out += (sL == 0 ? 'M' : 'X');
        } else if (fL > tL) {
            out += std::to_string(fL - tL);
            out += 'D';
        } else if (fL < tL) {
            out += std::to_string(tL - fL);
            out += 'I';
        } else {
            throw std::runtime_error("Unexpected edit");
        }
    }
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Alignment → per-node segments (packing later)
static void process_alignment(
    const std::string& raw,
    const std::unordered_set<uint32_t>& wanted_nodes,
    const std::string& chrom_filter,
    std::unordered_map<uint32_t, std::vector<Segment>>& out_map)
{
    vg::Alignment aln;
    if (!aln.ParseFromArray(raw.data(), (int)raw.size())) return;

    if (aln.mapping_quality() <= 10) return;

    if (!chrom_filter.empty()) {
        bool ok = false;
        for (const auto& pos : aln.refpos()) {
            if (pos.name() == chrom_filter) { ok = true; break; }
        }
        if (!ok) return;
    }

    const std::string& read_seq = aln.sequence();
    const std::string& read_qual = aln.quality();
    int16_t mapq = (int16_t)aln.mapping_quality();
    size_t read_offset = 0;

    for (const auto& mapping : aln.path().mapping()) {
        uint32_t nid = (uint32_t)mapping.position().node_id();

        // advance read_offset even if not wanted
        if (!wanted_nodes.count(nid)) {
            for (const auto& e : mapping.edit()) {
                read_offset += (size_t)e.to_length();
            }
            continue;
        }

        int32_t node_off = (int32_t)mapping.position().offset();
        char strand = mapping.position().is_reverse() ? '-' : '+';

        std::string seq_parts;
        std::string bq_parts;
        seq_parts.reserve(256);
        bq_parts.reserve(256);
        std::string cig; cig.reserve(64);

        for (const auto& e : mapping.edit()) {
            int fL = e.from_length();
            int tL = e.to_length();
            int sL = (int)e.sequence().size();

            if (fL == tL) {
                cig += std::to_string(fL);
                cig.push_back(sL == 0 ? 'M' : 'X');
            } else if (fL > 0 && tL == 0) {
                cig += std::to_string(fL - tL); cig.push_back('D');
            } else if (fL == 0 && tL > 0) {
                cig += std::to_string(tL - fL); cig.push_back('I');
            } else {
                throw std::runtime_error("Unexpected edit combo");
            }

            if (tL > 0) {
                const size_t begin = read_offset;
                const size_t end   = begin + (size_t)tL;
                if (end <= read_seq.size()) {
                    seq_parts.append(read_seq.data() + begin, (size_t)tL);
                    // upcase in place
                    for (size_t i = seq_parts.size() - tL; i < seq_parts.size(); ++i) {
                        char& c = seq_parts[i];
                        if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
                    }
                }
                if (end <= read_qual.size()) {
                    bq_parts.append(read_qual.data() + begin, (size_t)tL);
                }
                read_offset = end;
            }
        }

        Segment s;
        s.offset = node_off;
        s.seq    = std::move(seq_parts);
        s.bq     = std::move(bq_parts);
        s.cigar  = std::move(cig);
        s.rq     = mapq;
        s.strand = strand;

        out_map[nid].push_back(std::move(s));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GAM record iterator over the custom framed groups (tag "GAM")
class GamIterator {
    std::unique_ptr<Reader> rd_;
    bool good_ = true;
public:
    explicit GamIterator(const std::string& path) : rd_(make_reader(path)) {}

    // Fills 'out' with the next raw message payload; returns false on EOF
    bool next(std::string& out) {
        if (!good_) return false;

        for (;;) {
            uint64_t group_count = 0;
            if (!read_varint(*rd_, group_count)) { good_ = false; return false; }
            if (group_count == 0) continue;

            uint64_t tag_len = 0;
            if (!read_varint(*rd_, tag_len)) { good_ = false; return false; }
            std::string tag(tag_len, '\0');
            if (!rd_->readN(tag.data(), tag_len)) { good_ = false; return false; }

            // Only process "GAM"
            const bool want = (tag == "GAM");
            // There are (group_count - 1) messages
            for (uint64_t i = 0; i + 1 < group_count; ++i) {
                uint64_t msg_size = 0;
                if (!read_varint(*rd_, msg_size)) { good_ = false; return false; }
                if (!want) {
                    if (!rd_->skipN((size_t)msg_size)) { good_ = false; return false; }
                    continue;
                }
                out.resize((size_t)msg_size);
                if (!rd_->readN(out.data(), (size_t)msg_size)) { good_ = false; return false; }
                return true; // one message per call
            }
            // if tag not "GAM", loop again to read next group
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// IO helpers for .dat / .idx (initialization & reuse)
static void write_exact(int fd, const void* buf, size_t n) {
    const uint8_t* p = static_cast<const uint8_t*>(buf);
    size_t sent = 0;
    while (sent < n) {
        ssize_t w = ::write(fd, p + sent, n - sent);
        if (w <= 0) {
            perror("write");
            throw std::runtime_error("write_exact failed");
        }
        sent += size_t(w);
    }
}

static void pwrite_exact(int fd, const void* buf, size_t n, off_t off) {
    const uint8_t* p = static_cast<const uint8_t*>(buf);
    size_t sent = 0;
    while (sent < n) {
        ssize_t w = ::pwrite(fd, p + sent, n - sent, off + sent);
        if (w <= 0) {
            perror("pwrite");
            throw std::runtime_error("pwrite_exact failed");
        }
        sent += size_t(w);
    }
}

struct InitResult {
    std::unordered_map<uint32_t, BlockInfo> block_infos;
    std::string dat_path;
    std::unordered_set<uint32_t> wanted_nodes;
};

static size_t record_size(int R, int C) {
    return RecordLayout(R, C).size;
}

// From JSON stats: keep nodes where (perfect+not_perfect)>0 AND not_perfect>1 AND not_perfect/total>0.05
static InitResult initialize_from_stats_json(const std::string& stats_json_path,
                                             const std::string& out_prefix)
{
    std::ifstream in(stats_json_path);
    if (!in) throw std::runtime_error("Cannot open stats_json: " + stats_json_path);
    json j; in >> j;

    std::unordered_set<uint32_t> wanted;
    std::unordered_map<uint32_t, uint32_t> node_counts;
    std::unordered_map<uint32_t, std::pair<int,int>> maxima;

    for (auto it = j.begin(); it != j.end(); ++it) {
        uint32_t nid = (uint32_t)std::stoul(it.key());
        const json& stat = it.value();
        int perfect = stat.value("perfect", 0);
        int not_perfect = stat.value("not_perfect", 0);
        int total = perfect + not_perfect;
        if (total > 0 && not_perfect > 1 && (double)not_perfect / (double)total > 0.05) {
            wanted.insert(nid);
            node_counts[nid] = (uint32_t)total;
            int R = std::max(1, stat.value("max_read_length", 1));
            int C = std::max(1, stat.value("max_cigar_length", 1));
            maxima[nid] = {R, C};
        }
    }

    std::cerr << "Filtered " << wanted.size() << " nodes from " << j.size() << " total.\n";

    // Build block table and pre-allocate .dat
    std::map<uint32_t, BlockInfo> sorted; // deterministic order
    uint64_t current_off = GLOBAL_HEADER_SIZE;
    for (uint32_t nid : wanted) {
        // Ensure deterministic order by moving into std::map
        sorted[nid] = {};
    }
    for (auto& kv : sorted) {
        uint32_t nid = kv.first;
        uint32_t nrec = node_counts[nid];
        auto [R, C] = maxima[nid];
        uint32_t rec_sz = (uint32_t)record_size(R, C);
        uint32_t blk_sz = (uint32_t)(BLOCK_HDR_SIZE + (uint64_t)nrec * rec_sz);

        BlockInfo bi{};
        bi.offset      = current_off;
        bi.n_records   = nrec;
        bi.current_pos = 0;
        bi.R           = (uint32_t)R;
        bi.C           = (uint32_t)C;
        bi.record_size = rec_sz;
        bi.block_size  = blk_sz;

        kv.second = bi;
        current_off += blk_sz;
    }

    const std::string dat_path = out_prefix + ".dat";
    int fd = ::open(dat_path.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0644);
    if (fd < 0) { perror("open dat"); throw std::runtime_error("open .dat failed"); }

    // Global magic + version
    write_exact(fd, GLOBAL_MAGIC, 6);
    GlobalVerPack ver{GLOBAL_MAJOR, GLOBAL_MINOR, (uint32_t)sorted.size(), {0}};
    write_exact(fd, &ver, sizeof(ver));

    // Block headers
    for (auto& kv : sorted) {
        BlockHdr hdr{};
        hdr.node_id       = kv.first;
        hdr.n_records     = kv.second.n_records;
        hdr.flags         = 0;
        hdr.max_read_len  = kv.second.R;
        hdr.max_cigar_len = kv.second.C;
        write_exact(fd, &hdr, sizeof(hdr));
    }

    // Pre-allocate file: seek to end-1 and write 0
    if (current_off > GLOBAL_HEADER_SIZE) {
        if (lseek(fd, (off_t)current_off - 1, SEEK_SET) == (off_t)-1) {
            perror("lseek");
            throw std::runtime_error("preallocate lseek failed");
        }
        uint8_t zero = 0;
        write_exact(fd, &zero, 1);
    }
    ::close(fd);

    // Write .idx
    const std::string idx_path = out_prefix + ".idx";
    std::ofstream idx(idx_path, std::ios::binary | std::ios::trunc);
    if (!idx) throw std::runtime_error("open .idx failed");

    // first: uint32_t count
    uint32_t count = (uint32_t)sorted.size();
    idx.write(reinterpret_cast<const char*>(&count), sizeof(count));
    // entries: "<I Q I I H I I" = 4+8+4+4+2+4+4 = 30 bytes
    for (const auto& kv : sorted) {
        const auto& bi = kv.second;
        uint32_t nid = kv.first;
        uint64_t off = bi.offset;
        uint32_t blk_sz = bi.block_size;
        uint32_t nrec = bi.n_records;
        uint16_t flags = 0;
        uint32_t R = bi.R, C = bi.C;

        idx.write(reinterpret_cast<const char*>(&nid), 4);
        idx.write(reinterpret_cast<const char*>(&off), 8);
        idx.write(reinterpret_cast<const char*>(&blk_sz), 4);
        idx.write(reinterpret_cast<const char*>(&nrec), 4);
        idx.write(reinterpret_cast<const char*>(&flags), 2);
        idx.write(reinterpret_cast<const char*>(&R), 4);
        idx.write(reinterpret_cast<const char*>(&C), 4);
    }

    InitResult res;
    res.dat_path = dat_path;
    res.wanted_nodes = std::unordered_set<uint32_t>(wanted.begin(), wanted.end());
    // move sorted map to unordered_map
    for (auto& kv : sorted) res.block_infos.emplace(kv.first, kv.second);
    return res;
}

// Reuse existing .dat/.idx (verifies block headers)
static InitResult reuse_existing(const std::string& out_prefix) {
    const std::string dat_path = out_prefix + ".dat";
    const std::string idx_path = out_prefix + ".idx";

    // Read idx
    std::ifstream idx(idx_path, std::ios::binary);
    if (!idx) throw std::runtime_error("Missing .idx: " + idx_path);
    uint32_t count=0;
    idx.read(reinterpret_cast<char*>(&count), 4);
    struct IdxEntry { uint32_t nid; uint64_t off; uint32_t blk_sz; uint32_t nrec; uint16_t flags; uint32_t R; uint32_t C; };
    std::vector<IdxEntry> ent(count);
    for (uint32_t i=0;i<count;++i) {
        idx.read(reinterpret_cast<char*>(&ent[i].nid), 4);
        idx.read(reinterpret_cast<char*>(&ent[i].off), 8);
        idx.read(reinterpret_cast<char*>(&ent[i].blk_sz), 4);
        idx.read(reinterpret_cast<char*>(&ent[i].nrec), 4);
        idx.read(reinterpret_cast<char*>(&ent[i].flags), 2);
        idx.read(reinterpret_cast<char*>(&ent[i].R), 4);
        idx.read(reinterpret_cast<char*>(&ent[i].C), 4);
    }

    // Verify dat header & block headers
    int fd = ::open(dat_path.c_str(), O_RDONLY);
    if (fd < 0) { perror("open dat"); throw std::runtime_error("Missing .dat: " + dat_path); }

    char magic[6];
    if (::read(fd, magic, 6) != 6 || std::memcmp(magic, GLOBAL_MAGIC, 6)!=0) {
        ::close(fd);
        throw std::runtime_error("Invalid .dat magic");
    }
    GlobalVerPack ver{};
    if (::read(fd, &ver, sizeof(ver)) != (ssize_t)sizeof(ver)) {
        ::close(fd);
        throw std::runtime_error("Invalid .dat version block");
    }
    if (ver.block_count != count) {
        std::cerr << "[warn] .dat block_count (" << ver.block_count << ") != .idx (" << count << ")\n";
    }

    std::unordered_map<uint32_t, BlockInfo> block_infos;
    std::unordered_set<uint32_t> wanted;

    for (const auto& e : ent) {
        // seek to block header in dat and read it back
        if (lseek(fd, (off_t)e.off, SEEK_SET) == (off_t)-1) {
            ::close(fd); throw std::runtime_error("lseek .dat failed");
        }
        BlockHdr hdr{};
        if (::read(fd, &hdr, sizeof(hdr)) != (ssize_t)sizeof(hdr)) {
            ::close(fd); throw std::runtime_error("short read block hdr");
        }
        if (hdr.node_id != e.nid || hdr.n_records != e.nrec) {
            std::cerr << "[warn] .dat/.idx mismatch for node " << e.nid << "\n";
        }
        uint32_t R = hdr.max_read_len ? hdr.max_read_len : e.R;
        uint32_t C = hdr.max_cigar_len ? hdr.max_cigar_len : e.C;
        uint32_t rec_sz = (uint32_t)record_size(R, C);

        BlockInfo bi{};
        bi.offset      = e.off;
        bi.n_records   = e.nrec;
        bi.current_pos = 0;
        bi.R           = R;
        bi.C           = C;
        bi.record_size = rec_sz;
        bi.block_size  = e.blk_sz;

        block_infos.emplace(e.nid, bi);
        wanted.insert(e.nid);
    }
    ::close(fd);

    InitResult res;
    res.dat_path = dat_path;
    res.block_infos = std::move(block_infos);
    res.wanted_nodes = std::move(wanted);
    std::cerr << "Reusing existing output with " << res.block_infos.size() << " node blocks from " << idx_path << "\n";
    return res;
}

// ─────────────────────────────────────────────────────────────────────────────
// Buffer flushing (packs per node then single pwrite)
static void flush_segment_buffer(
    int dat_fd,
    std::unordered_map<uint32_t, BlockInfo>& block_infos,
    std::unordered_map<uint32_t, std::vector<Segment>>& segbuf)
{
    if (segbuf.empty()) return;

    while (!segbuf.empty()) {
        auto it = segbuf.begin();
        uint32_t nid = it->first;
        auto segs = std::move(it->second);
        segbuf.erase(it);
        if (segs.empty()) continue;

        auto bi_it = block_infos.find(nid);
        if (bi_it == block_infos.end()) continue;
        BlockInfo& bi = bi_it->second;

        const uint64_t base_offset = bi.offset + BLOCK_HDR_SIZE;
        const uint32_t R = bi.R, C = bi.C;
        const size_t rec_sz = bi.record_size;

        const size_t n = segs.size();
        std::vector<uint8_t> buf(n * rec_sz);

        // pack
        uint8_t* p = buf.data();
        for (const auto& s : segs) {
            // <h
            int16_t off = (int16_t)s.offset;
            std::memcpy(p, &off, sizeof(off)); p += sizeof(off);

            // seq[R]
            if (!s.seq.empty()) {
                size_t copy = std::min((size_t)R, s.seq.size());
                std::memcpy(p, s.seq.data(), copy);
                if (copy < (size_t)R) std::memset(p + copy, 0, R - copy);
            } else {
                std::memset(p, 0, R);
            }
            p += R;

            // bq[R]
            if (!s.bq.empty()) {
                size_t copy = std::min((size_t)R, s.bq.size());
                std::memcpy(p, s.bq.data(), copy);
                if (copy < (size_t)R) std::memset(p + copy, 0, R - copy);
            } else {
                std::memset(p, 0, R);
            }
            p += R;

            // cigar[C]
            if (!s.cigar.empty()) {
                size_t copy = std::min((size_t)C, s.cigar.size());
                std::memcpy(p, s.cigar.data(), copy);
                if (copy < (size_t)C) std::memset(p + copy, 0, C - copy);
            } else {
                std::memset(p, 0, C);
            }
            p += C;

            // rq <h
            int16_t rq = s.rq;
            std::memcpy(p, &rq, sizeof(rq)); p += sizeof(rq);

            // strand c
            char st = (s.strand == '+' || s.strand == '-') ? s.strand : '+';
            std::memcpy(p, &st, sizeof(st)); p += sizeof(st);
        }

        // write
        uint64_t pos = base_offset + (uint64_t)bi.current_pos * rec_sz;
        pwrite_exact(dat_fd, buf.data(), buf.size(), (off_t)pos);
        bi.current_pos += (uint32_t)n;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CLI + main pipeline
struct Args {
    std::string gam_path;
    std::string stats_json;
    std::string out_prefix;
    uint64_t milestone = 1'000'000;
    std::string chr_filter; // optional
    bool use_existing = false;
};

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i=1;i<argc;++i) {
        std::string s = argv[i];
        auto need = [&](const char* flag){
            if (i+1>=argc) throw std::runtime_error(std::string("Missing value for ")+flag);
            return std::string(argv[++i]);
        };
        if (s=="--gam") a.gam_path = need("--gam");
        else if (s=="--stats-json") a.stats_json = need("--stats-json");
        else if (s=="--out-prefix") a.out_prefix = need("--out-prefix");
        else if (s=="--milestone") a.milestone = std::stoull(need("--milestone"));
        else if (s=="--chr") a.chr_filter = need("--chr");
        else if (s=="--use-existing") a.use_existing = true;
        else if (s=="-h" || s=="--help") {
            std::cout <<
            "Usage: gam_segments --gam file.gam[.gz] --stats-json stats.json --out-prefix out\n"
            "       [--milestone 1000000] [--chr chr1] [--use-existing]\n";
            std::exit(0);
        }
        else {
            std::cerr << "Unknown arg: " << s << "\n";
            std::exit(2);
        }
    }
    if (!a.use_existing) {
        if (a.gam_path.empty() || a.stats_json.empty() || a.out_prefix.empty())
            throw std::runtime_error("Missing required args");
    } else {
        if (a.out_prefix.empty())
            throw std::runtime_error("--use-existing needs --out-prefix");
        if (a.gam_path.empty())
            throw std::runtime_error("--use-existing needs --gam");
    }
    return a;
}

int main(int argc, char** argv) {
    google::protobuf::ShutdownProtobufLibrary(); // ensure symbol linkage
    try {
        Args args = parse_args(argc, argv);

        InitResult init;
        if (args.use_existing) {
            std::cerr << "Reusing existing .dat/.idx...\n";
            init = reuse_existing(args.out_prefix);
        } else {
            std::cerr << "Initializing output files from JSON maxima...\n";
            init = initialize_from_stats_json(args.stats_json, args.out_prefix);
            std::cerr << "Output file created: " << init.dat_path << "\n";
        }

        // Open .dat for writing
        int dat_fd = ::open(init.dat_path.c_str(), O_RDWR);
        if (dat_fd < 0) { perror("open .dat"); return 1; }

        GamIterator it(args.gam_path);

        const uint64_t BUFFER_SEGMENTS = 500'000'000ull;
        uint64_t next_milestone = args.milestone;
        uint64_t total_reads = 0;
        uint64_t total_segments = 0;

        auto t0 = std::chrono::steady_clock::now();

        std::unordered_map<uint32_t, std::vector<Segment>> segment_buffer;
        segment_buffer.reserve(1 << 16);

        std::string raw;
        while (it.next(raw)) {
            std::unordered_map<uint32_t, std::vector<Segment>> segs_by_node;
            process_alignment(raw, init.wanted_nodes, args.chr_filter, segs_by_node);
            ++total_reads;

            for (auto& kv : segs_by_node) {
                auto& dst = segment_buffer[kv.first];
                auto& src = kv.second;
                total_segments += src.size();
                // append
                dst.insert(dst.end(),
                           std::make_move_iterator(src.begin()),
                           std::make_move_iterator(src.end()));
            }

            if (total_segments >= BUFFER_SEGMENTS) {
                flush_segment_buffer(dat_fd, init.block_infos, segment_buffer);
                total_segments = 0;
            }

            if (args.milestone && total_reads >= next_milestone) {
                auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                std::cerr << total_reads << " reads processed | " << dt << " s\n";
                next_milestone += args.milestone;
            }
        }

        // final flush
        flush_segment_buffer(dat_fd, init.block_infos, segment_buffer);
        ::close(dat_fd);

        auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "\nFinal Summary:\n";
        std::cout << "  Total reads processed: " << total_reads << "\n";
        std::cout << "  Nodes included: " << init.block_infos.size() << "\n";
        std::cout << "  Elapsed time: " << dt << " seconds\n";

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[error] " << ex.what() << "\n";
        return 1;
    }
}
