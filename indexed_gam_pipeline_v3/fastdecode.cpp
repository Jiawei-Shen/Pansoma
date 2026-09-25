// Native decoder: candidates.decode_alignment (with left_align_indels and indel_observations)
// in C++, reading the serialized vg Alignment directly. A line-by-line port of the Python
// reference in candidates.py -- it changes no rule; native.py wraps the result into the same
// Read / Observation / unsupported-event objects, checks it against the Python decoder at load
// time and falls back to Python for any record this code fails on.
//
// Portable by construction: standard C++17 and pybind11 only, no CPU-specific instructions,
// no floating-point shortcuts (the only float, an insertion's mean quality, is one exact IEEE
// division). Built by `python -m <package>.native compile`, which compiles in the
// SHA-256 of this file (FASTDECODE_SOURCE_SHA256) so a stale build is never used.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

#ifndef FASTDECODE_SOURCE_SHA256
#define FASTDECODE_SOURCE_SHA256 "unknown"
#endif

namespace {

struct Col {  // one alignment column; chars are ASCII codes, '-' for a gap
    uint8_t read, ref, op;
    uint8_t reverse, boundary;
    int16_t quality;
    int32_t pos;
    int32_t visit;
    int64_t node;
};

struct VisitRec {
    int64_t node;
    int64_t start, end;
    uint8_t reverse;
    int32_t index, first, last;
    int64_t node_length;
};

// ---------------------------------------------------------------- protobuf wire format

struct EditRaw { int64_t f = 0, t = 0; std::string seq; };
struct MappingRaw { int64_t node = 0, offset = 0; bool reverse = false; std::vector<EditRaw> edits; };
struct AlignmentRaw { std::string sequence, quality; std::vector<MappingRaw> mappings; int64_t mapq = 0; };

static uint64_t varint(const uint8_t*& p, const uint8_t* end) {
    uint64_t result = 0;
    for (int shift = 0; shift < 64; shift += 7) {
        if (p >= end) throw py::value_error("truncated varint");
        uint8_t b = *p++;
        result |= uint64_t(b & 0x7f) << shift;
        if (!(b & 0x80)) return result;
    }
    throw py::value_error("varint too long");
}

static void skip(const uint8_t*& p, const uint8_t* end, int wire) {
    switch (wire) {
        case 0: varint(p, end); break;
        case 1: p += 8; break;
        case 2: { uint64_t n = varint(p, end); p += n; break; }
        case 5: p += 4; break;
        default: throw py::value_error("unsupported wire type");
    }
    if (p > end) throw py::value_error("truncated field");
}

static std::string_view bytes_field(const uint8_t*& p, const uint8_t* end) {
    uint64_t n = varint(p, end);
    if (n > uint64_t(end - p)) throw py::value_error("truncated length-delimited field");
    std::string_view v(reinterpret_cast<const char*>(p), n);
    p += n;
    return v;
}

static int32_t as_int32(uint64_t v) { return static_cast<int32_t>(static_cast<uint32_t>(v)); }

static void parse_position(const uint8_t* p, const uint8_t* end, MappingRaw& m) {
    while (p < end) {
        uint64_t key = varint(p, end);
        int field = int(key >> 3), wire = int(key & 7);
        if (field == 1 && wire == 0) m.node = static_cast<int64_t>(varint(p, end));
        else if (field == 2 && wire == 0) m.offset = static_cast<int64_t>(varint(p, end));
        else if (field == 4 && wire == 0) m.reverse = varint(p, end) != 0;
        else skip(p, end, wire);
    }
}

static void parse_edit(const uint8_t* p, const uint8_t* end, EditRaw& e) {
    while (p < end) {
        uint64_t key = varint(p, end);
        int field = int(key >> 3), wire = int(key & 7);
        if (field == 1 && wire == 0) e.f = as_int32(varint(p, end));
        else if (field == 2 && wire == 0) e.t = as_int32(varint(p, end));
        else if (field == 3 && wire == 2) e.seq = std::string(bytes_field(p, end));
        else skip(p, end, wire);
    }
}

static void parse_mapping(const uint8_t* p, const uint8_t* end, MappingRaw& m) {
    while (p < end) {
        uint64_t key = varint(p, end);
        int field = int(key >> 3), wire = int(key & 7);
        if (field == 1 && wire == 2) {
            auto v = bytes_field(p, end);
            auto* q = reinterpret_cast<const uint8_t*>(v.data());
            parse_position(q, q + v.size(), m);
        } else if (field == 2 && wire == 2) {
            auto v = bytes_field(p, end);
            auto* q = reinterpret_cast<const uint8_t*>(v.data());
            m.edits.emplace_back();
            parse_edit(q, q + v.size(), m.edits.back());
        } else skip(p, end, wire);
    }
}

static void parse_path(const uint8_t* p, const uint8_t* end, AlignmentRaw& a) {
    while (p < end) {
        uint64_t key = varint(p, end);
        int field = int(key >> 3), wire = int(key & 7);
        if (field == 2 && wire == 2) {
            auto v = bytes_field(p, end);
            auto* q = reinterpret_cast<const uint8_t*>(v.data());
            a.mappings.emplace_back();
            parse_mapping(q, q + v.size(), a.mappings.back());
        } else skip(p, end, wire);
    }
}

static AlignmentRaw parse_alignment(std::string_view raw) {
    AlignmentRaw a;
    auto* p = reinterpret_cast<const uint8_t*>(raw.data());
    auto* end = p + raw.size();
    while (p < end) {
        uint64_t key = varint(p, end);
        int field = int(key >> 3), wire = int(key & 7);
        if (field == 1 && wire == 2) a.sequence = std::string(bytes_field(p, end));
        else if (field == 4 && wire == 2) a.quality = std::string(bytes_field(p, end));
        else if (field == 5 && wire == 0) a.mapq = as_int32(varint(p, end));
        else if (field == 2 && wire == 2) {
            auto v = bytes_field(p, end);
            auto* q = reinterpret_cast<const uint8_t*>(v.data());
            parse_path(q, q + v.size(), a);
        } else skip(p, end, wire);
    }
    return a;
}

// ---------------------------------------------------------------- helpers (Python semantics)

static inline char comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; case 'N': return 'N';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a'; case 'n': return 'n';
        default: return c;
    }
}
static std::string rc(std::string_view s) {
    std::string r(s.rbegin(), s.rend());
    for (auto& c : r) c = comp(c);
    return r;
}
static inline char up(char c) { return (c >= 'a' && c <= 'z') ? char(c - 32) : c; }
static std::string upper(std::string_view s) {
    std::string r(s);
    for (auto& c : r) c = up(c);
    return r;
}
static inline bool acgt(uint8_t c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

struct Run { int64_t s, e; uint8_t op; };

// candidates.indel_runs
static std::vector<Run> indel_runs(const std::vector<Col>& cols, bool join) {
    std::vector<Run> runs;
    int64_t n = cols.size(), i = 0;
    while (i < n) {
        const Col& c = cols[i];
        if (c.op == 'I' || c.op == 'D') {
            int64_t j = i + 1;
            while (j < n && cols[j].op == c.op && (cols[j].visit == c.visit || (join && cols[j].reverse == c.reverse))) ++j;
            runs.push_back({i, j, c.op});
            i = j;
        } else ++i;
    }
    return runs;
}

// candidates._shiftable
static inline bool shiftable(const Col& nb, bool reverse) {
    return nb.op == 'M' && !nb.boundary && bool(nb.reverse) == reverse && nb.read == nb.ref && acgt(nb.read);
}

static inline bool has_n(const std::vector<Col>& cols, int64_t lo, int64_t hi, uint8_t op) {
    for (int64_t k = lo; k < hi; ++k) if (up(char(op == 'I' ? cols[k].read : cols[k].ref)) == 'N') return true;
    return false;
}

// candidates.left_align_indels; returns moves, updates visit column ranges
static std::vector<std::pair<int64_t, int64_t>> left_align(std::vector<Col>& cols, std::vector<VisitRec>& visits,
                                                           int64_t max_indel) {
    int64_t n = cols.size();
    std::vector<Run> fwd, rev;
    for (const Run& r : indel_runs(cols, true)) {
        if (r.e - r.s > max_indel || has_n(cols, r.s, r.e, r.op)) continue;
        (cols[r.s].reverse ? rev : fwd).push_back(r);
    }
    std::stable_sort(fwd.begin(), fwd.end(), [](const Run& a, const Run& b) { return a.s < b.s; });
    std::stable_sort(rev.begin(), rev.end(), [](const Run& a, const Run& b) { return a.s > b.s; });
    std::vector<Run> order(fwd);
    order.insert(order.end(), rev.begin(), rev.end());
    std::vector<std::pair<int64_t, int64_t>> moves;
    std::vector<Col> tmp;
    for (Run run : order) {
        int64_t s = run.s, e = run.e;
        uint8_t op = run.op;
        bool r = cols[s].reverse;
        int64_t origin = cols[s].node;
        int step = r ? 1 : -1;
        while (true) {
            int64_t k = r ? e : s - 1;
            if (!(0 <= k && k < n)) break;
            Col nb = cols[k];
            if (nb.op == op && bool(nb.reverse) == r) {
                int64_t j = k;
                while (0 <= j + step && j + step < n) {
                    const Col& c = cols[j + step];
                    if (c.op != op || bool(c.reverse) != r) break;
                    j += step;
                }
                int64_t lo = r ? s : j, hi = r ? j + 1 : e;
                if (hi - lo > max_indel || has_n(cols, lo, hi, op)) break;
                s = lo; e = hi;
                continue;
            }
            if (!shiftable(nb, r)) break;
            if (op == 'I') {
                Col edge = r ? cols[s] : cols[e - 1];
                if (nb.read != edge.read) break;
                tmp.clear();
                auto inserted = [&](uint8_t b, int16_t q) {
                    Col c{b, '-', 'I', nb.reverse, 1, q, nb.pos, nb.visit, nb.node};
                    tmp.push_back(c);
                };
                Col matched{edge.read, nb.ref, 'M', nb.reverse, 0, edge.quality, nb.pos, nb.visit, nb.node};
                if (r) {
                    tmp.push_back(matched);
                    for (int64_t x = s + 1; x < e; ++x) inserted(cols[x].read, cols[x].quality);
                    inserted(nb.read, nb.quality);
                    std::copy(tmp.begin(), tmp.end(), cols.begin() + s);  // cols[s:e+1]
                    s += 1; e += 1;
                } else {
                    inserted(nb.read, nb.quality);
                    for (int64_t x = s; x < e - 1; ++x) inserted(cols[x].read, cols[x].quality);
                    tmp.push_back(matched);
                    std::copy(tmp.begin(), tmp.end(), cols.begin() + (s - 1));  // cols[s-1:e]
                    s -= 1; e -= 1;
                }
            } else {
                Col edge = r ? cols[s] : cols[e - 1];
                if (nb.ref != edge.ref) break;
                Col gap{'-', nb.ref, 'D', nb.reverse, 0, -1, nb.pos, nb.visit, nb.node};
                Col base{nb.read, edge.ref, 'M', edge.reverse, 0, nb.quality, edge.pos, edge.visit, edge.node};
                cols[k] = gap;
                if (r) { cols[s] = base; s += 1; e += 1; }
                else { cols[e - 1] = base; s -= 1; e -= 1; }
            }
        }
        if (op == 'I') {
            int64_t k = r ? s - 1 : e;
            Col owner = (0 <= k && k < n && !cols[k].boundary && bool(cols[k].reverse) == r) ? cols[k]
                                                                                            : (r ? cols[e - 1] : cols[s]);
            for (int64_t x = s; x < e; ++x) {
                cols[x].node = owner.node; cols[x].pos = owner.pos; cols[x].visit = owner.visit; cols[x].boundary = 1;
            }
        }
        if (cols[s].node != origin) moves.emplace_back(origin, cols[s].node);
    }
    std::vector<int64_t> first(visits.size(), -1), last(visits.size(), -1);
    for (int64_t i = 0; i < n; ++i) {
        int32_t v = cols[i].visit;
        if (first[v] < 0) first[v] = i;
        last[v] = i + 1;
    }
    int64_t edge = 0;
    for (auto& visit : visits) {
        if (first[visit.index] >= 0) { visit.first = int32_t(first[visit.index]); visit.last = int32_t(last[visit.index]); }
        else { visit.first = int32_t(edge); visit.last = int32_t(edge); }
        edge = visit.last;
    }
    return moves;
}

// ---------------------------------------------------------------- the decoder

struct Decoder {
    py::object sequences;  // node -> forward sequence (str); any mapping
    py::object targets;    // container supporting `in`, or None
    // per read: node -> (owned str, its UTF-8 view)
    std::unordered_map<int64_t, std::pair<py::object, std::string_view>> seq_cache;

    std::string_view seq(int64_t node) {
        auto it = seq_cache.find(node);
        if (it != seq_cache.end()) return it->second.second;
        py::int_ key(node);
        PyObject* value = PyObject_GetItem(sequences.ptr(), key.ptr());  // KeyError like the reference
        if (value == nullptr) throw py::error_already_set();
        py::object owned = py::reinterpret_steal<py::object>(value);
        Py_ssize_t size;
        const char* data = PyUnicode_AsUTF8AndSize(owned.ptr(), &size);  // cached inside the str object
        if (data == nullptr) throw py::error_already_set();
        std::string_view v(data, size_t(size));
        seq_cache.emplace(node, std::make_pair(owned, v));
        return v;
    }
    bool is_target(int64_t node) {
        if (targets.is_none()) return true;
        py::int_ key(node);
        int hit = PySequence_Contains(targets.ptr(), key.ptr());
        if (hit < 0) throw py::error_already_set();
        return hit == 1;
    }
    char base_at(int64_t node, int64_t pos) {
        auto s = seq(node);
        if (pos < 0 || pos >= int64_t(s.size())) throw py::index_error("string index out of range");
        return s[pos];
    }

    py::tuple decode(py::bytes raw_bytes, int64_t max_indel, bool left) {
        seq_cache.clear();
        char* buffer;
        Py_ssize_t length;
        if (PyBytes_AsStringAndSize(raw_bytes.ptr(), &buffer, &length) != 0) throw py::error_already_set();
        AlignmentRaw a = parse_alignment(std::string_view(buffer, size_t(length)));

        std::vector<Col> cols;
        std::vector<VisitRec> visits;
        py::list observations, unsupported;
        walk(a, max_indel, cols, visits, &observations, &unsupported);
        return finish(cols, visits, max_indel, left, observations, unsupported);
    }

    // The edits of every mapping -> columns and visits (candidates.decode_alignment's walk). SNP
    // observations and unsupported events are collected only when the lists are given.
    void walk(const AlignmentRaw& a, int64_t max_indel, std::vector<Col>& cols, std::vector<VisitRec>& visits,
              py::list* observations, py::list* unsupported) {
        std::string sequence = upper(a.sequence);
        const std::string& qualities = a.quality;
        bool has_q = !qualities.empty();
        if (has_q && qualities.size() != sequence.size()) throw py::value_error("Read quality length differs from sequence length");
        cols.reserve(sequence.size() + 64);
        visits.reserve(a.mappings.size());
        int64_t read_cursor = 0;
        struct E { int64_t f, t; std::string seq; int kind; int index; };
        std::vector<E> edits;
        for (size_t mi = 0; mi < a.mappings.size(); ++mi) {
            const MappingRaw& m = a.mappings[mi];
            int64_t nid = m.node;
            bool observe = observations != nullptr && is_target(nid);
            std::string_view forward = seq(nid);
            bool reverse = m.reverse;
            std::string reference = reverse ? rc(forward) : std::string(forward);
            int64_t cursor = m.offset, initial = cursor;
            int32_t first = int32_t(cols.size());
            edits.clear();
            for (size_t ei = 0; ei < m.edits.size(); ++ei) {
                const EditRaw& o = m.edits[ei];
                int kind = (o.f == 0 && o.t > 0) ? 1 : (o.t == 0 && o.f > 0) ? 2 : 0;
                if (!edits.empty() && kind != 0 && edits.back().kind == kind) {
                    edits.back().f += o.f; edits.back().t += o.t; edits.back().seq += o.seq;
                } else edits.push_back({o.f, o.t, o.seq, kind, int(ei)});
            }
            for (const E& ed : edits) {
                int64_t f = ed.f, t = ed.t;
                if (f < 0 || t < 0 || cursor < 0 || cursor + f > int64_t(reference.size()) || read_cursor + t > int64_t(sequence.size()))
                    throw py::value_error("Invalid edit bounds: node " + std::to_string(nid) + ", mapping " + std::to_string(mi) +
                                          ", edit " + std::to_string(ed.index));
                std::string_view ref(reference.data() + cursor, size_t(f));
                std::string_view alt(sequence.data() + read_cursor, size_t(t));
                std::string replacement = upper(ed.seq);
                if (!replacement.empty() && (int64_t(replacement.size()) != t || replacement != alt))
                    throw py::value_error("Edit replacement disagrees with read at node " + std::to_string(nid));
                if (!f && !t) continue;
                char op = (f == t && replacement.empty()) ? 'M' : (f == t) ? 'X' : (!f) ? 'I' : (!t) ? 'D' : 'C';
                int64_t pos = reverse ? int64_t(forward.size()) - cursor - f : cursor;
                auto q_at = [&](int64_t k) -> int { return has_q ? int(uint8_t(qualities[read_cursor + k])) : -1; };
                if (op == 'X') {
                    if (observe) {  // SNP observations only on target nodes
                        std::string cref = reverse ? rc(ref) : std::string(ref), calt = reverse ? rc(alt) : std::string(alt);
                        for (int64_t k = 0; k < f; ++k) {
                            char r = cref[k], b = calt[k];
                            if (r != b && up(r) != 'N' && up(b) != 'N')
                                observations->append(py::make_tuple("SNP", nid, pos + k, py::str(&r, 1), py::str(&b, 1),
                                                                    py::tuple(), int(mi), q_at(reverse ? f - 1 - k : k)));
                        }
                    }
                } else if (unsupported != nullptr && (op == 'I' || op == 'D') && std::max(f, t) > max_indel) {
                    std::string cref = reverse ? rc(ref) : std::string(ref), calt = reverse ? rc(alt) : std::string(alt);
                    unsupported->append(py::make_tuple(0, nid, pos, cref, calt, op == 'I' ? "INS" : "DEL", int(mi), ed.index));
                } else if (unsupported != nullptr && op == 'C') {
                    std::string cref = reverse ? rc(ref) : std::string(ref), calt = reverse ? rc(alt) : std::string(alt);
                    unsupported->append(py::make_tuple(1, nid, pos, cref, calt, "", int(mi), ed.index));
                }
                int64_t width = std::max(f, t), flen = int64_t(forward.size());
                for (int64_t k = 0; k < width; ++k) {
                    bool b = k >= f;
                    int64_t p = cursor + std::min(k, f);
                    if (reverse) p = flen - p - (b ? 0 : 1);
                    cols.push_back(Col{uint8_t(k < t ? alt[k] : '-'), uint8_t(k < f ? ref[k] : '-'), uint8_t(op),
                                       uint8_t(reverse), uint8_t(b), int16_t(k < t ? q_at(k) : -1), int32_t(p),
                                       int32_t(mi), nid});
                }
                cursor += f;
                read_cursor += t;
            }
            int64_t flen = int64_t(forward.size());
            int64_t lo = reverse ? flen - cursor : initial, hi = reverse ? flen - initial : cursor;
            visits.push_back(VisitRec{nid, lo, hi, uint8_t(reverse), int32_t(mi), first, int32_t(cols.size()), flen});
        }
        if (read_cursor != int64_t(sequence.size())) throw py::value_error("GAM edits do not consume the complete read sequence");
    }

    py::tuple finish(std::vector<Col>& cols, std::vector<VisitRec>& visits, int64_t max_indel, bool left,
                     py::list& observations, py::list& unsupported) {
        py::list moves;
        if (left)
            for (auto& mv : left_align(cols, visits, max_indel)) moves.append(py::make_tuple(mv.first, mv.second));

        indel_observations(cols, max_indel, observations);

        for (const Run& run : indel_runs(cols, true)) {
            if (run.e - run.s <= max_indel) continue;
            std::unordered_map<int32_t, int64_t> pieces;
            for (int64_t k = run.s; k < run.e; ++k) pieces[cols[k].visit] += 1;
            int64_t most = 0;
            for (auto& kv : pieces) most = std::max(most, kv.second);
            if (pieces.size() > 1 && most <= max_indel) {
                const Col& first = cols[cols[run.s].reverse ? run.e - 1 : run.s];
                unsupported.append(py::make_tuple(2, first.node, first.pos, run.op == 'I' ? "INS" : "DEL", run.e - run.s,
                                                  first.visit, int64_t(pieces.size())));
            }
        }

        py::array_t<Col> col_array(cols.size());
        if (!cols.empty()) std::memcpy(col_array.mutable_data(), cols.data(), cols.size() * sizeof(Col));
        py::array_t<VisitRec> visit_array(visits.size());
        if (!visits.empty()) std::memcpy(visit_array.mutable_data(), visits.data(), visits.size() * sizeof(VisitRec));
        return py::make_tuple(col_array, visit_array, observations, moves, unsupported);
    }

    // candidates._read_base_quality: -2 for None
    static int read_base_quality(const std::vector<Col>& cols, int64_t i, int step) {
        int64_t n = cols.size();
        while (0 <= i && i < n) {
            if (cols[i].read != '-') return cols[i].quality == -1 ? -2 : cols[i].quality;
            i += step;
        }
        return -2;
    }

    // candidates._repeat_extent
    std::vector<int64_t> repeat_extent(const std::vector<Col>& cols, int64_t s, int64_t e, uint8_t op) {
        bool r = cols[s].reverse;
        std::string bases;
        if (r) {
            for (int64_t k = e - 1; k >= s; --k) bases.push_back(op == 'I' ? comp(char(cols[k].read)) : base_at(cols[k].node, cols[k].pos));
        } else {
            for (int64_t k = s; k < e; ++k) bases.push_back(op == 'I' ? char(cols[k].read) : base_at(cols[k].node, cols[k].pos));
        }
        std::vector<int64_t> passed;
        int64_t n = cols.size(), k = r ? s - 1 : e;
        size_t head = 0;  // bases is a rotating window: bases[head:] + appended
        while (0 <= k && k < n && shiftable(cols[k], r)) {
            char base = base_at(cols[k].node, cols[k].pos);
            if (up(base) != up(bases[head])) break;
            bases.push_back(base);
            ++head;
            passed.push_back(k);
            k += r ? -1 : 1;
        }
        return passed;
    }

    // candidates.indel_observations (appends to `out`)
    void indel_observations(const std::vector<Col>& cols, int64_t max_indel, py::list& out) {
        for (const Run& run : indel_runs(cols, true)) {
            int64_t s = run.s, e = run.e;
            bool r = cols[s].reverse;
            const Col& first = cols[r ? e - 1 : s];
            if (e - s > max_indel || !is_target(first.node)) continue;
            std::vector<int64_t> passed = repeat_extent(cols, s, e, run.op);
            std::string_view forward = seq(first.node);
            if (run.op == 'I') {
                std::string inserted;
                for (int64_t k = s; k < e; ++k) inserted.push_back(char(cols[k].read));
                std::string alt = r ? rc(inserted) : inserted;
                int64_t a = std::max<int64_t>(0, int64_t(first.pos) - 1);
                bool anchor_n = a < int64_t(forward.size()) && up(forward[a]) == 'N';
                if (upper(alt).find('N') != std::string::npos || anchor_n) continue;
                int64_t total = 0, count = 0;
                for (int64_t k = s; k < e; ++k) { total += cols[k].quality; ++count; }
                for (int64_t k : passed) { total += cols[k].quality; ++count; }
                double quality = double(total) / double(count);
                out.append(py::make_tuple("INS", first.node, int64_t(first.pos), "", alt, py::tuple(), first.visit, quality));
            } else {
                std::vector<const Col*> ordered;
                if (r) for (int64_t k = e - 1; k >= s; --k) ordered.push_back(&cols[k]);
                else for (int64_t k = s; k < e; ++k) ordered.push_back(&cols[k]);
                std::string ref;
                for (auto* c : ordered) ref.push_back(base_at(c->node, c->pos));
                if (upper(ref).find('N') != std::string::npos) continue;
                std::vector<std::array<int64_t, 3>> segments;
                for (auto* c : ordered) {
                    if (!segments.empty() && segments.back()[0] == c->node && segments.back()[2] == c->pos) segments.back()[2] += 1;
                    else segments.push_back({c->node, int64_t(c->pos), int64_t(c->pos) + 1});
                }
                int64_t right;
                if (r && !passed.empty()) right = *std::min_element(passed.begin(), passed.end()) - 1;
                else if (!passed.empty()) right = *std::max_element(passed.begin(), passed.end()) + 1;
                else right = r ? s - 1 : e;
                std::vector<int> flank;
                int q1 = read_base_quality(cols, r ? e : s - 1, r ? 1 : -1);
                int q2 = read_base_quality(cols, right, r ? -1 : 1);
                if (q1 != -2) flank.push_back(q1);
                if (q2 != -2) flank.push_back(q2);
                if (flank.empty())
                    for (int64_t k : passed) if (cols[k].quality != -1) flank.push_back(cols[k].quality);
                int quality = flank.empty() ? -1 : *std::min_element(flank.begin(), flank.end());
                py::list path;
                for (size_t x = 1; x < segments.size(); ++x)
                    path.append(py::make_tuple(segments[x][0], segments[x][1], segments[x][2]));
                out.append(py::make_tuple("DEL", first.node, segments[0][1], ref, "", py::tuple(path), first.visit, quality));
            }
        }
    }
};

// ---------------------------------------------------------------- discovery counts

// Per node, in order of first appearance: mappings without / with an edit and the longest read
// (run.discover's node_stats). add_raw applies discover's rule to vg's edits as written;
// add_normalized first decodes the record and left-normalizes its indels exactly like the builder,
// so an indel counts on the node the builder will see it on. A record the normalization cannot
// decode (e.g. a node missing from `sequences`) is counted with the raw rule (`fallbacks`).
struct Discovery {
    struct Stat { int64_t perfect = 0, not_perfect = 0, max_len = 0; };
    std::unordered_map<int64_t, size_t> index;
    std::vector<int64_t> nodes;
    std::vector<Stat> stats;
    int64_t alignments = 0, used = 0, fallbacks = 0;

    void count(int64_t node, bool imperfect, int64_t length) {
        auto it = index.find(node);
        size_t k;
        if (it == index.end()) {
            k = nodes.size();
            index.emplace(node, k);
            nodes.push_back(node);
            stats.emplace_back();
        } else k = it->second;
        Stat& s = stats[k];
        (imperfect ? s.not_perfect : s.perfect) += 1;
        s.max_len = std::max(s.max_len, length);
    }
    void count_raw(const AlignmentRaw& a) {
        int64_t length = int64_t(a.sequence.size());
        for (const MappingRaw& m : a.mappings) {
            if (!m.node) continue;
            bool imperfect = false;
            for (const EditRaw& e : m.edits)
                if (e.f != e.t || !e.seq.empty()) { imperfect = true; break; }
            count(m.node, imperfect, length);
        }
    }
    static std::string_view view(py::handle raw) {
        char* data; Py_ssize_t size;
        if (PyBytes_AsStringAndSize(raw.ptr(), &data, &size) != 0) throw py::error_already_set();
        return std::string_view(data, size_t(size));
    }
    void add_raw(py::list messages, int64_t min_mapq) {
        for (py::handle raw : messages) {
            ++alignments;
            AlignmentRaw a = parse_alignment(view(raw));
            if (a.mapq <= min_mapq) continue;
            ++used;
            count_raw(a);
        }
    }
    void add_normalized(py::list messages, py::object sequences, int64_t min_mapq, int64_t max_indel) {
        std::vector<Col> cols;
        std::vector<VisitRec> visits;
        std::vector<char> edited;
        for (py::handle raw : messages) {
            ++alignments;
            AlignmentRaw a = parse_alignment(view(raw));
            if (a.mapq <= min_mapq) continue;
            ++used;
            cols.clear(); visits.clear();
            try {
                Decoder d{sequences, py::none(), {}};
                d.walk(a, max_indel, cols, visits, nullptr, nullptr);
                left_align(cols, visits, max_indel);
            } catch (const std::exception&) {
                ++fallbacks;
                count_raw(a);
                continue;
            }
            edited.assign(visits.size(), 0);
            for (const Col& c : cols)
                if (c.op != 'M') edited[size_t(c.visit)] = 1;
            int64_t length = int64_t(a.sequence.size());
            for (size_t v = 0; v < visits.size(); ++v)
                if (visits[v].node) count(visits[v].node, edited[v] != 0, length);
        }
    }
    py::tuple result() const {
        size_t n = nodes.size();
        py::array_t<int64_t> node(n), perfect(n), not_perfect(n), max_len(n);
        auto a = node.mutable_unchecked<1>(); auto b = perfect.mutable_unchecked<1>();
        auto c = not_perfect.mutable_unchecked<1>(); auto d = max_len.mutable_unchecked<1>();
        for (size_t k = 0; k < n; ++k) {
            a(k) = nodes[k]; b(k) = stats[k].perfect; c(k) = stats[k].not_perfect; d(k) = stats[k].max_len;
        }
        py::dict counters;
        counters["alignments"] = alignments; counters["used"] = used; counters["fallbacks"] = fallbacks;
        return py::make_tuple(node, perfect, not_perfect, max_len, counters);
    }
};

// Distinct node IDs visited by the records of one group (to look up their sequences).
py::array_t<int64_t> group_nodes(py::list messages) {
    std::vector<int64_t> out;
    for (py::handle raw : messages)
        for (const MappingRaw& m : parse_alignment(Discovery::view(raw)).mappings) out.push_back(m.node);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    py::array_t<int64_t> result(out.size());
    if (!out.empty()) std::memcpy(result.mutable_data(), out.data(), out.size() * sizeof(int64_t));
    return result;
}

}  // namespace

PYBIND11_MODULE(_fastdecode, m) {
    PYBIND11_NUMPY_DTYPE(Col, read, ref, op, reverse, boundary, quality, pos, visit, node);
    PYBIND11_NUMPY_DTYPE(VisitRec, node, start, end, reverse, index, first, last, node_length);
    m.attr("source_sha256") = FASTDECODE_SOURCE_SHA256;
    m.def("decode", [](py::bytes raw, py::object sequences, py::object targets, int64_t max_indel, bool left_align) {
        Decoder d{sequences, targets, {}};
        return d.decode(raw, max_indel, left_align);
    }, py::arg("raw"), py::arg("sequences"), py::arg("targets"), py::arg("max_indel") = 50, py::arg("left_align") = true,
       "Decode one serialized vg Alignment: (columns, visits, observations, moves, unsupported).");
    py::class_<Discovery>(m, "Discovery")
        .def(py::init<>())
        .def("add_raw", &Discovery::add_raw, py::arg("messages"), py::arg("min_mapq"))
        .def("add_normalized", &Discovery::add_normalized, py::arg("messages"), py::arg("sequences"),
             py::arg("min_mapq"), py::arg("max_indel") = 50)
        .def("result", &Discovery::result, "(nodes, perfect, not_perfect, max_read_length, counters), first-appearance order");
    m.def("group_nodes", &group_nodes, py::arg("messages"), "Sorted distinct node IDs of a group's records.");
}
