#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/bytes.h>
#include <pybind11/pytypes.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <stdexcept>

#if defined(_WIN32)
  #include <io.h>
  #define fseeko _fseeki64
  #define ftello _ftelli64
#else
  #include <unistd.h>
#endif

namespace py = pybind11;

struct SegView {
    int16_t offset;
    std::string seq;     // raw bytes
    std::string bq;      // raw bytes
    std::string cigar;   // raw bytes
    int16_t rq;
    char strand;
};

// Extract a segment from either a Python object with attributes or a dict.
static SegView get_segment(const py::handle& obj) {
    auto get_bytes = [&](const py::handle& h)->std::string {
        // Accept both bytes and bytearray; coerce to bytes
        if (PyBytes_Check(h.ptr())) {
            return std::string(PyBytes_AsString(h.ptr()), (size_t)PyBytes_Size(h.ptr()));
        } else if (PyByteArray_Check(h.ptr())) {
            return std::string(PyByteArray_AsString(h.ptr()), (size_t)PyByteArray_Size(h.ptr()));
        } else {
            // fallback: try __bytes__()
            py::object b = py::reinterpret_borrow<py::object>(h);
            if (py::hasattr(b, "__bytes__")) {
                py::bytes bb = b.attr("__bytes__")();
                return std::string(bb);
            }
            // or encode() if it's str
            if (PyUnicode_Check(h.ptr())) {
                py::bytes bb = py::reinterpret_borrow<py::object>(h).attr("encode")();
                return std::string(bb);
            }
            throw std::runtime_error("Expected bytes/bytearray/str-like for seq/bq/cigar/strand");
        }
    };

    SegView s{0, {}, {}, {}, 0, '+'};

    bool is_dict = PyDict_Check(obj.ptr());
    auto get_attr = [&](const char* name)->py::object {
        if (is_dict) {
            py::dict d = py::reinterpret_borrow<py::dict>(obj);
            if (!d.contains(name)) throw std::runtime_error(std::string("Missing key: ") + name);
            return py::reinterpret_borrow<py::object>(d[name]);
        } else {
            if (!py::hasattr(obj, name)) throw std::runtime_error(std::string("Missing attribute: ") + name);
            return py::reinterpret_borrow<py::object>(obj.attr(name));
        }
    };

    // offset (int)
    {
        py::object o = get_attr("offset");
        long long v = py::int_(o);
        if (v < INT16_MIN || v > INT16_MAX) throw std::runtime_error("offset out of int16 range");
        s.offset = static_cast<int16_t>(v);
    }

    // seq/bq/cigar (bytes-like)
    s.seq   = get_bytes(get_attr("seq"));
    s.bq    = get_bytes(get_attr("bq"));
    s.cigar = get_bytes(get_attr("cigar"));

    // rq (int)
    {
        py::object o = get_attr("rq");
        long long v = py::int_(o);
        if (v < INT16_MIN || v > INT16_MAX) throw std::runtime_error("rq out of int16 range");
        s.rq = static_cast<int16_t>(v);
    }

    // strand (bytes-like, take first byte if present)
    {
        std::string st = get_bytes(get_attr("strand"));
        s.strand = st.empty() ? '+' : st[0];
        if (s.strand != '+' && s.strand != '-') s.strand = '+';
    }

    return s;
}

// pad/truncate src to exactly L bytes (null-pad on the right)
static inline void pad_trunc_copy(char* dst, const std::string& src, size_t L) {
    size_t n = src.size();
    if (n >= L) {
        std::memcpy(dst, src.data(), L);
    } else {
        std::memcpy(dst, src.data(), n);
        std::memset(dst + n, 0, L - n);
    }
}

// Compute record size for given R, C  (matches Python struct "<h R s R s C s h c")
static inline size_t record_size(int R, int C) {
    // < : little-endian, standard packing
    // h (2) + R + R + C + h (2) + c (1)
    return 2 + (size_t)R + (size_t)R + (size_t)C + 2 + 1;
}

// Pack a vector of segments into a contiguous byte buffer
static std::string pack_segments(int R, int C, const std::vector<SegView>& segs) {
    if (R <= 0 || C <= 0) throw std::runtime_error("R and C must be positive");
    const size_t rec_sz = record_size(R, C);
    const size_t n = segs.size();
    std::string buf;
    buf.resize(rec_sz * n);

    size_t off = 0;
    for (const auto& s : segs) {
        char* base = buf.data() + off;

        // int16 offset (little-endian)
        int16_t off16 = s.offset;
        std::memcpy(base, &off16, sizeof(int16_t));
        base += 2;

        // seq[R], bq[R], cigar[C]
        pad_trunc_copy(base, s.seq,   (size_t)R); base += R;
        pad_trunc_copy(base, s.bq,    (size_t)R); base += R;
        pad_trunc_copy(base, s.cigar, (size_t)C); base += C;

        // int16 rq
        int16_t rq16 = s.rq;
        std::memcpy(base, &rq16, sizeof(int16_t));
        base += 2;

        // char strand
        *base = s.strand;

        off += rec_sz;
    }
    return buf;
}

// Write at absolute position in a file (like pwrite, but portable).
static void write_at_absolute_pos(FILE* f, uint64_t pos, const char* data, size_t nbytes) {
    if (fseeko(f, static_cast<off_t>(pos), SEEK_SET) != 0) {
        throw std::runtime_error("fseeko failed");
    }
    size_t written = std::fwrite(data, 1, nbytes, f);
    if (written != nbytes) {
        throw std::runtime_error("fwrite failed (short write)");
    }
    // do not fflush here; let the caller control fsync if needed
}

// Python: pack_and_write(dat_path, write_pos, R, C, segments)
static size_t pack_and_write_py(const std::string& dat_path,
                                uint64_t write_pos,
                                int R, int C,
                                py::sequence py_segments) {
    // Collect segments
    std::vector<SegView> segs;
    segs.reserve(py_segments.size());
    for (auto item : py_segments) {
        segs.push_back(get_segment(item));
    }
    // Pack
    std::string buf = pack_segments(R, C, segs);

    // Write
    FILE* f = std::fopen(dat_path.c_str(), "r+b");
    if (!f) throw std::runtime_error("Failed to open dat file for r+b: " + dat_path);
    try {
        write_at_absolute_pos(f, write_pos, buf.data(), buf.size());
    } catch (...) {
        std::fclose(f);
        throw;
    }
    std::fclose(f);
    return buf.size();
}

// Python: pack_only(R, C, segments) -> bytes
static py::bytes pack_only_py(int R, int C, py::sequence py_segments) {
    std::vector<SegView> segs;
    segs.reserve(py_segments.size());
    for (auto item : py_segments) {
        segs.push_back(get_segment(item));
    }
    std::string buf = pack_segments(R, C, segs);
    return py::bytes(buf);
}

PYBIND11_MODULE(npu_writer, m) {
    m.doc() = "Node-segment packer/writer matching Python struct layout <h R s R s C s h c (little-endian)>";

    m.def("pack_and_write", &pack_and_write_py,
          py::arg("dat_path"),
          py::arg("write_pos"),
          py::arg("R"),
          py::arg("C"),
          py::arg("segments"),
          R"(Pack segments (offset, seq, bq, cigar, rq, strand) with fixed field sizes (R,C)
and write them at absolute file position write_pos in dat_path. Returns bytes written.)");

    m.def("pack_only", &pack_only_py,
          py::arg("R"),
          py::arg("C"),
          py::arg("segments"),
          R"(Pack segments and return the raw bytes (no I/O).)");
}
