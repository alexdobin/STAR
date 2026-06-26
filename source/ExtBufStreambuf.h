#ifndef CODE_ExtBufStreambuf
#define CODE_ExtBufStreambuf

#include <streambuf>
#include <ios>

// A std::streambuf backed by a caller-owned, fixed-size char buffer (zero-copy).
//
// This replaces the non-portable idiom
//     someStringStream.rdbuf()->pubsetbuf(externalBuffer, size);
// which STAR used to read/write reads and SAM/BAM records directly from/to a
// pre-allocated char array without copying.
//
// std::basic_stringbuf::setbuf() (what pubsetbuf() calls on a string stream) has
// implementation-defined behavior: libstdc++ (GCC/Linux) adopts the external
// buffer, but libc++ (clang/macOS) ignores it entirely. With libc++ the streams
// therefore stayed empty, which made STAR read 0 input reads and write empty
// SAM/BAM output. This class implements the buffer adoption explicitly so the
// behavior is identical on every standard library.
class ExtBufStreambuf : public std::streambuf {
public:
    // Use the external buffer as the get (read) area: [base, base+n).
    void setReadBuffer(char* base, std::streamsize n) {
        setg(base, base, base + n);
    }

    // Use the external buffer as the put (write) area: [base, base+n).
    void setWriteBuffer(char* base, std::streamsize n) {
        setp(base, base + n);
    }

protected:
    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                     std::ios_base::openmode which) override {
        if (which & std::ios_base::in) {
            char* base = eback();
            char* end  = egptr();
            char* newp = (dir == std::ios_base::beg) ? base + off
                       : (dir == std::ios_base::cur) ? gptr() + off
                       :                               end  + off;
            if (newp < base) newp = base;
            if (newp > end)  newp = end;
            setg(base, newp, end);
            return pos_type(newp - base);
        }
        if (which & std::ios_base::out) {
            char* base = pbase();
            char* end  = epptr();
            char* newp = (dir == std::ios_base::beg) ? base + off
                       : (dir == std::ios_base::cur) ? pptr() + off
                       :                               end  + off;
            if (newp < base) newp = base;
            if (newp > end)  newp = end;
            setp(base, end);
            pbump(static_cast<int>(newp - base));
            return pos_type(newp - base);
        }
        return pos_type(off_type(-1));
    }

    pos_type seekpos(pos_type sp, std::ios_base::openmode which) override {
        return seekoff(off_type(sp), std::ios_base::beg, which);
    }
};

#endif
