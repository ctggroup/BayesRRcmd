#include "marker.h"

#include "compression.h"

Marker::~Marker()
{

}

CompressedMarker Marker::compress() const
{
    // Prepare a stream to write into
    const auto bufferSize = static_cast<unsigned int>(size());
    std::unique_ptr<char[]> buffer;
    buffer.reset(new char[bufferSize]);

    std::ostringstream stream;
    stream.rdbuf()->pubsetbuf(buffer.get(), size());

    // Write the marker to the stream
    write(&stream);

    // Compress the stream
    const auto maxCompressedOutputSize = maxCompressedDataSize<char>(bufferSize);

    CompressedMarker compressed;
    compressed.buffer.reset(new unsigned char[maxCompressedOutputSize]);
    compressed.index.originalSize = static_cast<unsigned long>(size());
    compressed.index.compressedSize = compressData(buffer.get(),
                                                   bufferSize,
                                                   compressed.buffer.get(),
                                                   maxCompressedOutputSize);
    return compressed;
}

void Marker::decompress(unsigned char *data, const IndexEntry &index)
{
    // Prepare a buffer to decompress into
    const auto bufferSize = static_cast<unsigned int>(index.originalSize);
    std::unique_ptr<char[]> buffer;
    buffer.reset(new char[bufferSize]);

    // Decompress into the buffer
    extractData(data + index.pos,
                static_cast<unsigned int>(index.compressedSize),
                reinterpret_cast<unsigned char*>(buffer.get()),
                bufferSize);

    // Prepare a stream to read from
    std::istringstream stream;
    stream.rdbuf()->pubsetbuf(buffer.get(),
                              static_cast<std::streamsize>(index.originalSize));

    // Read the marker from the stream
    read(&stream);
}
