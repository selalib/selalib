/**
 * Wrapper functions providing access to ZFP via an LZ4-like interface.
 * If the macro USE_ZFP is not defined (-> ZFPConfig.cmake),
 * the functions are stub functions, only.
 * 2016, Klaus Reuter, khr@rzg.mpg.de
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#ifdef USE_ZFP
#include <zfp.h>
#endif
#include "zfp_wrapper.h"

// return value: maximum length of the byte stream after compression
int zfp_maximum_size(double* source, int n_doubles_source, unsigned int precision) {
    int n_bytes = 0;
#ifdef USE_ZFP
    zfp_field* field;
    zfp_stream* zfp;
    const int nx = 4;
    const int ny = 4;
    int nz;
    if ((n_doubles_source % 64) == 0) {
        nz = n_doubles_source / (nx * ny);
        field = zfp_field_3d(source, zfp_type_double, nx, ny, nz);

        zfp = zfp_stream_open(NULL);
        zfp_stream_set_precision(zfp, precision);
        n_bytes = zfp_stream_maximum_size(zfp, field);

        zfp_field_free(field);
        zfp_stream_close(zfp);
    }
#endif
    return n_bytes;
}

// return value: length of the compressed byte stream
int zfp_compress_default(double* source, char* dest, int n_doubles_source, int max_n_bytes_dest) {
    unsigned int precision;
    precision = default_precision;
    return zfp_compress_safe(source, dest, n_doubles_source, max_n_bytes_dest, precision);
}

// return value: length of the compressed byte stream
int zfp_compress_safe(double* source, char* dest, int n_doubles_source, int max_n_bytes_dest, unsigned int precision) {
    int n_bytes = 0;
#ifdef USE_ZFP
    zfp_field* field;
    zfp_stream* zfp;
    bitstream* bit_stream;
    int idx;
    uint bits;
    const int nx = 4;
    const int ny = 4;
    int nz;
    if ((n_doubles_source % 64) == 0) {
        nz = n_doubles_source / (nx * ny);
        field = zfp_field_3d(source, zfp_type_double, nx, ny, nz);

        zfp = zfp_stream_open(NULL);
        bit_stream = stream_open(dest, max_n_bytes_dest);
        zfp_stream_set_bit_stream(zfp, bit_stream);
        zfp_stream_set_precision(zfp, precision);
        zfp_stream_rewind(zfp);

        zfp_compress(zfp, field);

        zfp_stream_flush(zfp);
        n_bytes = zfp_stream_compressed_size(zfp);

        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(bit_stream);
    }
#endif
    return n_bytes;
}

// return value: nonzero upon success
int zfp_decompress_default(char* source, double* dest, int n_bytes_source, int max_n_doubles_dest) {
    unsigned int precision;
    precision = default_precision;
    return zfp_decompress_safe(source, dest, n_bytes_source, max_n_doubles_dest, precision);
}

// return value: nonzero upon success
int zfp_decompress_safe(char* source, double* dest, int n_bytes_source, int max_n_doubles_dest, unsigned int precision) {
    int success = 0;
#ifdef USE_ZFP
    zfp_field* field;
    zfp_stream* zfp;
    bitstream* bit_stream;
    const int nx = 4;
    const int ny = 4;
    int nz;

    nz = max_n_doubles_dest / (nx * ny);
    field = zfp_field_3d(dest, zfp_type_double, nx, ny, nz);

    zfp = zfp_stream_open(NULL);
    bit_stream = stream_open(source, n_bytes_source);
    zfp_stream_set_bit_stream(zfp, bit_stream);
    zfp_stream_set_precision(zfp, precision);
    zfp_stream_rewind(zfp);

    success = zfp_decompress(zfp, field);

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(bit_stream);
#endif
    return success;
}
