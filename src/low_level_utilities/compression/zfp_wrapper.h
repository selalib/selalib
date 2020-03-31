/**
 * Wrapper functions providing access to ZFP via an LZ4-like interface.
 * 2016, Klaus Reuter, khr@rzg.mpg.de
 */

unsigned int default_precision = 32;

int zfp_maximum_size(double* source, int n_doubles_source, unsigned int precision);

int zfp_compress_default(double* source, char* dest, int n_doubles_source, int max_n_bytes_dest);
int zfp_compress_safe(double* source, char* dest, int n_doubles_source, int max_n_bytes_dest, unsigned int precision);

int zfp_decompress_default(char* source, double* dest, int n_bytes_source, int max_n_doubles_dest);
int zfp_decompress_safe(char* source, double* dest, int n_bytes_source, int max_n_doubles_dest, unsigned int precision);
