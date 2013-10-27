#ifndef SAX_H
#define SAX_H

#define SAX_ALPHABET "ACGTVLIMFWPSYNQDEKRHX"

char *sax_from_array_force_window(double *array, size_t length, size_t alphabet_size, size_t window_size, double mean, double std);
/* given an array, the length of the array, and the desired SAX discretization size, */
/* make a new version of the array that is the SAX symbolic version. this one also allows the */
/* window size to be fixed */

char *sax_from_array(double *array, size_t length, size_t alphabet_size);
/* given an array, the length of the array, and the desired SAX discretization size, */
/* make a new version of the array that is the SAX symbolic version. */

#endif
