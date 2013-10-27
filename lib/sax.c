#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <float.h>
#include <math.h>
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "sax.h"

static double *sax_make_cut_points_array(int alphabet_size)
/* Incredibly pathetic function to make the all-important array. */
{
    double *cut_points;
    AllocArray(cut_points, 20);
    cut_points[0] = -DBL_MAX;
    switch (alphabet_size)
    {
    case 3:
	cut_points[1] = -0.43;
	cut_points[2] = 0.43;
	break;
    case 4:
	cut_points[1] = -0.67;
	cut_points[2] = 0;
	cut_points[3] = 0.67;
	break;
    case 5:
	cut_points[1] = -0.84;
	cut_points[2] = -0.25;
	cut_points[3] = 0.25;
	cut_points[4] = 0.84;
	break;
    case 6:
	cut_points[1] = -0.97;
	cut_points[2] = -0.43;
	cut_points[3] = 0;
	cut_points[4] = 0.43;
	cut_points[5] = 0.97;
	break;
    case 7:
	cut_points[1] = -1.07;
	cut_points[2] = -0.57;
	cut_points[3] = -0.18;
	cut_points[4] = 0.18;
	cut_points[5] = 0.57;
	cut_points[6] = 1.07;
	break;
    case 8:
	cut_points[1] = -1.15;
	cut_points[2] = -0.67;
	cut_points[3] = -0.32;
	cut_points[4] = 0;
	cut_points[5] = 0.32;
	cut_points[6] = 0.67;
	cut_points[7] = 1.15;
	break;
    case 9:
	cut_points[1] = -1.22;
	cut_points[2] = -0.76;
	cut_points[3] = -0.43;
	cut_points[4] = -0.14;
	cut_points[5] = 0.14;
	cut_points[6] = 0.43;
	cut_points[7] = 0.76;
	cut_points[8] = 1.22;
	break;
    case 10:
	cut_points[1] = -1.28;
	cut_points[2] = -0.84;
	cut_points[3] = -0.52;
	cut_points[4] = -0.25;
	cut_points[5] = 0.;
	cut_points[6] = 0.25;
	cut_points[7] = 0.52;
	cut_points[8] = 0.84;
	cut_points[9] = 1.28;
	break;
    case 11:
	cut_points[1] = -1.34;
	cut_points[2] = -0.91;
	cut_points[3] = -0.6;
	cut_points[4] = -0.35;
	cut_points[5] = -0.11;
	cut_points[6] = 0.11;
	cut_points[7] = 0.35;
	cut_points[8] = 0.6;
	cut_points[9] = 0.91;
	cut_points[10] = 1.34;
	break;
    case 12:
	cut_points[1] = -1.38;
	cut_points[2] = -0.97;
	cut_points[3] = -0.67;
	cut_points[4] = -0.43;
	cut_points[5] = -0.21;
	cut_points[6] = 0;
	cut_points[7] = 0.21;
	cut_points[8] = 0.43;
	cut_points[9] = 0.67;
	cut_points[10] = 0.97;
	cut_points[11] = 1.38;
	break;
    case 13:
	cut_points[1] = -1.43;
	cut_points[2] = -1.02;
	cut_points[3] = -0.74;
	cut_points[4] = -0.5;
	cut_points[5] = -0.29;
	cut_points[6] = -0.1;
	cut_points[7] = 0.1;
	cut_points[8] = 0.29;
	cut_points[9] = 0.5;
	cut_points[10] = 0.74;
	cut_points[11] = 1.02;
	cut_points[12] = 1.43;
	break;
    case 14:
	cut_points[1] = -1.47;
	cut_points[2] = -1.07;
	cut_points[3] = -0.79;
	cut_points[4] = -0.57;
	cut_points[5] = -0.37;
	cut_points[6] = -0.18;
	cut_points[7] = 0;
	cut_points[8] = 0.18;
	cut_points[9] = 0.37;
	cut_points[10] = 0.57;
	cut_points[11] = 0.79;
	cut_points[12] = 1.07;
	cut_points[13] = 1.47;
	break;
    case 15:
	cut_points[1] = -1.5;
	cut_points[2] = -1.11;
	cut_points[3] = -0.84;
	cut_points[4] = -0.62;
	cut_points[5] = -0.43;
	cut_points[6] = -0.25;
	cut_points[7] = -0.08;
	cut_points[8] = 0.08;
	cut_points[9] = 0.25;
	cut_points[10] = 0.43;
	cut_points[11] = 0.62;
	cut_points[12] = 0.84;
	cut_points[13] = 1.11;
	cut_points[14] = 1.5;
	break;
    case 16:
	cut_points[1] = -1.53;
	cut_points[2] = -1.15;
	cut_points[3] = -0.89;
	cut_points[4] = -0.67;
	cut_points[5] = -0.49;
	cut_points[6] = -0.32;
	cut_points[7] = -0.16;
	cut_points[8] = 0;
	cut_points[9] = 0.16;
	cut_points[10] = 0.32;
	cut_points[11] = 0.49;
	cut_points[12] = 0.67;
	cut_points[13] = 0.89;
	cut_points[14] = 1.15;
	cut_points[15] = 1.53;
	break;
    case 17:
	cut_points[1] = -1.56;
	cut_points[2] = -1.19;
	cut_points[3] = -0.93;
	cut_points[4] = -0.72;
	cut_points[5] = -0.54;
	cut_points[6] = -0.38;
	cut_points[7] = -0.22;
	cut_points[8] = -0.07;
	cut_points[9] = 0.07;
	cut_points[10] = 0.22;
	cut_points[11] = 0.38;
	cut_points[12] = 0.54;
	cut_points[13] = 0.72;
	cut_points[14] = 0.93;
	cut_points[15] = 1.19;
	cut_points[16] = 1.56;
	break;
    case 18:
	cut_points[1] = -1.59;
	cut_points[2] = -1.22;
	cut_points[3] = -0.97;
	cut_points[4] = -0.76;
	cut_points[5] = -0.59;
	cut_points[6] = -0.43;
	cut_points[7] = -0.28;
	cut_points[8] = -0.14;
	cut_points[9] = 0;
	cut_points[10] = 0.14;
	cut_points[11] = 0.28;
	cut_points[12] = 0.43;
	cut_points[13] = 0.59;
	cut_points[14] = 0.76;
	cut_points[15] = 0.97;
	cut_points[16] = 1.22;
	cut_points[17] = 1.59;
	break;
    case 19:
	cut_points[1] = -1.62;
	cut_points[2] = -1.25;
	cut_points[3] = -1;
	cut_points[4] = -0.8;
	cut_points[5] = -0.63;
	cut_points[6] = -0.48;
	cut_points[7] = -0.34;
	cut_points[8] = -0.2;
	cut_points[9] = -0.07;
	cut_points[10] = 0.07;
	cut_points[11] = 0.2;
	cut_points[12] = 0.34;
	cut_points[13] = 0.48;
	cut_points[14] = 0.63;
	cut_points[15] = 0.8;
	cut_points[16] = 1;
	cut_points[17] = 1.25;
	cut_points[18] = 1.62;
	break;
    case 20:
	cut_points[1] = -1.64;
	cut_points[2] = -1.28;
	cut_points[3] = -1.04;
	cut_points[4] = -0.84;
	cut_points[5] = -0.67;
	cut_points[6] = -0.52;
	cut_points[7] = -0.39;
	cut_points[8] = -0.25;
	cut_points[9] = -0.13;
	cut_points[10] = 0;
	cut_points[11] = 0.13;
	cut_points[12] = 0.25;
	cut_points[13] = 0.39;
	cut_points[14] = 0.52;
	cut_points[15] = 0.67;
	cut_points[16] = 0.84;
	cut_points[17] = 1.04;
	cut_points[18] = 1.28;
	cut_points[19] = 1.64;
	break;
    default:
	errAbort("bad alphabet size.  need 2-20");
    }
    return cut_points;
}

static void sax_z_norm(double *array, size_t len, double mean, double std)
/* z-normalize an array in place. (needs GSL) */
{
    double stdd = (std > 0) ? std : 1;
    double meen = (std > 0) ? mean : 0;
    int i;
    for (i = 0; i < len; i++)
	array[i] = (array[i] - meen)/stdd;
}

static int sax_sum_ints(int * const array, size_t len)
/*  */
{
    int sum = 0;
    int i;
    for (i = 0; i < len; i++)
	sum += array[i];
    return sum;
}

static int *sax_thresh_cuts(double * const cut_points, size_t len, double val)
/* make an array of ones and zeroes the same size as the cut points array */
/* depending on the value given */
/* this one could be improved because it allocates a new array at each */
/* datapoint in the loop of the calling function. */
{
    int *array;
    int i;
    AllocArray(array, len);
    for (i = 0; i < len; i++)
	array[i] = (cut_points[i] <= val) ? 1 : 0;
    return array;
}

static void sax_make_ints(int *array, double * const cut_points, size_t len, double * const z_normed_data, size_t data_len)
/* pretty simple it just sums threshold array */
{
    int i;
    for (i = 0; i < data_len; i++)
    {
	int *tc = sax_thresh_cuts(cut_points, len, z_normed_data[i]);
	array[i] = sax_sum_ints(tc, len);
	freeMem(tc);
    }
}

static int *sax_overlapping(double *array, size_t length, double * const cut_points, size_t alphabet_size, 
			    size_t window_size, double mean, double std)
/* return an array of ints to be converted into the sax chars.  this one does much of the work because */
/* it runs SAX more than once and combines the overlapping arrays.  */
{
    int max_sax_len = pow(2, (int)log2((int)length));
    int i = 0;
    double *z_norm;
    int *sax_ints;
    int *ret;
    if (window_size > 0)
	max_sax_len = window_size;
    AllocArray(z_norm, max_sax_len);
    AllocArray(sax_ints, max_sax_len);
    AllocArray(ret, length);
    while (i + max_sax_len <= length)
    {
	int j, k;
	for (j = 0; j < max_sax_len; j++)
	    z_norm[j] = array[i + j];
	sax_z_norm(z_norm, max_sax_len, mean, std);
	sax_make_ints(sax_ints, cut_points, alphabet_size, z_norm, max_sax_len);
	/* copy into the new array */
	k = 0;
	if (i > 0)
	    k = max_sax_len/4;
	for (j = i + k; j < i + max_sax_len; j++)
	    ret[j] = sax_ints[k++];
	/* set i for the next iteration... if i+max is past the end of length, */
	/* make it exactly the length */
	if (i + max_sax_len == length)
	    break;
	i += max_sax_len/2;
	if (i + max_sax_len > length)
	    i = length - max_sax_len;
    }
    freez(&z_norm);
    freez(&sax_ints);
    return ret;
}

boolean is_pow2(int number)
/* determine whether number is a valid power of 2. */
{
    int x = number;
    while (x > 1)
    {
	if ((x & 1) == 1)
	    return FALSE;
	x = x >> 1;
    }
    return TRUE;
}

/** "public" functions **/

char *sax_from_array_force_window(double *array, size_t length, size_t alphabet_size, size_t window_size, double mean, double std)
/* given an array, the length of the array, and the desired SAX discretization size, */
/* make a new version of the array that is the SAX symbolic version. this one also allows the */
/* window size to be fixed */
{
    char *alphabet = SAX_ALPHABET;
    char *sax;
    double *cut_points = sax_make_cut_points_array(alphabet_size);
    int *sax_array;
    int i;
    /* check window_size (if given) is a power of 2 */
    if ((window_size > 0) && (!is_pow2((int)window_size)))
	errAbort("sax error: %d not a valid window size, must be a power of 2.", (int)window_size);
    sax_array = sax_overlapping(array, length, cut_points, alphabet_size, window_size, mean, std);
    AllocArray(sax, length + 1);
    for (i = 0; i < length; i++)
	sax[i] = alphabet[sax_array[i]-1];
    sax[length] = '\0';
    freez(&cut_points);
    freez(&sax_array);
    return sax;
}

char *sax_from_array(double *array, size_t length, size_t alphabet_size)
/* given an array, the length of the array, and the desired SAX discretization size, */
/* make a new version of the array that is the SAX symbolic version. */
{
    return sax_from_array_force_window(array, length, alphabet_size, 0, 0, -1);
}
