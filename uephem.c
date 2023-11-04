/*
 * Copyright (c) 2022-2023 C. J. Howard
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#define DE_OFFSET_TIME 0xA5C
#define DE_OFFSET_TABLE1 0xA88
#define DE_OFFSET_DENUM 0xB18
#define DE_OFFSET_TABLE2 0xB1C
#define DE_ENDIAN_SWAP(denum) (denum&0xFFFF0000L)
#define DE_MIN_ITEM_ID 0
#define DE_MAX_ITEM_ID 14
#define DE_MAX_NCONST 400
#define DE_CNAME_LENGTH 6

/// Number of components for items 0-14.
const int DE_NCOMP[15] = {3,3,3,3,3,3,3,3,3,3,3,2,3,3,1};

const char* err_fopen = "file open failed";
const char* err_fread = "file read failed";
const char* err_bad_arg = "bad argument";
const char* err_no_item = "item not found";
const char* err_date = "date out of range";

/// Prints an error string and exits.
void error(const char* msg)
{	
	printf("error: %s\n", msg);
	exit(EXIT_FAILURE);
}

/// Swaps the byte order of a 32-bit number.
uint32_t swap32(uint32_t x)
{
	x = ((x << 8) & 0xFF00FF00) | ((x >> 8) & 0xFF00FF);
	return (x << 16) | (x >> 16);
}

/// Swaps the byte order of a 64-bit number.
uint64_t swap64(uint64_t x)
{
	x = ((x << 8) & 0xFF00FF00FF00FF00) | ((x >> 8) & 0x00FF00FF00FF00FF);
	x = ((x << 16) & 0xFFFF0000FFFF0000) | ((x >> 16) & 0x0000FFFF0000FFFF);
	return (x << 32) | (x >> 32);
}

/// Calls fread then swaps the byte order of 32 and 64-bit data.
size_t fread_swap(void* ptr, size_t size, size_t count, FILE* stream)
{
	size_t status = fread(ptr, size, count, stream);
	if (size == sizeof(uint64_t))
		for (uint64_t* ptr64 = (uint64_t*)ptr; count; --count, ++ptr64)
			*ptr64 = swap64(*ptr64);
	else if (size == sizeof(uint32_t))
		for (uint32_t* ptr32 = (uint32_t*)ptr; count; --count, ++ptr32)
			*ptr32 = swap32(*ptr32);
	return status;
}

/// Evaluates a Chebyshev polynomial.
double chebyshev(const double* a, int n, double x)
{
	double y = *(a++);
	y += *(a++) * x;
	
	double n2 = 1.0, n1 = x, n0;
	x *= 2.0;
	
	for (n -= 2; n; --n, n2 = n1, n1 = n0)
	{
		n0 = x * n1 - n2;
		y += *(a++) * n0;
	}
	
	return y;
}

/// Evaluates the derivative of a Chebyshev polynomial.
double chebyshev_derivative(const double* a, int n, double x)
{
	double y = 0.0;
	double n2 = 0.0, n1 = 1.0, n0;
	
	for (int i = 1; i < n; ++i)
	{
		y += *(a + i) * i * n1;
		n0 = 2.0 * x * n1 - n2;
		n2 = n1;
		n1 = n0;
	}
	
	return y;
}

int main(int argc, char* argv[])
{
	if (argc != 4 && argc != 6)
	{
		printf("usage: uephem <file> <item ID> <t0> [<t1> <resolution>]\n");
		return EXIT_FAILURE;
	}
	
	// Open ephemeris file
	FILE* file;
	file = fopen(argv[1], "rb");
	if (!file)
		error(err_fopen);
	
	// Parse item ID
	char* endptr;
	int32_t item_id = strtol(argv[2], &endptr, 0);
	if (*endptr != '\0')
		error(err_bad_arg);
	if (item_id < DE_MIN_ITEM_ID || item_id > DE_MAX_ITEM_ID)
		error(err_no_item);
	
	// Parse JD start time
	double jd_start = strtod(argv[3], &endptr);
	if (*endptr != '\0')
		error(err_bad_arg);
	
	// Init JD end time and time resolution to single time output
	double jd_end = jd_start;
	int resolution = 1;
	
	if (argc == 6)
	{
		// Parse JD end time
		jd_end = strtod(argv[4], &endptr);
		if (*endptr != '\0')
			error(err_bad_arg);
		
		// Parse time resolution
		resolution = strtol(argv[5], &endptr, 0);
		if (*endptr != '\0' || resolution <= 0)
			error(err_bad_arg);
	}
	
	// Determine JD step size
	double jd_step = (jd_end - jd_start) / (resolution - 1);
	if (jd_start == jd_end)
	{
		jd_step = 0.0;
		resolution = 1;
	}
	else if (resolution == 1)
	{
		jd_start = (jd_start + jd_end) * 0.5;
		jd_step = 0.0;
	}
	
	// Read DE version number
	int32_t denum;
	fseek(file, DE_OFFSET_DENUM, SEEK_SET);
	fread(&denum, sizeof(int32_t), 1, file);
	if (ferror(file))
		error(err_fread);
	
	// If the DE number has any data in its most significant word, swap the byte order of following fread calls
	size_t (*fread_ntoh)(void*, size_t, size_t, FILE*) = &fread;
	if (DE_ENDIAN_SWAP(denum))
		fread_ntoh = &fread_swap;
	
	// Read file time
	double file_time[3];
	fseek(file, DE_OFFSET_TIME, SEEK_SET);
	fread_ntoh(file_time, sizeof(double), 3, file);
	
	// Check if time parameters are within the file's timespan
	if (jd_start < file_time[0] || jd_start > file_time[1] ||
		jd_end < file_time[0] || jd_end > file_time[1])
		error(err_date);
	
	// Read and combine coefficient tables
	int32_t nconst;
	int32_t table[15][3];
	fread_ntoh(&nconst, sizeof(int32_t), 1, file);
	fseek(file, DE_OFFSET_TABLE1, SEEK_SET);
	fread_ntoh(table, sizeof(int32_t), 12 * 3, file);
	fseek(file, DE_OFFSET_TABLE2, SEEK_SET);
	fread_ntoh(&table[12][0], sizeof(int32_t), 3, file);
	if (nconst > DE_MAX_NCONST)
		fseek(file, (nconst - DE_MAX_NCONST) * DE_CNAME_LENGTH, SEEK_CUR);
	fread_ntoh(&table[13][0], sizeof(int32_t), 2 * 3, file);
	if (ferror(file))
		error(err_fread);
	
	// Check if the specified item has any coefficients
	if (!table[item_id][2])
		error(err_no_item);
	
	// Determine number of coefficients per record
	long rec_ncoeff = 0;
	for (int i = 0; i < 15; ++i)
	{
		long ncoeff = table[i][0] + table[i][1] * table[i][2] * DE_NCOMP[i] - 1;
		if (ncoeff > rec_ncoeff)
			rec_ncoeff = ncoeff;
	}
	
	// Allocate record buffer
	long rec_size = rec_ncoeff * sizeof(double);
	double* rec_buf = malloc(rec_size);
	
	// Get item's number of components and number of coefficients per component
	int32_t item_ncomp = DE_NCOMP[item_id];
	int32_t item_ncoeff = table[item_id][1];
	
	// Seek to first coefficient record
	fseek(file, rec_size * 2, SEEK_SET);
	long rec_index[2] = {-1, -1};
	
	// Init record-skipping constants
	const long max_rec_skip = LONG_MAX / rec_size;
	const long max_rec_skip_size = rec_size * max_rec_skip;
	
	// For each time point
	for (int i = 0; i < resolution; ++i)
	{
		// Calculate JD time
		double jd = jd_start + jd_step * i;
		
		// Determine coefficient record containing current time
		rec_index[1] = (long)((jd - file_time[0]) / file_time[2]);
		
		// If record index changed
		if (rec_index[0] != rec_index[1])
		{
			// Determine number of records to skip
			long rec_skip = rec_index[1] - rec_index[0] - 1;
			
			// Handle 32-bit limits when seeking through large files
			while (rec_skip >= max_rec_skip)
			{
				fseek(file, max_rec_skip_size, SEEK_CUR);
				rec_skip -= max_rec_skip;
			}
			
			// Seek to record
			if (rec_skip)
				fseek(file, rec_size * rec_skip, SEEK_CUR);
			
			// Read record into buffer
			fread_ntoh(rec_buf, sizeof(double), rec_ncoeff, file);
			if (ferror(file))
				error(err_fread);
			
			rec_index[0] = rec_index[1];
		}
		
		// Get index of the subinterval for the item at the given JD
		double subinterval_duration = file_time[2] / (double)table[item_id][2];
		int subinterval_index = (int)((jd - rec_buf[0]) / subinterval_duration);
		
		// Remap JD to the Chebyshev domain [-1, 1]
		double subinterval_start = rec_buf[0] + subinterval_index * subinterval_duration;
		double t = (jd - subinterval_start) / subinterval_duration * 2.0 - 1.0;
		
		// Pointer to the first coefficient of the first property
		const double* coeffs_start = rec_buf + (table[item_id][0] - 1) + subinterval_index * item_ncoeff * item_ncomp;
		
		// Print Julian date
		printf("%.*f", DBL_DECIMAL_DIG, jd);
		
		// For each component in the item
		const double* coeffs = coeffs_start;
		for (int j = 0; j < item_ncomp; ++j, coeffs += item_ncoeff)
		{
			// Evaluate the Chebyshev polynomial and output the result
			double x = chebyshev(coeffs, item_ncoeff, t);
			printf(",%.*e", DBL_DECIMAL_DIG, x);
		}
		
		if (item_id < 13)
		{
			// For each component in the item
			coeffs = coeffs_start;
			for (int j = 0; j < item_ncomp; ++j, coeffs += item_ncoeff)
			{
				// Evaluate the derivative of the Chebyshev polynomial w.r.t. time and output the result
				double dx = chebyshev_derivative(coeffs, item_ncoeff, t) / subinterval_duration * 2.0;
				printf(",%.*e", DBL_DECIMAL_DIG, dx);
			}
		}
		
		printf("\n");
	}
	
	// Clean up
	free(rec_buf);
	fclose(file);
	
	return EXIT_SUCCESS;
}
