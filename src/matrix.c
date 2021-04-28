#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high)
{
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, unsigned int seed, double low, double high)
{
    srand(seed);
    for (int i = 0; i < result->rows; i++)
    {
        for (int j = 0; j < result->cols; j++)
        {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols)
{
    if (rows <= 0 || cols <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "error");
        return -1;
    }

    (*mat) = (matrix *)malloc(sizeof(matrix));

    if (!(*mat))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to allocate space for new matrix structure");
        return -1;
    }

    (*mat)->data = (double *)calloc(rows * cols, sizeof(double));

    if (!(*mat)->data)
    {
        PyErr_SetString(PyExc_ValueError, "Failed to allocate space for data within matrix struct");
        return -1;
    }

    (*mat)->parent = NULL;
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->ref_cnt = 1;

    return 0;
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`.
 * You should return -1 if either `rows` or `cols` or both are non-positive or if any
 * call to allocate memory in this function fails.
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols)
{
    /* TODO: YOUR CODE HERE */
    if (rows <= 0 || cols <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "error");
        return -1;
    }

    (*mat) = (matrix *)malloc(sizeof(matrix));

    if (!(*mat))
    {
        PyErr_SetString(PyExc_ValueError, "Failed to allocate space for new matrix structure");
        return -1;
    }

    (*mat)->data = from->data + offset;
    (*mat)->parent = from;
    from->ref_cnt = from->ref_cnt + 1;

    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->ref_cnt = 1;
    return 0;
}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or if `mat` is the last existing slice of its parent matrix and its parent matrix has no other references
 * (including itself). You cannot assume that mat is not NULL.
 */
void deallocate_matrix(matrix *mat)
{
    /* TODO: YOUR CODE HERE */
    if (mat == NULL)
    {
        PyErr_SetString(PyExc_ValueError, "error");
        return -1;
    }
    if (!(mat->parent) && (mat->ref_cnt == 1))
    {
        free(mat->data);
    }

    if (mat->parent && mat->parent->ref_cnt > 1)
    {
        mat->parent->ref_cnt -= 1;
    }

    if (mat->parent && mat->parent->ref_cnt == 1)
    {
        free(mat->parent->data);
    }

    free(mat);
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col)
{
    /* TODO: YOUR CODE HERE */
    if (!mat)
    {
        PyErr_SetString(PyExc_ValueError, "error");
        return -1;
    }
    double *mat_data = mat->data;
    int total_offset = ((row) * (mat->cols)) + (col);

    return *(mat_data + total_offset);
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val)
{
    /* TODO: YOUR CODE HERE */
    if (!mat)
    {
        PyErr_SetString(PyExc_ValueError, "error");
        return -1;
    }
    double *mat_data = mat->data;
    int total_offset = ((row) * (mat->cols)) + (col);

    *(mat_data + total_offset) = val;
}

/*
 * Sets all entries in mat to val
 */
void fill_matrix(matrix *mat, double val)
{
    double *data_ptr = mat->data;
    /* TODO: YOUR CODE HERE */
    for (int i = 0; i < (mat->rows * mat->cols) / 4 * 4; i += 4)
    {
        *(data_ptr + i) = val;
        *(data_ptr + i + 1) = val;
        *(data_ptr + i + 2) = val;
        *(data_ptr + i + 3) = val;
    }
    for (int i = (mat->rows * mat->cols) / 4 * 4;
         i < (mat->rows * mat->cols); i++)
    {
        *(data_ptr + i) = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2)
{

    /* TODO: YOUR CODE HERE */
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    double *result_ptr = result->data;
    double *mat1_ptr = mat1->data;
    double *mat2_ptr = mat2->data;

    __m256d mat1_simd = _mm256_setzero_pd();
    __m256d mat2_simd = _mm256_setzero_pd();
    __m256d add_result = _mm256_setzero_pd();

#pragma omp parallel for
    for (int i = thread_id; i < (mat1->rows * mat1->cols) / 16 * 16; i += (16 * num_threads))
    {

        __m256d *start_0_mat1 = (__m256d *)(mat1_ptr + i);
        __m256d *start_0_mat2 = (__m256d *)(mat2_ptr + i);
        mat1_simd = _mm256_loadu_pd(start_0_mat1);
        mat2_simd = _mm256_loadu_pd(start_0_mat2);
        add_result = _mm256_add_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i, add_result);

        __m256d *start_1_mat1 = (__m256d *)(mat1_ptr + i + 4);
        __m256d *start_1_mat2 = (__m256d *)(mat2_ptr + i + 4);
        mat1_simd = _mm256_loadu_pd(start_1_mat1);
        mat2_simd = _mm256_loadu_pd(start_1_mat2);
        add_result = _mm256_add_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 4, add_result);

        __m256d *start_2_mat1 = (__m256d *)(mat1_ptr + i + 8);
        __m256d *start_2_mat2 = (__m256d *)(mat2_ptr + i + 8);
        mat1_simd = _mm256_loadu_pd(start_2_mat1);
        mat2_simd = _mm256_loadu_pd(start_2_mat2);
        add_result = _mm256_add_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 8, add_result);

        __m256d *start_3_mat1 = (__m256d *)(mat1_ptr + i + 12);
        __m256d *start_3_mat2 = (__m256d *)(mat2_ptr + i + 12);
        mat1_simd = _mm256_loadu_pd(start_3_mat1);
        mat2_simd = _mm256_loadu_pd(start_3_mat2);
        add_result = _mm256_add_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 12, add_result);
    }

#pragma omp parallel for
    for (int i = (mat1->rows * mat1->cols) / 16 * 16 + thread_id;
         i < (mat1->rows * mat1->cols); i += num_threads)
    {
        *(result_ptr + i) = *(mat1_ptr + i) + *(mat2_ptr + i);
    }

    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2)
{

    /* TODO: YOUR CODE HERE */
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    double *result_ptr = result->data;
    double *mat1_ptr = mat1->data;
    double *mat2_ptr = mat2->data;

    __m256d mat1_simd = _mm256_setzero_pd();
    __m256d mat2_simd = _mm256_setzero_pd();
    __m256d add_result = _mm256_setzero_pd();

#pragma omp parallel for
    for (int i = thread_id; i < (mat1->rows * mat1->cols) / 16 * 16; i += (16 * num_threads))
    {

        __m256d *start_0_mat1 = (__m256d *)(mat1_ptr + i);
        __m256d *start_0_mat2 = (__m256d *)(mat2_ptr + i);
        mat1_simd = _mm256_loadu_pd(start_0_mat1);
        mat2_simd = _mm256_loadu_pd(start_0_mat2);
        add_result = _mm256_sub_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i, add_result);

        __m256d *start_1_mat1 = (__m256d *)(mat1_ptr + i + 4);
        __m256d *start_1_mat2 = (__m256d *)(mat2_ptr + i + 4);
        mat1_simd = _mm256_loadu_pd(start_1_mat1);
        mat2_simd = _mm256_loadu_pd(start_1_mat2);
        add_result = _mm256_sub_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 4, add_result);

        __m256d *start_2_mat1 = (__m256d *)(mat1_ptr + i + 8);
        __m256d *start_2_mat2 = (__m256d *)(mat2_ptr + i + 8);
        mat1_simd = _mm256_loadu_pd(start_2_mat1);
        mat2_simd = _mm256_loadu_pd(start_2_mat2);
        add_result = _mm256_sub_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 8, add_result);

        __m256d *start_3_mat1 = (__m256d *)(mat1_ptr + i + 12);
        __m256d *start_3_mat2 = (__m256d *)(mat2_ptr + i + 12);
        mat1_simd = _mm256_loadu_pd(start_3_mat1);
        mat2_simd = _mm256_loadu_pd(start_3_mat2);
        add_result = _mm256_sub_pd(mat1_simd, mat2_simd);
        _mm256_storeu_pd(result_ptr + i + 12, add_result);
    }

#pragma omp parallel for
    for (int i = (mat1->rows * mat1->cols) / 16 * 16 + thread_id;
         i < (mat1->rows * mat1->cols); i += num_threads)
    {
        *(result_ptr + i) = *(mat1_ptr + i) - *(mat2_ptr + i);
    }

    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */

int mul_matrix(matrix *result, matrix *mat1, matrix *mat2)
{
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    double *result_ptr = result->data;
    int M = mat1->rows;
    int N = mat2->cols;
    int L = mat1->cols;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            *(result->data + (i * N) + j) = 0;
            for (int k = 0; k < L; k++)
            {
                *(result->data + (i * N) + j) += *(mat1->data + (i * L) + k) * *(mat2->data + (k * N) + j);
            }
        }
    }
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow)
{
    matrix *result_matrix = NULL;
    allocate_matrix(&result_matrix, mat->rows, mat->cols);

    int count = 0;
    if (pow == 0)
    {
        for (int i = 0; i < mat->rows; i++)
        {
            for (int j = 0; j < mat->rows; j++)
            {
                if (count == j)
                {
                    *(result_matrix->data + (i * mat->rows) + j) = 1;
                }
                else
                {
                    *(result_matrix->data + (i * mat->rows) + j) = 0;
                }
            }
            count++;
        }
        memcpy(result->data, result_matrix->data, sizeof(double) * mat->rows * mat->cols);
        return result;
    }

    memcpy(result_matrix->data, mat->data, sizeof(double) * mat->rows * mat->cols);

    for (int i = 1; i < pow; i++)
    {
        matrix *old_result = NULL;
        allocate_matrix(&old_result, mat->rows, mat->cols);
        memcpy(old_result->data, result_matrix->data, sizeof(double) * mat->rows * mat->cols);
        mul_matrix(result_matrix, mat, old_result);
    }
    memcpy(result->data, result_matrix->data, sizeof(double) * mat->rows * mat->cols);
    return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat)
{
    /* TODO: YOUR CODE HERE */
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    double *result_ptr = result->data;
    double *mat_ptr = mat->data;

    __m256d mat_simd = _mm256_setzero_pd();
    __m256d neg_result = _mm256_setzero_pd();
    __m256d all_zero = _mm256_setzero_pd();

    int total_size = (mat->rows * mat->cols);

#pragma omp parallel for
    for (int i = thread_id; i < (total_size / 16) * 16; i += (16 * num_threads))
    {

        __m256d *start_0_mat = (__m256d *)(mat_ptr + i);
        mat_simd = _mm256_loadu_pd(start_0_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0.0), mat_simd);
        _mm256_storeu_pd(result_ptr + i, neg_result);

        __m256d *start_1_mat = (__m256d *)(mat_ptr + i + 4);
        mat_simd = _mm256_loadu_pd(start_1_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0.0), mat_simd);
        _mm256_storeu_pd(result_ptr + i + 4, neg_result);

        __m256d *start_2_mat = (__m256d *)(mat_ptr + i + 8);
        mat_simd = _mm256_loadu_pd(start_2_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0.0), mat_simd);
        _mm256_storeu_pd(result_ptr + i + 8, neg_result);

        __m256d *start_3_mat = (__m256d *)(mat_ptr + i + 12);
        mat_simd = _mm256_loadu_pd(start_3_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0.0), mat_simd);
        _mm256_storeu_pd(result_ptr + i + 12, neg_result);
    }

#pragma omp parallel for
    for (int i = (total_size / 16) * 16 + thread_id;
         i < total_size; i += num_threads)
    {
        *(result_ptr + i) = -(*(mat_ptr + i));
    }

    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */

int abs_matrix(matrix *result, matrix *mat)
{
    /* TODO: YOUR CODE HERE */

    double *result_ptr = result->data;
    double *mat_ptr = mat->data;

    __m256d mat_simd = _mm256_setzero_pd();
    __m256d neg_result = _mm256_setzero_pd();
    __m256d abs_result = _mm256_setzero_pd();
    __m256d all_zero = _mm256_setzero_pd();

    int total_size = (mat->rows * mat->cols);
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

#pragma omp parallel for

    for (int i = thread_id; i < (total_size / 16) * 16; i += (16 * num_threads))
    {

        __m256d *start_0_mat = (__m256d *)(mat_ptr + i);
        mat_simd = _mm256_loadu_pd(start_0_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0), mat_simd);
        abs_result = _mm256_max_pd(neg_result, mat_simd);
        _mm256_storeu_pd(result_ptr + i, abs_result);

        __m256d *start_1_mat = (__m256d *)(mat_ptr + i + 4);
        mat_simd = _mm256_loadu_pd(start_1_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0), mat_simd);
        abs_result = _mm256_max_pd(neg_result, mat_simd);
        _mm256_storeu_pd(result_ptr + i + 4, abs_result);

        __m256d *start_2_mat = (__m256d *)(mat_ptr + i + 8);
        mat_simd = _mm256_loadu_pd(start_2_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0), mat_simd);
        abs_result = _mm256_max_pd(neg_result, mat_simd);
        _mm256_storeu_pd(result_ptr + i + 8, abs_result);

        __m256d *start_3_mat = (__m256d *)(mat_ptr + i + 12);
        mat_simd = _mm256_loadu_pd(start_3_mat);
        neg_result = _mm256_sub_pd(_mm256_set1_pd(0), mat_simd);
        abs_result = _mm256_max_pd(neg_result, mat_simd);
        _mm256_storeu_pd(result_ptr + i + 12, abs_result);
    }

    for (int i = (total_size / 16) * 16;
         i < total_size; i++)
    {
        *(result_ptr + i) = abs(*(mat_ptr + i));
    }

    return 0;
}

