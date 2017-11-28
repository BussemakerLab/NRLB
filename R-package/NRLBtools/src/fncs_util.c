/*<<<
 * ------------------------------------------------------------------------
 * Copyright (C) 2006-2017 The Trustees of Columbia University in the City
 * of New York. All rights reserved.
 *
 * The REDUCE Suite of software programs has been developed and actively
 * maintained by Dr. Xiang-Jun Lu at the Bussemaker Laboratory at the
 * Department of Biological Sciences, Columbia University. The suite has
 * its origin from the linear regression-based REDUCE algorithm of
 * Dr. Harmen Bussemaker to discover cis-regulatory elements. Dr. Barret
 * Foat extended REDUCE by adding an optimization procedure to get the
 * so-called Position Specific Affinity Matrix, and implemented his
 * algorithms in the MatrixREDUCE pipeline using Perl and GNU Scientific
 * Library. Dr. Lu has completely rewritten and significantly enhanced the
 * REDUCE/MatrixREDUCE codebase using ANSI C to make the programs efficient
 * and self-contained. The REDUCE Suite (now at v2.2) consists of a dozen
 * standalone, yet interconnected programs.
 *
 * This suite of software programs is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Click
 * Software License file (at the $REDUCE_SUITE/license directory, in either
 * Word, text, or PDF format) for more details.
 *
 * Please check http://reducesuite.bussemakerlab.org/ for more information.
 * ------------------------------------------------------------------------
>>>*/

#include "REDUCE_Suite.h"

void nrerror(char *error_text)
{
    fatal("%s\n", error_text);
}

void vector_boundary_check(long nl, long nh, char *fun_name)
{
    if (nl > nh)
        fatal("wrong vector boundary for %s: [%ld to %ld]\n", fun_name, nl, nh);
}

void matrix_boundary_check(long nrl, long nrh, long ncl, long nch, char *fun_name)
{
    if (nrl > nrh || ncl > nch)
        fatal("wrong matrix boundary for %s: [%ld to %ld; %ld to %ld]\n",
              fun_name, nrl, nrh, ncl, nch);
}

/* ------------------------------------------------------------------ */
void init_cvector(char *v, long nl, long nh, char c)
{
    long i;

    for (i = nl; i <= nh; i++)
        v[i] = c;
}

void init_cmatrix(char **m, long nrl, long nrh, long ncl, long nch, char c)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_cvector(m[i], ncl, nch, c);
}

static char *cvector_nr(long nl, long nh)
/* allocate a char vector with subscript range [nl..nh] */
{
    char *v;

    v = (char *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(char)));
    if (!v)
        nrerror("allocation failure in cvector()");

    return v - nl + NR_END;
}

char *cvector(long nl, long nh)
{
    char *v;

    vector_boundary_check(nl, nh, "cvector()");
    v = cvector_nr(nl, nh);
    init_cvector(v, nl, nh, '\0');

    return v;
}

void free_cvector(char *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);
}

static char **cmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    char **m;
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

    /* allocate pointers to rows */
    m = (char **) malloc((size_t) ((nrow + NR_END) * sizeof(char *)));
    if (!m)
        nrerror("allocation failure 1 in cmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(char)));
    if (!m[nrl])
        nrerror("allocation failure 2 in cmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
{
    char **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "cmatrix()");
    m = cmatrix_nr(nrl, nrh, ncl, nch);
    init_cmatrix(m, nrl, nrh, ncl, nch, '\0');

    return m;
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);
}

/* ------------------------------------------------------------------ */
static double *dvector_nr(long nl, long nh)
{
    double *v;

    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in dvector()");

    return v - nl + NR_END;
}

double *dvector(long nl, long nh)
{
    double *v;

    vector_boundary_check(nl, nh, "dvector()");
    v = dvector_nr(nl, nh);
    init_dvector(v, nl, nh, 0.0);

    return v;
}

void free_dvector(double *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);
}

static double **dmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    /* allocate pointers to rows */
    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure 1 in dmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in dmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "dmatrix()");
    m = dmatrix_nr(nrl, nrh, ncl, nch);
    init_dmatrix(m, nrl, nrh, ncl, nch, 0.0);

    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);
}

void diff_dvector(double *df, double *d1, double *d2, long nl, long nh)
{
    long i;

    for (i = nl; i <= nh; i++)
        df[i] = d1[i] - d2[i];
}

double ave_dvector(double *d, long nl, long nh)
{
    long i, n = nh - nl + 1;
    double dsum = 0.0;

    if (n <= 0)
        fatal("array indices [%ld --> %ld] for <ave_dvector()>\n", nl, nh);

    for (i = nl; i <= nh; i++)
        dsum += d[i];

    return dsum / n;
}

double var_dvector(double *d, long nl, long nh)
{
    long i, n = nh - nl + 1;
    double aved, temp, dsum = 0.0;

    if (n < 2)
        fatal("less than 2 items for <var_dvector()>\n");

    aved = ave_dvector(d, nl, nh);
    for (i = nl; i <= nh; i++) {
        temp = d[i] - aved;
        dsum += temp * temp;
    }

    return dsum / (n - 1);
}

double max_dvector(double *d, long nl, long nh)
{
    long i;
    double dmax = -XBIG;

    for (i = nl; i <= nh; i++)
        dmax = dval_max(dmax, d[i]);

    return dmax;
}

double min_dvector(double *d, long nl, long nh)
{
    long i;
    double dmin = XBIG;

    for (i = nl; i <= nh; i++)
        dmin = dval_min(dmin, d[i]);

    return dmin;
}

double norm_dvector(double *d, long nl, long nh)
{
    long i;
    double dsum = 0.0;

    for (i = nl; i <= nh; i++)
        dsum += d[i] * d[i];

    return sqrt(dsum);
}

double norm_inf_dvector(double *d, long nl, long nh)
{
    long i;
    double dmax = -XBIG;

    for (i = nl; i <= nh; i++)
        dmax = dval_max(dmax, fabs(d[i]));

    return dmax;
}

void copy_dvector(double *src, long nl, long nh, double *dst)
{
    long i;

    for (i = nl; i <= nh; i++)
        dst[i] = src[i];
}

void init_dvector(double *v, long nl, long nh, double d)
{
    long i;

    for (i = nl; i <= nh; i++)
        v[i] = d;
}

void init_dmatrix(double **m, long nrl, long nrh, long ncl, long nch, double d)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_dvector(m[i], ncl, nch, d);
}

void identity_dmatrix(double **m, long n)
{
    long i;

    for (i = 1; i <= n; i++) {
        init_dvector(m[i], 1, n, 0.0);
        m[i][i] = 1.0;
    }
}

void copy_dmatrix(double **a, long nrl, long nrh, long ncl, long nch, double **o)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        copy_dvector(a[i], ncl, nch, o[i]);
}

/* ------------------------------------------------------------------ */
void copy_lvector(long *src, long nl, long nh, long *dst)
{
    long i;

    for (i = nl; i <= nh; i++)
        dst[i] = src[i];
}

void init_lvector(long *v, long nrl, long nrh, long l)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        v[i] = l;
}

void init_lmatrix(long **m, long nrl, long nrh, long ncl, long nch, long l)
{
    long i;

    for (i = nrl; i <= nrh; i++)
        init_lvector(m[i], ncl, nch, l);
}

static long *lvector_nr(long nl, long nh)
{
    long *v;

    v = (long *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in lvector()");

    return v - nl + NR_END;
}

long *lvector(long nl, long nh)
{
    long *v;

    vector_boundary_check(nl, nh, "lvector()");
    v = lvector_nr(nl, nh);
    init_lvector(v, nl, nh, 0);

    return v;
}

void free_lvector(long *v, long nl, long nh)
{
    free((FREE_ARG) (v + nl - NR_END));
    UNUSED_PARAMETER(nh);
}

static long **lmatrix_nr(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    long **m;

    /* allocate pointers to rows */
    m = (long **) malloc((size_t) ((nrow + NR_END) * sizeof(long *)));
    if (!m)
        nrerror("allocation failure 1 in lmatrix()");
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (long *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(long)));
    if (!m[nrl])
        nrerror("allocation failure 2 in lmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return pointer to array of pointers to rows */
    return m;
}

long **lmatrix(long nrl, long nrh, long ncl, long nch)
{
    long **m;

    matrix_boundary_check(nrl, nrh, ncl, nch, "lmatrix()");
    m = lmatrix_nr(nrl, nrh, ncl, nch);
    init_lmatrix(m, nrl, nrh, ncl, nch, 0);

    return m;
}

void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
    UNUSED_PARAMETER(nrh);
    UNUSED_PARAMETER(nch);
}

/* ------------------------------------------------------------------ */
double dval_sqr(double dval)
{
    return (dval == 0.0) ? 0.0 : dval * dval;
}

double dval_max(double a, double b)
{
    return (a > b) ? a : b;
}

double dval_min(double a, double b)
{
    return (a < b) ? a : b;
}

long lval_max(long a, long b)
{
    return (a > b) ? a : b;
}

long lval_min(long a, long b)
{
    return (a < b) ? a : b;
}

long dval_in_range(double dval, double dlow, double dhigh)
{
    return dval >= dlow && dval <= dhigh;
}

long lval_in_range(long lval, long llow, long lhigh)
{
    return lval >= llow && lval <= lhigh;
}

/* only suitable from small array */
long lval_in_array(long lval, long ib, long ie, long *lvec)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (lval == lvec[i])
            return TRUE;

    return FALSE;
}

void dval_swap(double *pa, double *pb)
{
    double temp;

    temp = *pa;
    *pa = *pb;
    *pb = temp;
}

void lval_swap(long *pa, long *pb)
{
    long temp;

    temp = *pa;
    *pa = *pb;
    *pb = temp;
}

void cval_swap(char *pa, char *pb)
{
    int c;

    c = *pa;
    *pa = *pb;
    *pb = c;
}

void str_swap(char **src, char **dst)
{
    char *ptr;

    ptr = *dst;
    *dst = *src;
    *src = ptr;
}

void dvec_swap(double **src, double **dst)
{
    double *ptr;

    ptr = *dst;
    *dst = *src;
    *src = ptr;
}

int dval_compare(const void *v1, const void *v2)
{
    const double *p1, *p2;

    p1 = (const double *) v1;
    p2 = (const double *) v2;

    return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

int lval_compare(const void *v1, const void *v2)
{
    const long *p1, *p2;

    p1 = (const long *) v1;
    p2 = (const long *) v2;

    return (*p1 > *p2) ? 1 : (*p1 < *p2) ? -1 : 0;
}

int cstr_compare(const void *v1, const void *v2)
{
    const char **p1, **p2;

    p1 = (const char **) v1;
    p2 = (const char **) v2;

    return strcmp(*p1, *p2);
}

/* first by length, then content (case-insensitive) */
int strtags_compare(const void *v1, const void *v2)
{
    long n1, n2;
    const struct_tag *p1, *p2;

    p1 = (const struct_tag *) v1;
    p2 = (const struct_tag *) v2;

    n1 = strlen(p1->str);
    n2 = strlen(p2->str);

    if (n1 > n2)
        return 1;
    else if (n1 < n2)
        return -1;
    else
        return case_strcmp(p1->str, p2->str);
}

int idstr_compare(const void *v1, const void *v2)
/* case-insensitive comparison by 'idstr' from struct_id struct */
{
    const struct_id *p1, *p2;

    p1 = (const struct_id *) v1;
    p2 = (const struct_id *) v2;

    return case_strcmp(p1->idstr, p2->idstr);
}

int ntop_compare(const void *v1, const void *v2)
/* by b_tval (POSITIVE value) from struct_ntopArr struct */
{
    const struct_ntop *p1, *p2;

    p1 = (const struct_ntop *) v1;
    p2 = (const struct_ntop *) v2;

    /* in descending order */
    return (p2->b_tval > p1->b_tval) ? 1 : (p2->b_tval < p1->b_tval) ? -1 : 0;
}

void call_system(char *cmd)
{
    int retval;

    if (system(NULL)) {
        retval = system(cmd);
        if (retval)
            fatal("system('%s') returns nonzero (%d)\n", cmd, retval);
    } else
        fatal("Error: your OS does not support system() call\n");
}

long decompose_list_ids_to_file(char *ids_str, char *sep_chars, char *filename)
{
    char *p;
    long num_ids = 0;
    FILE *fp;

    fp = open_file(filename, "w");

    for (p = strtok(ids_str, sep_chars); p != NULL; p = strtok(NULL, sep_chars)) {
        num_ids++;
        fprintf(fp, "%s\n", trim(p));
    }

    close_file(fp);

    return num_ids;
}

long extract_numlist(char *num_str, char *sep_chars, long lmin, long lmax, long *vnum)
/* Extract a list of integers from a number string. The items must be
 * uniq, within range, and are separated by *sep_chars */
{
    char msg[BUF512], *items[BUF512];
    long i, k, num_items, num = 0;

    if (is_empty_string(num_str))
        return num;

    num_items = item_list(num_str, items, BUF512, sep_chars);

    for (i = 1; i <= num_items; i++) {
        k = cvt2long(items[i]);

        if (k == LONG_MAX) {
            sprintf(msg, "\t<%s> is not a valid integer", items[i]);
            log_msg(msg);

        } else if (k < lmin || k > lmax) {
            sprintf(msg, "\t<%ld> is out of range [%ld, %ld]", k, lmin, lmax);
            log_msg(msg);

        } else if (lval_in_array(k, 1, num, vnum)) {
            sprintf(msg, "\t<%ld> already exists", k);
            log_msg(msg);

        } else {
            num++;
            vnum[num] = k;
        }
    }

    if (num > 1)
        qsort(vnum + 1, num, sizeof(vnum[1]), lval_compare);

    return num;
}
