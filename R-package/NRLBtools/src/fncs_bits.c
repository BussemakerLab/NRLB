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

void convert_qua2seq(long *c, long num, char *seq)
/* decode quaternary number array (0123) to base sequence (ACGT) */
{
    long i;

    for (i = 0; i < num; i++)
        seq[num - i - 1] = BASES[c[i]];
    seq[i] = '\0';

    if (Gvars.RNA)
        cvtstr_c1toc2(seq, 'T', 'U');
}

void convert_dec2seq(long idx, long num, char *seq)
/* decode a decimal number to base sequence (ACGT) */
{
    long c[QUALEN16];

    fast_d2q(idx, c);
    convert_qua2seq(c, num, seq);
}

void fast_d2q(unsigned long n, long *c)
/* fast decimal to quaternary conversion */
{
    long i;

    for (i = 0; i < 2 * QUALEN16; i += 2)  /* 32 bits for 4-byte long */
        *(c++) = (n >> i) & 0x3;
}

unsigned long fast_q2d(long *c)
/* fast quaternary to decimal conversion */
{
    long i = QUALEN16;
    unsigned long n = 0;

    while (i--) {
        n <<= 2;
        n += *(c + i);
    }

    return n;
}

/* fast quaternary to decimal from an ACGT string */
long fast_str_q2d(char *p)
{
    char *pq;
    long i, idx = 0, num = strlen(p);

    pq = my_strdup(p);

    seq2number(pq);

    for (i = 0; i < num; i++) {
        idx <<= 2;
        idx += pq[i];
    }

    free_cvector(pq, 0, DUMMY);

    return idx;
}

void reverse_cmpl(char *seq, long debug)
{
    char *p, msg[BUF512];
    unsigned long i;

    reverse_string(seq);  /* in-place reverse */

    for (i = 0; i < strlen(seq); i++) {
        p = strchr(CASE_BASES, seq[i]);
        if (p)
            seq[i] = CASE_RC_BASES[p - CASE_BASES];

        else if (debug) {
            sprintf(msg, "\tnonstandard base '%c' ... [%s]", seq[i], CASE_BASES);
            log_msg(msg);
        }
    }
}

unsigned long fast_cidx(unsigned long n, long Xcount)
/* decimal to its reverse complementary decimal: compacted version */
{
    long i, j = 0, k = 2 * Xcount;
    unsigned long c, idx = 0;

    for (i = 0; i < Xcount; i++) {
        c = TOP_BIDX3 - ((n >> j) & 0x3);
        idx += (c << (k - j - 2));
        j += 2;
    }

    return idx;
}

unsigned long turnon_bits(long startbit, long numbits)
{
    unsigned long word;

    word = ~((unsigned long) ~0 << numbits) << startbit;

    return word;
}

unsigned long turnoff_bits(long startbit, long numbits)
{
    return ~turnon_bits(startbit, numbits);
}

unsigned long toggle_bits(unsigned long word, long startbit, long numbits)
{
    word ^= turnon_bits(startbit, numbits);

    return word;
}

long testbit(unsigned long word, long bit_to_test)
{
    word >>= bit_to_test;
    word &= 1;

    return word;
}

void printbits(unsigned long word)
{
    long i, USL_SIZE = 8 * sizeof(unsigned long);

    for (i = USL_SIZE - 1; i >= 0; i--)
        fprintf(stderr, "%1ld", testbit(word, i));

    fprintf(stderr, "\n");
}
