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

long get_psam_strand(long cmd_strand, long psam_strand)
/* for PSAM: command-line-setting; PSAM specification; +1 */
{
    long strand;

    strand = (cmd_strand == UNSET_LVAL) ? psam_strand : cmd_strand;

    if (!strand)
        strand = 1;

    return strand;
}

long get_motif_strand(long cmd_strand)
/* for motif: if unset from command-line, use leading strand */
{
    return (cmd_strand == UNSET_LVAL) ? 1 : cmd_strand;
}

void log_strand_msg(long strand)
{
    display_strand_msg(Gvars.RUNLOG, "PSAM based on", "", strand);
    display_strand_msg(Gvars.PRGLOG, "PSAM based on", "", strand);
}

void display_strand_msg(FILE * fp, char *prefix, char *suffix, long strand)
{
    char *stnd_msg;

    stnd_msg = strand2msg(strand);

    if (strand == 2)
        fprintf(fp, "%s %s strands%s\n", prefix, stnd_msg, suffix);
    else
        fprintf(fp, "%s %s strand%s\n", prefix, stnd_msg, suffix);

    free_cvector(stnd_msg, 0, DUMMY);
}

char *strand2msg(long strand)
{
    char msg[BUF512];

    if (strand == 1)
        strcpy(msg, "forward");
    else if (strand == 2)
        strcpy(msg, "both");
    else if (strand == -1)
        strcpy(msg, "reverse complementary");
    else if (strand == 0)
        strcpy(msg, "dynamically determined");
    else
        fatal("wrong value for strand <%d>: must be [-1, 0, 1, 2]\n", strand);

    return my_strdup(msg);
}

long msg2strand(char *msg)
{
    long strand;

    lowerstr(msg);

    if (strstr(msg, "lead") || strstr(msg, "forward"))
        strand = 1;
    else if (strstr(msg, "both") || strstr(msg, "double"))
        strand = 2;
    else if (strstr(msg, "reverse") || strstr(msg, "comp"))
        strand = -1;
    else  /* undetermined: see also get_psam_strand() */
        strand = 0;

    return strand;
}

void set_residual_file(char *resid_file, char *measfile, char *outdir)
/* derive residual file name: default to {measfile}.resid [NOT used]*/
{
    char msg[BUF512], str[BUF512], *p;

    if (is_empty_string(resid_file))  /* name unspecified from command line */
        sprintf(resid_file, "%s.resid", basename(measfile));

    else if (strchr(resid_file, '/') != NULL) {  /* remove directory */
        p = basename(resid_file);
        sprintf(msg, "\treset residual <%s> to <%s/%s>\n", resid_file, outdir, p);
        log_msg(msg);
        strcpy(str, p);  /* need a temporary variable: 'str' */
        strcpy(resid_file, str);
    }
}

void get_residuals(struct_data * tdat, long pnum, double **X, struct_fitpars * lfpars)
{
    long i, j, k;
    double dval;

    for (i = 1; i <= tdat->ncol; i++)
        for (j = 1; j <= tdat->nrow; j++) {
            dval = tdat->data[i][j];  /* initial value */

            if (dval < XBIG_CUTOFF) {
                dval -= lfpars[i].fval[0];  /* intercept, with x = 1.0 */
                for (k = 1; k <= pnum; k++)
                    dval -= lfpars[i].fval[k] * X[k][j];
                tdat->resid[i][j] = dval;
            }
        }
}

void get_predicted(struct_data * tdat)
/* calculated as 'data[] - resid[]' and replaced resid[] */
{
    long i, j;
    double dval;

    if (!tdat->resid)
        fatal("no residual info yet -- can't run get_predicted()\n");

    for (i = 1; i <= tdat->ncol; i++)
        for (j = 1; j <= tdat->nrow; j++) {
            dval = tdat->data[i][j];
            if (dval < XBIG_CUTOFF)
                tdat->resid[i][j] = dval - tdat->resid[i][j];
        }
}

double get_slope_abs_tvalue(long n, double sX, double sX2, double sXY, double sY,
                            double sY2)
/* calculate the abs(t-value) based on FIVE sums linear regression.
 * note: r^2 will gives parallel results with ONE experiment */
{
    double SSX, SSY, SSXY, SSE, aveX, aveY, b, b_se;
    long df = n - 2;

    if (df <= 0)
        fatal("df = %ld: degree of freedom must be positive\n", df);

    aveX = sX / n;
    aveY = sY / n;

    SSX = sX2 - aveX * sX;
    if (SSX <= TOLX)  /* special case: no variation in X, as in Ron's case */
        return 0.0;

    SSY = sY2 - aveY * sY;
    SSXY = sXY - aveX * sY;

/*
  r^2 and t-value are not strictly parallel with multiple experiments;
  as in the spellman-alpha data set distributed with the REDUCE Suite.
  for example:

  34,35c34,35 with r^2
  2  ACGCG    +14.077   0.03470  1008[849/5514]  forward     alpha_factor_release_sample016[4]
  3  ACGCGT   +14.099   0.03435  347[309/5591] palindrome  alpha_factor_release_sample024[12]
  --------- with t-value
  2  ACGCGT   +14.099   0.03435  347[309/5591] palindrome  alpha_factor_release_sample024[12]
  3  ACGCG    +14.077   0.03470  1008[849/5514]  forward     alpha_factor_release_sample016[4]
*/
    /* return (SSXY * SSXY) / (SSX * SSY); */

    b = SSXY / SSX;
    if (Gvars.seed_criterion == SEED_BY_SLOPE_ITSELF)
        return fabs(b);

    SSE = SSY - b * SSXY;
    b_se = sqrt(SSE / df / SSX);

    return (b_se > TOLX) ? (fabs(b) / b_se) : 0.0;
}

void lsfit_parameters(long n, double sX, double sX2, double sXY, double sY, double sY2,
                      struct_ls5fit * fitpars)
/* linear fitted parameters for [y = a + bx] are stored in struct fitpars */
{
    double SSX, SSY, SSXY, SSR, SSE, aveX, aveY;
    double a, b, s2, a_se, b_se, a_tval, b_tval;
    long df = n - 2;

    if (df <= 0)
        fatal("df = %ld: degree of freedom must be positive\n", df);

    fitpars->num = n;
    fitpars->sumX = sX;

    aveX = sX / n;
    aveY = sY / n;

    SSX = sX2 - aveX * sX;
    if (SSX <= TOLX)
        return;

    SSY = sY2 - aveY * sY;
    fitpars->SSY = SSY;

    SSXY = sXY - aveX * sY;

    b = SSXY / SSX;
    fitpars->b = b;  /* slope */

    a = aveY - b * aveX;
    fitpars->a = a;  /* intercept */

    SSR = b * SSXY;
    fitpars->SSR = SSR;
    fitpars->r2 = SSR / SSY;  /* also: (SSXY * SSXY) / (SSX * SSY) */

    SSE = SSY - SSR;
    fitpars->SSE = SSE;

    s2 = SSE / df;
    if (s2 > TOLX)
        fitpars->F = SSR / s2;  /* F value */

    b_se = sqrt(s2 / SSX);  /* std. error for slope */
    if (b_se > TOLX) {
        b_tval = b / b_se;
        fitpars->b_tval = b_tval;  /* t-value for slope: with the same sign as slope */
        fitpars->b_pval = NRC_betai(0.5 * df, 0.5, df / (df + b_tval * b_tval));
    }

    a_se = b_se * sqrt(sX2 / n);  /* std. error for intercept */
    if (a_se > TOLX) {
        a_tval = a / a_se;
        fitpars->a_tval = a_tval;  /* t-value for the intercept */
        fitpars->a_pval = NRC_betai(0.5 * df, 0.5, df / (df + a_tval * a_tval));
    }
}

static long get_y_vs_x_5sums(long nrow, double *x, double *y, double *five_sums)
{
    long i, num = 0;

    for (i = 1; i <= 5; i++)
        five_sums[i] = 0.0;

    for (i = 1; i <= nrow; i++) {
        if (y[i] < XBIG_CUTOFF) {
            num++;
            five_sums[1] += x[i];
            five_sums[2] += x[i] * x[i];
            five_sums[3] += x[i] * y[i];
            five_sums[4] += y[i];
            five_sums[5] += y[i] * y[i];
        }
    }

    return num;
}

double get_y_vs_x_slope_abs_tvalue(long nrow, double *x, double *y)
{
    long num;
    double tval, fsums[6];

    num = get_y_vs_x_5sums(nrow, x, y, fsums);
    tval = get_slope_abs_tvalue(num, fsums[1], fsums[2], fsums[3], fsums[4], fsums[5]);

    return tval;
}

void get_y_vs_x_lsfitpars(long nrow, double *x, double *y, struct_ls5fit * fitpars)
/* y[] vs x[] linear ls-fitting: y = a + b * x */
{
    long num;
    double fsums[6];

    num = get_y_vs_x_5sums(nrow, x, y, fsums);
    init_fitpars(fitpars);
    lsfit_parameters(num, fsums[1], fsums[2], fsums[3], fsums[4], fsums[5], fitpars);
}

long get_num_mids(long num, double *counts, double *y)
/* number of matched ids: one per sequence maximum */
{
    long i, num_mids = 0;

    for (i = 1; i <= num; i++)
        if (y[i] < XBIG_CUTOFF && counts[i] > TOLX)
            num_mids++;

    return num_mids;
}

void print_fitpars(char *str, struct_ls5fit * fitpars)
{
    char msg[BUF512];

    sprintf(msg, "%s linear fit statistics:", str);
    log_msg(msg);

    output_linefit_pars(Gvars.RUNLOG, fitpars);
    output_linefit_pars(Gvars.PRGLOG, fitpars);
}

void get_psam_fitpars(struct_data * tdat, struct_seqArr * seqs, struct_seed * seed,
                      struct_ls5fit * fitpars, long debug)
{
    double *y, *counts;

    counts = dvector(1, seqs->num_seqs);
    y = tdat->resid[seed->col_idx];

    get_psam_counts(seqs, seed->topo, seed->strand, seed->w, counts);
    get_y_vs_x_lsfitpars(tdat->nrow, counts, y, fitpars);

    seed->num_mids = get_num_mids(tdat->nrow, counts, y);

    free_dvector(counts, 1, DUMMY);

    if (debug)
        print_fitpars("PSAM", fitpars);
}

/* if t->Xcount = 0 then return 1 (single strand) or 2 (both strands) */
double get_count_per_window(char *p0, struct_topo * t, long strand, double *w)
{
    double Nw1 = 1.0, Nw2 = 1.0, Nw = 1.0;
    long i, eidx = t->Tcount - 1;  /* end-index */

    if (strand == 1) {
        for (i = 0; i < t->Xcount; i++)
            Nw *= w[i * NUM_BASE4 + p0[t->offset[i]]];

    } else if (strand == -1) {
        for (i = 0; i < t->Xcount; i++)
            Nw *= w[i * NUM_BASE4 + TOP_BIDX3 - p0[eidx - t->offset[i]]];

    } else if (strand == 2) {
        for (i = 0; i < t->Xcount; i++) {
            Nw1 *= w[i * NUM_BASE4 + p0[t->offset[i]]];
            Nw2 *= w[i * NUM_BASE4 + TOP_BIDX3 - p0[eidx - t->offset[i]]];
        }
        Nw = Nw1 + Nw2;

    } else
        fatal("strand=%ld not in [+1, -1, 2]\n", strand);

    return Nw;
}

double get_max_pcount_seq(struct_seq * seq, struct_topo * t, long strand, double *w)
{
    char *p;  /* as a short-hand notation */
    long m, nb, nseq;
    double Nw, pmax = -XBIG;

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {  /* loop over each fragment */
        nb = seq->num_nb[nseq];  /* number of bases per fragment; */
        p = seq->seq_fragments[nseq];  /* pointer to starting base */

        for (m = 0; m <= nb - t->Tcount; m++) {  /* per window */
            if (t->num_parens && !with_matched_parens_bin(t, p + m))
                continue;

            Nw = get_count_per_window(p + m, t, strand, w);
            pmax = dval_max(pmax, Nw);
        }
    }

    return pmax;  /* max. per sliding window */
}

double csum_per_seq(struct_seq * seq, struct_topo * t, long strand, double *w,
                    long normalize, double pmax, double threshold)
{
    char *p;  /* as a short-hand notation */
    long m, nb, nseq;
    double Nw, Nsum = 0.0;

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {  /* loop over each fragment */
        nb = seq->num_nb[nseq];  /* number of bases per fragment; */
        p = seq->seq_fragments[nseq];  /* pointer to starting base */

        for (m = 0; m <= nb - t->Tcount; m++) {  /* per window */
            if (t->num_parens && !with_matched_parens_bin(t, p + m))
                continue;

            Nw = get_count_per_window(p + m, t, strand, w);

            /* check for normalization first; also check 'pmax' to
             * avoid pmax = 0.0 in case of no match (same for below) */
            if (normalize && pmax)
                Nw /= pmax;

            if (Nw >= threshold)  /* add only those above threshold */
                Nsum += Nw;
        }
    }

    return Nsum;  /* per sequence */
}

void count_per_sequence(struct_seq * seq, struct_topo * t, long strand, double *w,
                        long normalize, double pmax, double threshold, double Nsum,
                        long column_wise, FILE * fp)
{
    char *p;
    double Nw;
    long k, m, nb, nseq, idx0;

    if (column_wise) {
        fprintf(fp, "### %s\t%g\n", seq->id, Nsum);  /* total per sequence */
        if (threshold > 1.0)  /* only for the cumulative */
            return;

    } else {
        fprintf(fp, "%s\t%g", seq->id, Nsum);
        if (threshold > 1.0) {  /* only for the cumulative */
            fprintf(fp, "\n");
            return;
        }
    }

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {
        nb = seq->num_nb[nseq];  /* number of bases per fragment; */
        p = seq->seq_fragments[nseq];  /* pointer to starting base */
        idx0 = p - seq->seq;  /* offset of the current fragment */

        for (m = 0; m <= nb - t->Tcount; m++) {  /* per window */
            if (t->num_parens && !with_matched_parens_bin(t, p + m))
                continue;

            Nw = get_count_per_window(p + m, t, strand, w);

            /* check for normalization first; also check 'pmax' to
             * avoid pmax = 0.0 in case of no match (see also above) */
            if (normalize && pmax)
                Nw /= pmax;

            if (Nw >= threshold && Nw > 0) {
                k = idx0 + m + 1;  /* 1-based index */
                if (column_wise)
                    fprintf(fp, "%ld\t%g\n", k, Nw);
                else
                    fprintf(fp, "\t%ld:%g", k, Nw);
            }
        }
    }

    if (!column_wise)
        fprintf(fp, "\n");
}

double get_pcount_seq(struct_seq * seq, struct_topo * t, long strand, double *w)
{
    char *p;  /* as a short-hand notation */
    long i, m, nb, nseq, eidx;
    double W, N = 0.0;

    eidx = t->Tcount - 1;  /* end-index */

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {  /* loop over each fragment */
        nb = seq->num_nb[nseq];  /* number of bases per fragment; */
        p = seq->seq_fragments[nseq];  /* pointer to starting base */

        if (strand == 1 || strand == 2) {  /* leading strand or both strands */
            for (m = 0; m <= nb - t->Tcount; m++) {  /* per window */
                if (t->num_parens && !with_matched_parens_bin(t, p + m))
                    continue;

                W = 1.0;
                for (i = 0; i < t->Xcount; i++)
                    W *= w[i * NUM_BASE4 + p[m + t->offset[i]]];
                N += W;
            }
        }

        if (strand == -1 || strand == 2) {  /* reverse complementary strand or both */
            for (m = 0; m <= nb - t->Tcount; m++) {
                if (t->num_parens && !with_matched_parens_bin(t, p + m))
                    continue;  /* to be checked -- maybe not relevant */

                W = 1.0;
                for (i = 0; i < t->Xcount; i++)
                    W *= w[i * NUM_BASE4 + TOP_BIDX3 - p[m + eidx - t->offset[i]]];
                N += W;
            }
        }
    }

    return N;
}

void get_psam_counts(struct_seqArr * seqs, struct_topo * t, long strand, double *w,
                     double *counts)
{
    long i;

    for (i = 1; i <= seqs->num_seqs; i++)
        counts[i] = get_pcount_seq(&seqs->sequences[i], t, strand, w);
}

/* get integer counts of 'motif_seq' within all the sequences */
void get_whole_counts(char *motif_seq, long strand, struct_seqArr * seqs,
                      double *motif_count)
{
    char *p;
    long i, j, k, nseq, nb, N, eidx, Xcount = strlen(motif_seq);
    double *motif_mtx;
    struct_seq *s;

    motif_mtx = dvector(0, 4 * Xcount - 1);

    iupac2psam(motif_seq, motif_mtx);
    eidx = Xcount - 1;

    for (i = 1; i <= seqs->num_seqs; i++) {
        N = 0;  /* initialize it */
        s = &seqs->sequences[i];
        for (nseq = 1; nseq <= s->num_fragments; nseq++) {
            nb = s->num_nb[nseq];
            p = s->seq_fragments[nseq];
            for (j = 0; j <= nb - Xcount; j++) {
                if (strand == 1 || strand == 2) {  /* leading strand or both strands */
                    for (k = 0; k < Xcount; k++)
                        if (!motif_mtx[k * NUM_BASE4 + p[j + k]])
                            break;
                    if (k == Xcount)
                        N++;
                }

                if (strand == -1 || strand == 2) {  /* reverse cmpl strand or both */
                    for (k = 0; k < Xcount; k++)
                        if (!motif_mtx[k * NUM_BASE4 + TOP_BIDX3 - p[j + eidx - k]])
                            break;
                    if (k == Xcount)
                        N++;
                }
            }
        }

        motif_count[i] = N;
    }

    free_dvector(motif_mtx, 0, DUMMY);
}

long svd_check_psam_counts(long ndat, long npar, char *type, double **X)
/* here X is 'psam_counts' which needs to be transposed */
{
    long i, k = 0;
    double **u, **v, *w;

    u = dmatrix(1, ndat, 1, npar);
    v = dmatrix(1, npar, 1, npar);
    w = dvector(1, npar);

    transpose_matrix(X, npar, ndat, u);

    NRC_svdcmp(u, ndat, npar, w, v);
    NRC_svdzwi(w, npar);

    for (i = 1; i <= npar; i++)
        if (w[i] == 0.0)
            k++;

    free_dmatrix(u, 1, DUMMY, 1, DUMMY);
    free_dmatrix(v, 1, DUMMY, 1, DUMMY);
    free_dvector(w, 1, DUMMY);

    if (k) {
        char msg[BUF512];

        sprintf(msg, "\n%ld degeneracy within %ld %s(s)\n", k, npar, type);
        log_msg(msg);
    }

    return k;
}

void write_fit_yX(long icol, long nrow, long ncol, double *y, double **X, long chk)
{
    long i, j;
    char filename[BUF512];
    FILE *fp;

    if (!chk)
        return;

    sprintf(filename, "fit_yX_%3.3ld.dat", icol);
    fp = open_file(filename, "w");

    fprintf(fp, "## Note: the 1st column is y; following %ld columns are X\n", ncol);
    for (i = 1; i <= nrow; i++) {
        fprintf(fp, "%+g", y[i]);
        for (j = 1; j <= ncol; j++)  /* 1st column in X is all 1s for intercept */
            fprintf(fp, "\t%+g", X[i][j]);
        fprintf(fp, "\n");
    }

    close_file(fp);
}

long get_number_of_matched_ids(long num, double *counts)
{
    long i, k = 0;

    for (i = 1; i <= num; i++)
        if (counts[i] > 0.0)
            k++;

    return k;
}

/* Output model summary from a MotifREDUCE run, following REDUCE.
   Note: multiple linear fit parameters for each experiment are
   directly available in 'lfpars' */
void output_model_summary(args_reduce * args, long num_motifs, struct_seed * smotifs,
                          struct_data * tdat, struct_seqArr * seqs,
                          struct_fitpars * lfpars)
{
    char *motif, filename[BUF512], *model_file = "model_summary.txt";
    double *counts;
    long i, j, k, num_seqs = seqs->num_seqs;
    FILE *fp;

    struct_ls5fit fitpars;
    struct_seed *seed;

    sprintf(filename, "%s/%s", args->outdir, model_file);
    fp = open_file(filename, "w");

    fprintf(fp, "Number_of_motifs=%ld\tNumber_of_experiments=%ld\n\n", num_motifs,
            tdat->ncol);

    counts = dvector(1, num_seqs);

    for (i = 1; i <= tdat->ncol; i++) {
        fprintf(fp, "experiment=%s\tcolumn=%ld\ttotal_valid_ids=%ld\n\n",
                tdat->col_names[i], i, tdat->oknum[i]);
        fprintf(fp, " no.  motif\tr2(univariate)\tF(univariate)\tF(multivariate)\tmatches"
                "\tmatched_ids\n");
        print_sep(fp, '-', 83);

        for (j = 1; j <= num_motifs; j++) {
            seed = &smotifs[j];
            motif = seed->full;

            get_motif_counts(motif, seed->strand, seqs, counts);

            k = get_number_of_matched_ids(num_seqs, counts);
            /* Note: smotifs[j].num_mids is experiment-specific; thus
             * not work with multiple experiments */

            get_y_vs_x_lsfitpars(tdat->nrow, counts, tdat->data[i], &fitpars);

            fprintf(fp, "%3ld   %s\t%.6f\t%+.6f\t%+.6f\t  %ld\t  %ld\n", j,
                    smotifs[j].full, fitpars.r2, fitpars.b, lfpars[i].fval[j],
                    (long) fitpars.sumX, k);

            init_dvector(counts, 1, num_seqs, 0.0);
        }

        print_sep(fp, '-', 83);
    }

    free_dvector(counts, 1, DUMMY);

    close_file(fp);
}

void perform_multifit(struct_data * tdat, long na, double **psam_counts,
                      struct_fitpars * lfpars, long bonf, long chk)
{
    char msg[BUF512];
    long df, i, j, k, num, ma = na + 1;
    double t, dval, p_val, sigma, chisq, *a, **X, *y, **covar;

    a = dvector(1, ma);  /* number of parameters to be fitted */
    covar = dmatrix(1, ma, 1, ma);

    for (i = 1; i <= tdat->ncol; i++) {
        num = tdat->oknum[i];  /* # of valid entries: w/ seq & non-NaN */

        X = dmatrix(1, num, 1, ma);
        y = dvector(1, num);

        /* populate matrix X and array y for ls-fitting */
        num = 0;
        for (j = 1; j <= tdat->nrow; j++) {
            dval = tdat->data[i][j];

            if (dval < XBIG_CUTOFF) {
                num++;
                y[num] = dval;
                X[num][1] = 1.0;  /* for intercept */
                for (k = 1; k <= na; k++)
                    X[num][k + 1] = psam_counts[k][j];
            }
        }

        dval = var_dvector(y, 1, num);
        if (dval <= TOLX) {
            sprintf(msg, "NO VARIATION in y [%g]: are you sure?", dval);
            log_msg(msg);
        }

        write_fit_yX(i, num, ma, y, X, chk);

        /* NRC_lfit(X, y, num, ma, a, covar, &chisq); */
        NRC_my_svdfit(X, y, num, ma, a, covar, &chisq);

        lfpars[i].num = num;
        lfpars[i].ma = ma;  /* number of fitted parameters including the intercept */
        lfpars[i].bonf = bonf;
        lfpars[i].chisq = chisq;

        df = num - ma;  /* degree of freedom for t-distribution */
        sigma = sqrt(chisq / df);
        for (j = 1; j <= ma; j++) {
            t = a[j] / sqrt(covar[j][j]) / sigma;
            p_val = NRC_betai(0.5 * df, 0.5, df / (df + t * t));

            k = j - 1;  /* starting from 0 for the intercept */
            lfpars[i].fval[k] = a[j];
            lfpars[i].tval[k] = t;

            /* lfpars[i].pval[k] = 1.0 - pow(1.0 - p_val, bonf); */
            lfpars[i].pval[k] = dval_min(p_val * bonf, 1.0);
        }

        free_dmatrix(X, 1, DUMMY, 1, DUMMY);
        free_dvector(y, 1, DUMMY);
    }

    free_dvector(a, 1, DUMMY);
    free_dmatrix(covar, 1, DUMMY, 1, DUMMY);
}

void process_sig_seed(long pnum, struct_data * tdat, struct_seqArr * seqs,
                      struct_seed * seed, double **psam_counts, struct_fitpars * lfpars)
{
    output_html_spar_fitpars(pnum, tdat, seqs, seed);

    get_psam_counts(seqs, seed->topo, seed->strand, seed->w, psam_counts[pnum]);

    perform_multifit(tdat, pnum, psam_counts, lfpars, seed->bonf, FALSE);

    get_residuals(tdat, pnum, psam_counts, lfpars);
}

long psam_motif_sig_output(char *type, double E_value, double cutoff)
{
    char msg[BUF512];
    long isok = E_value < cutoff;

    sprintf(msg, "\tE-value=%g", E_value);
    log_msg(msg);

    if (isok) {
        sprintf(msg, "\tThis %s is significant (E-value smaller than specific "
                "cutoff of %g)", type, cutoff);
        log_msg(msg);

    } else {
        sprintf(msg, "\tThis %s is NOT significant (E-value larger than specific "
                "cutoff of %g)", type, cutoff);
        log_msg(msg);
    }

    return isok;
}

long is_sig_motif(struct_seed * seed, double p_value)
{
    char msg[BUF512];
    double pval, b_pval = seed->fitpars.b_pval;

    /* the following way of correcting p-value leads to 0 for small b_pval */
    /* pval = 1.0 - pow(1.0 - b_pval, seed->bonf); */
    pval = dval_min(b_pval * seed->bonf, 1.0);

    log_msg("Checking seed motif significance:");

    sprintf(msg, "\tnumber of tested candidate motifs: %ld", seed->bonf);
    log_msg(msg);

    return (seed->fitpars.SSY > Gvars.misc.SSY) &&
        psam_motif_sig_output("motif", pval, p_value);
}

void free_unsig_seed(double *counts, struct_seed * seed)
{
    free_dvector(counts, 1, DUMMY);
    free_seed(seed);
}

long isok_motif(double p_value, struct_data * tdat, struct_seqArr * seqs,
                struct_seed * seed, long *psam_num, double **psam_counts,
                struct_fitpars * lfpars)
{
    if (is_sig_motif(seed, p_value)) {
        process_sig_seed(*psam_num, tdat, seqs, seed, psam_counts, lfpars);
        return TRUE;
    }

    /* this motif is not significant */
    free_unsig_seed(psam_counts[*psam_num], seed);
    (*psam_num)--;

    return FALSE;
}

long is_sig_ron_correction(double N, long Xcount, double abs_r, double p_value)
/* using Ron's two empirical equations to decide if 'abs_r' is significant */
{
    char msg[BUF512];
    long df = N - 2;  /* based on t-distribution */
    double sqrt_N = sqrt(N), r0, sigma, t_val, pval_onesided, pval;

    r0 = (1.63 + 0.58 * Xcount) / sqrt_N;
    sigma = 0.66 / sqrt_N;

    t_val = (abs_r - r0) / sigma;

    /* 'pval_onesided' corresponds to R: pt(abs(t_val), df, lower.tail = FALSE) */
    pval_onesided = 0.5 * NRC_betai(0.5 * df, 0.5, df / (df + t_val * t_val));  /* t-distribution */
    pval = (t_val > 0) ? pval_onesided : (1.0 - pval_onesided);

    log_msg("Checking PSAM significance:");

    sprintf(msg, "\t|r|=%g\tr0=%g\tsigma=%g\tt_value=%g", abs_r, r0, sigma, t_val);
    log_msg(msg);

    return psam_motif_sig_output("PSAM", pval, p_value);
}

long isok_psam(double p_value, struct_data * tdat, struct_seqArr * seqs,
               struct_seed * seed, long *psam_num, double **psam_counts,
               struct_fitpars * lfpars)
{
    double abs_r;
    struct_ls5fit fitpars;

    get_psam_fitpars(tdat, seqs, seed, &fitpars, TRUE);
    abs_r = sqrt(fitpars.r2);

    if ((fitpars.SSY > Gvars.misc.SSY) &&
        is_sig_ron_correction(seed->fitpars.num, seed->topo->Xcount, abs_r, p_value)) {
        process_sig_seed(*psam_num, tdat, seqs, seed, psam_counts, lfpars);
        return TRUE;
    }

    /* this PSAM is not significant */
    free_unsig_seed(psam_counts[*psam_num], seed);
    (*psam_num)--;

    return FALSE;
}

void reset_Xmtx(long num_ok, long num0, long *num, long *is_okay, long num_seqs,
                double **counts, struct_tag * str_tags)
/* Note: the memory space for the redundant ones are not freed here. */
{
    long i = 0, j;

    if (num_ok == 0) {
        log_msg("no matchs!");
        fatal("\n");
    }

    if (num_ok == num0)
        return;

    while (++i <= num_ok) {
        if (is_okay[i])
            continue;

        for (j = i + 1; j <= num0; j++) {
            if (!is_okay[j])  /* skip bad one */
                continue;

            copy_dvector(counts[j], 1, num_seqs, counts[i]);
            copy_strtag(&str_tags[j], &str_tags[i]);

            is_okay[j] = 0;

            break;  /* one row per iteration */
        }
    }

    *num = num_ok;
}

void check_for_zeros(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx)
/* handle cases where Xmtx[i] are all zeros, e.g., for a very long motif */
{
    char msg[BUF512];
    long i, j, num_ok = 0, num0 = *num, *is_okay;

    if (num0 < 1)
        return;

    is_okay = lvector(1, num0);
    init_lvector(is_okay, 1, num0, 1);  /* assuming all are Okay */

    for (i = 1; i <= num0; i++) {
        for (j = 1; j <= num_Xrow; j++)
            if (Xmtx[i][j])
                break;

        if (j > num_Xrow) {  /* all are zeros */
            is_okay[i] = 0;  /* as a marker */
            sprintf(msg, "\tskip %ld %s: NO match in sequence", i, str_tags[i].str);
            log_msg(msg);

        } else
            num_ok++;
    }

    reset_Xmtx(num_ok, num0, num, is_okay, num_Xrow, Xmtx, str_tags);

    free_lvector(is_okay, 1, DUMMY);
}

void check_for_degeneracy(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx)
/* handle degenerate cases such as A=T & C=G for single base Xmtx with both strands */
{
    char msg[BUF512];
    long i, j, num_ok = 0, num0 = *num, nmsg = 100, *is_okay;
    time_t time0;

    struct_ls5fit fitpars;

    if (num0 < 2)
        return;

    if (num0 >= nmsg) {
        sprintf(msg, "\npairwise checking for degeneracy (%ldx%ld)", num0, num0);
        log_msg(msg);
        time(&time0);
    }

    is_okay = lvector(1, num0);
    init_lvector(is_okay, 1, num0, 1);

    for (i = 1; i < num0; i++) {
        if (!is_okay[i])
            continue;

        for (j = i + 1; j <= num0; j++) {
            if (!is_okay[j])
                continue;

            get_y_vs_x_lsfitpars(num_Xrow, Xmtx[i], Xmtx[j], &fitpars);

            if (fitpars.SSE <= TOLX) {
                sprintf(msg, "\tignore %s[%ld] ==> redundant of %s[%ld]",
                        str_tags[j].str, j, str_tags[i].str, i);
                log_msg(msg);
                is_okay[j] = 0;  /* set a marker in later rows */
            }
        }
    }

    for (i = 1; i <= num0; i++)
        if (is_okay[i])
            num_ok++;

    reset_Xmtx(num_ok, num0, num, is_okay, num_Xrow, Xmtx, str_tags);

    free_lvector(is_okay, 1, DUMMY);

    if (num0 >= nmsg)
        print_used_time(time0);
}

void check_for_multifit(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx)
/* handle degenerate cases such as single nucleotide fit A/C/G/T with equal
 * sequence length: only three are independent. This function is called
 * recursively to enumerate all such possibilities. */
{
    char msg[BUF512];
    long i, j, k, all_is_fine = TRUE, num_valid_rows, num0 = *num, *is_okay;
    double chisq, *a, **X, *y, **covar;

    if (num0 < 2)
        return;

    is_okay = lvector(1, num0);
    init_lvector(is_okay, 1, num0, 1);

    a = dvector(1, num0);
    covar = dmatrix(1, num0, 1, num0);

    num_valid_rows = get_number_of_valid_rows(num0, num_Xrow, Xmtx);

    X = dmatrix(1, num_valid_rows, 1, num0);
    y = dvector(1, num_valid_rows);

    for (i = 1; i <= num_valid_rows; i++)
        X[i][1] = 1.0;  /* column for intercept */

    k = num0;  /* loop over each combination */
    while (k > 1) {
        num_valid_rows = 0;  /* reset the counter */

        for (i = 1; i <= num_Xrow; i++) {
            for (j = 1; j <= num0; j++)
                if (Xmtx[j][i] > XBIG_CUTOFF)
                    break;

            if (j <= num0)  /* with NA item in i-th row */
                continue;

            num_valid_rows++;
            y[num_valid_rows] = Xmtx[num0][i];  /* the last one against all previous ones */
            for (j = 1; j < num0; j++)
                X[num_valid_rows][j + 1] = Xmtx[j][i];  /* note the index order */
        }

        /* NRC_lfit(X, y, num_valid_rows, num0, a, covar, &chisq); */
        NRC_my_svdfit(X, y, num_valid_rows, num0, a, covar, &chisq);

        if (chisq > TOLX) {
            k--;
            swap_strtag(&str_tags[num0], &str_tags[k]);
            dvec_swap(&Xmtx[num0], &Xmtx[k]);

        } else {
            sprintf(msg, "\tignore %s[%ld] ==> redundant per linear combination of"
                    " the rest", str_tags[num0].str, k);
            log_msg(msg);

            sprintf(msg, "\t\t%ld\t%ld\t%g\n", k, num_valid_rows, chisq);
            log_msg(msg);

            all_is_fine = 0;
            is_okay[num0] = 0;
            reset_Xmtx(num0 - 1, num0, num, is_okay, num_Xrow, Xmtx, str_tags);

            break;
        }
    }

    free_dvector(a, 1, DUMMY);
    free_dmatrix(covar, 1, DUMMY, 1, DUMMY);
    free_dmatrix(X, 1, DUMMY, 1, DUMMY);
    free_dvector(y, 1, DUMMY);

    free_lvector(is_okay, 1, DUMMY);

    if (!all_is_fine)  /* repeat elimination */
        check_for_multifit(num, str_tags, num_Xrow, Xmtx);

    else {  /* put into original order */
        for (i = num0; i > 1; i--) {
            swap_strtag(&str_tags[i], &str_tags[i - 1]);
            dvec_swap(&Xmtx[i], &Xmtx[i - 1]);
        }
    }
}

void write_design_matrix(struct_data * tdat, long num, double **psam_counts,
                         struct_tag * str_tags, char *type, char *outdir)
/* write out ONLY the design matrix for verification purpose */
{
    char filename[BUF512];
    long i, j, idx;
    FILE *fp;

    sprintf(filename, "%s/design_matrix.tsv", outdir);
    fp = open_file(filename, "w");

    fprintf(fp, "## there are %ld data columns\n", num);
    for (i = 1; i <= num; i++)
        fprintf(fp, "## %s%ld -- %s\n", type, i, str_tags[i].str);

    for (i = 1; i <= num; i++)
        fprintf(fp, "\t%s%ld", type, i);
    fprintf(fp, "\n");

    for (i = 1; i <= tdat->nrow; i++) {
        idx = tdat->seqidx[i];

        if (idx) {  /* with matching fitmtx ID (as in sequence) */
            fprintf(fp, "%s", tdat->ids[i]);

            for (j = 1; j <= num; j++) {
                if (psam_counts[j][idx] > XBIG_CUTOFF)
                    fatal("NA for <%s> in fitmtx\n", tdat->ids[i]);
                fprintf(fp, "\t%g", psam_counts[j][idx]);
            }
            fprintf(fp, "\n");
        }
    }

    close_file(fp);
}

void calculate_fitpars(struct_data * tdat, struct_ls5fit * fitpars)
{
    long i, nrow = tdat->nrow;
    double *fitval;

    fitval = dvector(1, nrow);

    for (i = 1; i <= tdat->ncol; i++) {
        diff_dvector(fitval, tdat->data[i], tdat->resid[i], 1, nrow);
        get_y_vs_x_lsfitpars(nrow, fitval, tdat->data[i], &fitpars[i]);
    }

    free_dvector(fitval, 1, DUMMY);
}

void validate_motif(long *num, struct_tag * str_tags)
{
    char str0[BUF512], full[BUF512], tag[BUF512];
    long i, num_ok = 0;

    for (i = 1; i <= *num; i++) {
        strcpy(str0, str_tags[i].str);

        expand_char_num_to_full(str0, full);
        validate_expanded_pattern(full);

        if (!is_empty_string(full)) {
            num_ok++;
            strcpy(tag, str_tags[i].tag);
            free_p_strtag(&str_tags[num_ok]);
            init_strtag_with_strings(full, tag, &str_tags[num_ok]);
        }
    }

    get_unique_strtag(&num_ok, str_tags);

    if (num_ok == 0)
        fatal("NO valid motifs with the -motif option\n");

    *num = num_ok;
}

void cat_dir_fname(char *outdir, char *fname, char *fullname)
{
    if (str_pmatch(fname, "/") || str_pmatch(fname, "~/") || str_pmatch(fname, "./")
        || str_pmatch(fname, "../") || is_std_out_err(basename(fname)))
        strcpy(fullname, fname);  /* fname with absolute address or stdout/stderr */

    else
        sprintf(fullname, "%s/%s", outdir, basename(fname));
}

void log_prg(char *msg)
{
    fprintf(Gvars.PRGLOG, "%s\n", msg);
}

void log_run(char *msg)
{
    fprintf(Gvars.RUNLOG, "%s\n", msg);
}

void log_msg(char *msg)
{
    fprintf(Gvars.RUNLOG, "%s\n", msg);
    fprintf(Gvars.PRGLOG, "%s\n", msg);
}

void log_msg_vchk(char *msg)
{
    if (Gvars.VERBOSITY == VERBOSE)
        fprintf(Gvars.RUNLOG, "%s\n", msg);
    fprintf(Gvars.PRGLOG, "%s\n", msg);
}

void log_msg_exit(char *msg)
{
    log_msg(msg);
    exit(1);
}

void get_topo_or_dict(char *outdir, char *dictfile, char *topo_one, char *topo_list,
                      struct_topoArr * topo, struct_dic * mlist)
{
    init_mlist(mlist);  /* always initialize it for later free_* */

    if (is_empty_string(dictfile)) {
        if (!is_empty_string(topo_list))
            read_topos(topo_list, TRUE, topo);
        else
            read_topos(topo_one, FALSE, topo);
        write_topos(topo, outdir, "clean_topo.dat");

    } else {
        read_dictfile(dictfile, mlist);
        categorize_mlist(topo, mlist);
    }
}

void find_seed(char *dictfile, struct_dic * mlist, struct_topoArr * topo,
               struct_data * tdat, struct_seqArr * seqs, long ntop, long nsel,
               struct_seed * seed)
{
    if (is_empty_string(dictfile))
        find_seed_by_topos(tdat, seqs, topo, seed, ntop, nsel);

    else
        find_seed_by_dictfile(tdat, seqs, topo, mlist, seed, ntop, nsel);
}

void set_reduce_defaults(args_reduce * args)
{
    strcpy(args->seqfile, "");
    strcpy(args->measfile, "");
    strcpy(args->outdir, ".");
    strcpy(args->runlog, "stderr");

    strcpy(args->topo, "");
    strcpy(args->topo_list, "");
    strcpy(args->dictfile, "");

    args->p_value = P_VALUE;

    args->max_motif = MAX_MOTIF;
    args->strand = 1;  /* default to leading strand */

    args->iupac_pos = 0;  /* default to not check for IUPAC degeneracy */
    strcpy(args->iupac_sym, "ALL");

    args->ntop = NTOP;
    args->nsel = FALSE;  /* not manually picked up, but use the top seed */
}

void write_reduce_options(args_reduce * args)
{
    char filename[BUF512], options[BUF512];
    FILE *fp;

    bname_ext(Gvars.PROGNAME, "opt", options);

    sprintf(filename, "%s/%s", args->outdir, options);  /* to get the file name */

    fp = open_file(filename, "w");

    fprintf(fp, "## %s\n\n", getenv("PWD"));  /* current working directory */

    fprintf(fp, "%s \\\n", Gvars.PROGNAME);
    fprintf(fp, "\t-sequence=%s \\\n", args->seqfile);
    fprintf(fp, "\t-measurement=%s \\\n", args->measfile);

    write_option_outdir_runlog(fp, args->outdir, args->runlog);

    if (is_empty_string(args->dictfile))
        write_option_topo(fp, args->topo, args->topo_list);
    else
        fprintf(fp, "\t-dictfile=%s \\\n", args->dictfile);

    fprintf(fp, "\t-p_value=%g \\\n", args->p_value);
    fprintf(fp, "\t-max_motif=%ld \\\n", args->max_motif);
    fprintf(fp, "\t-strand=%ld \\\n", args->strand);
    fprintf(fp, "\t-ntop=%ld\n", args->ntop);

    close_file(fp);
}

static void parse_iupac_symbol(char *symbols)
{
    char *iupac2 = "KMRSWY";  /* 6 */
    char *iupac3 = "KMRSWYBDHV";  /* 6 + 4 = 10 */
    char *iupac4 = "KMRSWYBDHVN";  /* 6 + 4 + 1 = 11 */

    if (str_pmatch(symbols, "ALL"))
        strcpy(symbols, iupac4);
    else if (*symbols == '2' || str_pmatch(symbols, "IUPAC2"))
        strcpy(symbols, iupac2);
    else if (*symbols == '3' || str_pmatch(symbols, "IUPAC3")
             || str_pmatch(symbols, "NON"))
        strcpy(symbols, iupac3);

    if (!string_contains_only_those_characters(symbols, iupac4))
        fatal("Invalid IUPAC symbol(s) in [%s] vs [%s]\n", symbols, iupac4);
}

void check_reduce_cmdline(args_reduce * args)
{
    check_required_file(args->seqfile, "", "-sequence=seqfile");
    check_required_file(args->measfile, "", "-measurement=measfile");

    if (is_empty_string(args->dictfile))
        check_topo_option(args->topo_list, args->topo);
    else
        check_required_file(args->dictfile, "", "-dictfile=motif_list");

    log_strand_msg(args->strand);

    if (args->iupac_pos) {
        char str[BUF512];

        parse_iupac_symbol(args->iupac_sym);
        sprintf(str, "Check for IUPAC degeneracy: position=%ld\tsymbols=%s",
                args->iupac_pos, args->iupac_sym);
        log_prg(str);
    }
}

void reduce_cmdline(int argc, char *argv[], args_reduce * args)
{
    char helpfile[BUF512];
    long i;

    set_reduce_defaults(args);
    bname_ext(Gvars.PROGNAME, "hlp", helpfile);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != DASH)
            break;

        if (check_common_options(argv[i], helpfile, args->runlog, args->outdir))
            continue;

        if (extract_option_topo(argv[i], args->topo, args->topo_list))
            continue;

        if (str_pmatch(argv[i], "-s=") || str_pmatch(argv[i], "-se"))  /* sequence file name */
            get_strvalue(argv[i], args->seqfile, TRUE);

        /* measurements file; -ex for backwards compatibility */
        else if (str_pmatch(argv[i], "-m=") || str_pmatch(argv[i], "-me") ||
                 str_pmatch(argv[i], "-ex"))
            get_strvalue(argv[i], args->measfile, TRUE);

        else if (str_pmatch(argv[i], "-d"))
            get_strvalue(argv[i], args->dictfile, TRUE);

        else if (str_pmatch(argv[i], "-p=") || str_pmatch(argv[i], "-pv") ||
                 str_pmatch(argv[i], "-p_v"))
            args->p_value = get_dvalue(argv[i], 0.0, 1.0);

        else if (str_pmatch(argv[i], "-ma"))
            args->max_motif = get_lvalue(argv[i], 0, BUF512);

        else if (str_pmatch(argv[i], "-ip") || str_pmatch(argv[i], "-iupac_p"))
            args->iupac_pos = get_lvalue(argv[i], 0, BUF512);

        else if (str_pmatch(argv[i], "-is") || str_pmatch(argv[i], "-iupac_s")) {
            if (has_no_equal_sign(argv[i]))
                strcpy(args->iupac_sym, "ALL");
            else {
                get_strvalue(argv[i], args->iupac_sym, FALSE);
                if (is_empty_string(args->iupac_sym))  /* for case: -iupac_sym= */
                    strcpy(args->iupac_sym, "ALL");
                else
                    upperstr(args->iupac_sym);
            }

        } else if (str_pmatch(argv[i], "-nt"))
            args->ntop = get_lvalue(argv[i], 0, BUF512);

        else if (str_pmatch(argv[i], "-ns"))
            args->nsel = set_switch_with_dft_true(argv[i]);

        else if (str_pmatch(argv[i], "-st"))
            args->strand = extract_strand_option(argv[i]);

        else
            display_unrecognized_option(argv[i]);
    }

    set_logs_and_check_option(args->runlog, args->outdir, argc, i, argv[i]);

    check_reduce_cmdline(args);
    write_reduce_options(args);
}

void output_reduce_results(args_reduce * args, long psam_num, struct_seed * smotifs,
                           struct_data * tdat, struct_fitpars * lfpars, char *type)
{
    char filename[BUF512], raw[BUF512], xml[BUF512], *outdir = args->outdir;

    if (case_strstr(Gvars.PROGNAME, "matrix") != NULL) {
        strcpy(raw, MATRIXREDUCE_RAW);
        strcpy(xml, MATRIXREDUCE_XML);

    } else {
        strcpy(raw, MOTIFREDUCE_RAW);
        strcpy(xml, MOTIFREDUCE_XML);
    }

    write_FtP(type, outdir, basename(args->measfile), tdat, psam_num, lfpars);

    sprintf(filename, "%s/%s", outdir, raw);
    write2raw(filename, args, tdat, psam_num, smotifs, lfpars);

    sprintf(filename, "%s/%s", outdir, xml);
    write2xml(filename, args, tdat, psam_num, smotifs, lfpars);

    write_list_resid_predict(psam_num, smotifs, outdir, tdat);
}
