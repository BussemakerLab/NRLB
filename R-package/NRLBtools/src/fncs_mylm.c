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

void populate_A_g_chisq(long ma, double *a, long *ia, struct_lm * lmdata, double *g,
                        double **A, double *chisq)
{
    long i, nWs = ma - 2;
    double dy, ymod, dval;
    double a2[BUF512] = { 0.0 };  /* 0-indexed */
    double da2[BUF512] = { 0.0 };  /* 0-indexed */
    double dyda[BUF512] = { 0.0 };  /* 1-indexed, as a[] and ia[] */
    double N, C, F;

    *chisq = 0.0;

    for (i = 1; i <= nWs; i++) {
        dy = a[i];
        a2[i - 1] = dy * dy;  /* a2[] and da2[] are 0-indexed */
        da2[i - 1] = 2.0 / dy;
    }

    F = a[ma - 1];
    C = a[ma];

    for (i = 1; i <= lmdata->nrow; i++) {
        dval = lmdata->yval[i];

        if (dval < XBIG_CUTOFF) {
            N = my_fdf(lmdata, i, ma, a2, da2, F, dyda);
            ymod = C + F * N;
            dyda[ma - 1] = N;  /* dF -- slope */
            dyda[ma] = 1.0;  /* dC -- intercept */

            dy = dval - ymod;

            update_A_g_chisq_per_data_point(ma, ia, dyda, dy, g, A, chisq);
        }
    }
}

double my_fdf(struct_lm * lmdata, long idx, long nWs, double *a2, double *da2,
              double F, double *dyda)
{
    char *p;
    long i, j, m, nb, nseq, eidx, didx[BUF512], strand = lmdata->strand;
    double W, df[BUF512] = { 0.0 };
    double N = 0.0;

    struct_seq *seq = &lmdata->seqs->sequences[idx];
    struct_topo *t = lmdata->topo;

    eidx = t->Tcount - 1;

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {
        nb = seq->num_nb[nseq];
        p = seq->seq_fragments[nseq];

        if (strand == 1 || strand == 2)  /* leading strand or both strands */
            for (m = 0; m <= nb - t->Tcount; m++) {
                if (t->num_parens && !with_matched_parens_bin(t, p + m))
                    continue;

                W = 1.0;
                for (i = 0; i < t->Xcount; i++) {
                    j = i * NUM_BASE4 + p[m + t->offset[i]];
                    W *= a2[j];
                    didx[i] = j;
                }

                N += W;
                for (i = 0; i < t->Xcount; i++)
                    df[didx[i]] += W * da2[didx[i]];
            }

        if (strand == -1 || strand == 2)  /* reverse cmpl strand or both */
            for (m = 0; m <= nb - t->Tcount; m++) {
                if (t->num_parens && !with_matched_parens_bin(t, p + m))
                    continue;  /* to be checked -- maybe not relevant */

                W = 1.0;
                for (i = 0; i < t->Xcount; i++) {
                    j = i * NUM_BASE4 + TOP_BIDX3 - p[m + eidx - t->offset[i]];
                    W *= a2[j];
                    didx[i] = j;
                }

                N += W;
                for (i = 0; i < t->Xcount; i++)
                    df[didx[i]] += W * da2[didx[i]];
            }
    }

    for (i = 0; i < nWs; i++)
        dyda[i + 1] = F * df[i];

    return N;
}

void update_A_g_chisq_per_data_point(long ma, long *ia, double *dyda, double dy,
                                     double *g, double **A, double *chisq)
{
    long i, j, k, m;

    i = 0;
    for (k = 1; k <= ma; k++) {
        if (!ia[k])
            continue;
        i++;
        g[i] += dy * dyda[k];

        j = 0;
        for (m = 1; m <= k; m++) {  /* only lower-half */
            if (!ia[m])
                continue;
            j++;
            A[i][j] += dyda[k] * dyda[m];
        }
    }

    *chisq += dy * dy;
}

void my_mrqcof(long ma, long mfit, double *a, long *ia, struct_lm * lmdata, double *g,
               double **A, double *chisq, double *ng)
{
    long i, j;

    for (i = 1; i <= mfit; i++) {
        g[i] = 0.0;
        for (j = 1; j <= mfit; j++)
            A[i][j] = 0.0;
    }

    populate_A_g_chisq(ma, a, ia, lmdata, g, A, chisq);

    for (i = 2; i <= mfit; i++)  /* get the other half of symmetry-related A */
        for (j = 1; j < i; j++)
            A[j][i] = A[i][j];

    *ng = norm_inf_dvector(g, 1, mfit);
}

static double get_mu(long mfit, double tau, double **A)
{
    long i;
    double dval = -XBIG;

    for (i = 1; i <= mfit; i++)
        if (A[i][i] > dval)
            dval = A[i][i];

    return tau * dval;
}

static double get_dL(long mfit, double mu, double *h, double *g)
{
    long i;
    double dL = 0.0;

    for (i = 1; i <= mfit; i++)
        dL += h[i] * (mu * h[i] + g[i]);

    return dL;
}

static void update_a(long ma, long *ia, double *a, double *h, double *a_new)
{
    long i, j = 0;

    for (i = 1; i <= ma; i++)
        if (ia[i]) {
            j++;
            a_new[i] = a[i] + h[j];
        } else
            a_new[i] = a[i];
}

static double norm_ia_dvector(long ma, double *a, long *ia)
{
    long i;
    double dsum = 0.0;

    for (i = 1; i <= ma; i++)
        if (ia[i])
            dsum += a[i] * a[i];

    return sqrt(dsum);
}

static void update_pars(struct_lm * lmdata, long k, long ma, long mfit, long *ia,
                        double *a, double *g, double *h, double **A, double *chisq,
                        double *ng, double *mu, double *nu, long num_print)
{
    char msg[BUF512];
    long i, j;
    double chisq_new, dchisq, ng_new, dL;
    double *a_new, *g_new, **A_new;

    a_new = dvector(1, ma);
    g_new = dvector(1, ma);
    A_new = dmatrix(1, ma, 1, ma);

    update_a(ma, ia, a, h, a_new);

    dL = get_dL(mfit, *mu, h, g);
    my_mrqcof(ma, mfit, a_new, ia, lmdata, g_new, A_new, &chisq_new, &ng_new);

    dchisq = *chisq - chisq_new;

    if (dL > 0 && dchisq > 0) {  /* update a[] and modify mu */
        sprintf(msg, "\t%4ld\t%10.4f\t%10.4f\tng=%g", k, *chisq, dchisq, ng_new);

        if (k % num_print == 0)
            fprintf(Gvars.RUNLOG, "%s\n", msg);
        fprintf(Gvars.PRGLOG, "%s\tmu=%g\n", msg, *mu);

        for (i = 1; i <= ma; i++)  /* ma, not mfit */
            a[i] = a_new[i];

        for (i = 1; i <= mfit; i++) {  /* mfit, not ma */
            g[i] = g_new[i];
            for (j = 1; j <= mfit; j++)
                A[i][j] = A_new[i][j];
        }

        *chisq = chisq_new;
        *ng = ng_new;

        *mu *= dval_max(1.0 / 3.0, 1 - pow(2.0 * dchisq / dL - 1.0, 3.0));
        *nu = 2.0;

    } else {
        *mu *= *nu;
        *nu *= 2.0;
    }

    free_dvector(a_new, 1, DUMMY);
    free_dvector(g_new, 1, DUMMY);
    free_dmatrix(A_new, 1, DUMMY, 1, DUMMY);
}

static void set_levmar_pars(struct_levmar * levmar_pars)
{
    char msg[BUF512], parfile[BUF512];
    FILE *fp;

    get_sysfile(parfile, "config", PKG_CFG);

    sprintf(msg, " ...... reading file: %s ......", parfile);
    log_msg_vchk(msg);

    fp = open_file(parfile, "r");

    levmar_pars->delta_g = extract_xml_line_double(fp, "delta_g", 0.05, FALSE);
    levmar_pars->delta_x = extract_xml_line_double(fp, "delta_x", 1.0e-10, FALSE);
    levmar_pars->max_iter = extract_xml_line_long(fp, "max_iter", 6000, FALSE);
    levmar_pars->min_oknum = extract_xml_line_long(fp, "min_oknum", 6, FALSE);
    levmar_pars->num_print = extract_xml_line_long(fp, "num_print", 50, FALSE);

    close_file(fp);

    fprintf(Gvars.PRGLOG, "Levenberg-Marquardt nonlinear ls-fitting parameters\n");
    fprintf(Gvars.PRGLOG, "\tdelta_g: %g\n", levmar_pars->delta_g);
    fprintf(Gvars.PRGLOG, "\tdelta_x: %g\n", levmar_pars->delta_x);
    fprintf(Gvars.PRGLOG, "\tmax_iter: %ld\n", levmar_pars->max_iter);
    fprintf(Gvars.PRGLOG, "\tmin_oknum: %ld\n", levmar_pars->min_oknum);
    fprintf(Gvars.PRGLOG, "\tnum_print: %ld\n", levmar_pars->num_print);
}

static void check_oknum(double dval, double dcut, long *oknum)
{
    if (dval > dcut + NR_EPS)  /* to account for round-off errors */
        *oknum = 0;  /* condition must be satisfied consecutively */
    else
        ++ * oknum;
}

static double get_position_Wmax_idx(double *pos, long *idx)
{
    long i, m = -1;
    double Wmax = -XBIG;

    for (i = 0; i < NUM_BASE4; i++)
        if (Wmax < pos[i]) {  /* not <=, thus A > C > G > T when equal */
            Wmax = pos[i];
            m = i;
        }

    *idx = m;

    return Wmax;
}

double get_position_Wmax(double *pos)
{
    long i;
    double Wmax = -XBIG;

    for (i = 0; i < NUM_BASE4; i++)
        if (Wmax < pos[i])
            Wmax = pos[i];

    return Wmax;
}

static long set_ia(long ma, double *a, long *ia)
{
    long i, j, k, m, fixed, mfit = ma, idx = 1, nW = ma - 2;

    init_lvector(ia, 1, ma, 1);  /* default to optimize all parameters */

    if (Gvars.misc.FIT4) {
        if (Gvars.misc.POS[BUF510]) {
            for (i = 0; i < nW / NUM_BASE4; i++) {
                fixed = Gvars.misc.POS[i] == -1;
                for (j = 0; j < NUM_BASE4; j++) {
                    k = idx + j;
                    if (fixed) {
                        ia[k] = 0;  /* ia[] is 1-indexed */
                        mfit--;
                    }
                }
                idx += NUM_BASE4;  /* move to the next position */
            }
        }

    } else {
        double Wmax, *a2;

        a2 = dvector(1, nW);  /* only Ws, excluding F & C */
        for (i = 1; i <= nW; i++)  /* a[] is 1-indexed */
            a2[i] = a[i] * a[i];

        for (i = 0; i < nW / NUM_BASE4; i++) {
            Wmax = get_position_Wmax_idx(&a2[idx], &m);
            UNUSED_PARAMETER(Wmax);

            for (j = 0; j < NUM_BASE4; j++) {
                k = idx + j;
                if (j == m) {
                    ia[k] = 0;  /* keep the maximum for current position fixed */
                    mfit--;
                }
            }
            idx += NUM_BASE4;  /* move to the next position */
        }

        free_dvector(a2, 1, DUMMY);
    }

    return mfit;
}

/* Rescale W for a position where its maximum is over a preset limit
 * to avoid wrongly chosen initial Wmax. F is scaled accordingly. */
static long rescale_posW_and_F(long ma, double *a)
{
    char msg[BUF512];
    long i, j, k = FALSE, idx = 1, nW = ma - 2;  /* number of Ws */
    double Wmax, *a2;

    a2 = dvector(1, nW);  /* a[] is 1-indexed */
    for (i = 1; i <= nW; i++) {
        a2[i] = a[i] * a[i];
        if (a2[i] > WMAX_CUTOFF)
            k = TRUE;  /* need to rescale */
    }

    if (!k) {  /* no need for rescaling */
        free_dvector(a2, 1, DUMMY);
        return FALSE;
    }

    /* rescale PSAM by making the maximum W = 1.0 for each position */
    for (i = 0; i < nW / NUM_BASE4; i++) {
        Wmax = get_position_Wmax(&a2[idx]);

        if (Wmax != 1.0) {
            sprintf(msg, "\t\t>>rescale position %ld, maximum W = %g<<", i + 1, Wmax);
            log_msg(msg);

            for (j = 0; j < NUM_BASE4; j++) {
                k = idx + j;
                a[k] = sqrt(a2[k] / Wmax);  /* rescaled Ws */
            }

            a[ma - 1] *= Wmax;  /* update F accordingly */
        }

        idx += NUM_BASE4;  /* move to the next position */
    }

    free_dvector(a2, 1, DUMMY);

    return TRUE;
}

static void lm_covar(long ma, long mfit, long num_data, double *a, long *ia, double **A,
                     double chisq)
{
    long i;
    double se, sigma, tval, **covar;

    if (Gvars.VERBOSITY < VDEBUG)
        return;

    sigma = sqrt(chisq / (num_data - mfit));

    covar = dmatrix(1, ma, 1, ma);
    get_covar(mfit, A, covar);  /* the inverse of A */

    if (ia != NULL)
        NRC_covsrt(covar, ma, ia, mfit);

    for (i = 1; i <= ma; i++) {
        if (covar[i][i] <= TOLX)
            se = tval = 0.0;

        else {
            se = sqrt(covar[i][i]) * sigma;
            tval = a[i] / se;
        }

        fprintf(Gvars.PRGLOG, "%3ld\t%12.6f\t%12.6f\t%12.6f\n", i, a[i], se, tval);
    }

    free_dmatrix(covar, 1, DUMMY, 1, DUMMY);
}

static void initialize_lm(long ma, long mfit, double *a, long *ia, struct_lm * lmdata,
                          double *g, double **A, double *chisq, double *ng, double *mu)
{
    double tau = 1.0e-3;

    my_mrqcof(ma, mfit, a, ia, lmdata, g, A, chisq, ng);

    *mu = get_mu(mfit, tau, A);
}

void do_mylm(long ma, double *a, struct_lm * lmdata)
{
    char msg[BUF512];
    long k, mfit, num_print, oknum_dg = 0, oknum_dx = 0, *ia;
    double nu = 2.0, mu_min = 1.0e-14;
    double dg, dx, dx_chk, mu, chisq, na, ng, nh, *h, *g, **A;
    struct_levmar levmar_pars;

    set_levmar_pars(&levmar_pars);
    dg = levmar_pars.delta_g;
    dx = levmar_pars.delta_x;
    num_print = levmar_pars.num_print;

    h = dvector(1, ma);
    g = dvector(1, ma);
    A = dmatrix(1, ma, 1, ma);

    ia = lvector(1, ma);
    mfit = set_ia(ma, a, ia);
    initialize_lm(ma, mfit, a, ia, lmdata, g, A, &chisq, &ng, &mu);

    for (k = 1; k <= levmar_pars.max_iter; k++) {
        if (rescale_posW_and_F(ma, a)) {
            mfit = set_ia(ma, a, ia);
            initialize_lm(ma, mfit, a, ia, lmdata, g, A, &chisq, &ng, &mu);
        }

        check_oknum(ng, dg, &oknum_dg);
        if (ng <= dg && oknum_dg >= levmar_pars.min_oknum) {
            sprintf(msg, "\t%4ld (%g): converged with gradient: %g <= %g",
                    k, chisq, ng, dg);
            log_msg(msg);
            break;
        }

        mu = dval_max(mu, get_mu(mfit, mu_min, A));
        lux_Axsolution(mfit, A, g, mu, h);

        nh = norm_dvector(h, 1, mfit);
        na = norm_ia_dvector(ma, a, ia);
        dx_chk = dx * (dx + na);

        check_oknum(nh, dx_chk, &oknum_dx);
        if (nh <= dx_chk && oknum_dx >= levmar_pars.min_oknum) {
            sprintf(msg, "\t%ld: stopped with dx: %g\t%g", k, nh, dx_chk);
            log_msg(msg);
            break;
        }

        update_pars(lmdata, k, ma, mfit, ia, a, g, h, A, &chisq, &ng, &mu, &nu,
                    num_print);
    }

    lm_covar(ma, mfit, lmdata->num, a, ia, A, chisq);

    free_lvector(ia, 1, DUMMY);
    free_dvector(h, 1, DUMMY);
    free_dvector(g, 1, DUMMY);
    free_dmatrix(A, 1, DUMMY, 1, DUMMY);
}

static void init_lmseed_x0(long num_pars, double W0, double intercept, double slope,
                           double *x0, double *x)
/* x[] is 1-indexed, with minimum set to W0; x0 is 0-based */
{
    double dval;
    long i, j, k, fixed, idx = 0, nW = num_pars - 2;

    x[num_pars - 1] = slope;
    x[num_pars] = intercept;  /* last position */

    for (i = 0; i < nW / NUM_BASE4; i++) {
        fixed = Gvars.misc.POS[BUF510] && Gvars.misc.POS[i] == -1;
        for (j = 0; j < NUM_BASE4; j++) {
            k = idx + j;
            dval = x0[k];
            x[k + 1] = sqrt(dval);  /* x[] is 1-indexed */
            if (!fixed && dval == 0.0)  /* optimized position */
                x[k + 1] = W0;
        }
        idx += NUM_BASE4;
    }
}

void init_lmseed(long num_pars, double W0, struct_seed * seed)
{
    double *x0;

    x0 = dvector(0, num_pars - 3);  /* w is 0-indexed, excluding F & C */

    if (Gvars.misc.ALLPOS)
        iupac2psam(seed->full, x0);
    else
        iupac2psam(seed->motif, x0);

    init_lmseed_x0(num_pars, W0, seed->fitpars.a, seed->fitpars.b, x0, seed->x);

    free_dvector(x0, 0, DUMMY);
}

void check_opt_allpos(struct_topo * t)
/* Tentative option to optimize all seed positions regardless of
 * topological pattern. */
{
    long i, k = t->Tcount;

    if (!Gvars.misc.ALLPOS)
        return;  /* default follows the topo-pattern */

    log_msg("\tall positions\n");

    t->Xcount = k;
    t->num_parens = 0;

    for (i = 0; i < k; i++)
        t->offset[i] = i;

    t->msize = 1 << 2 * k;
}

void lm_optimize_seed(struct_seed * seed, struct_data * tdat, struct_seqArr * seqs)
{
    long num_pars;
    double W0 = sqrt(WEPS);
    struct_lm lmdata;  /* passing parameters to LM optimizer */

    if (seed->topo->Xcount == 0)
        return;

    log_msg("Optimizing:");

    check_opt_allpos(seed->topo);

    num_pars = NUM_BASE4 * seed->topo->Xcount + 2;  /* plus F & C */
    seed->x = dvector(1, num_pars);  /* parameters to be optimized */

    init_lmseed(num_pars, W0, seed);

    lmdata.nrow = tdat->nrow;
    lmdata.num = seed->fitpars.num;
    lmdata.strand = seed->strand;
    lmdata.topo = seed->topo;
    lmdata.yval = tdat->resid[seed->col_idx];
    lmdata.seqs = seqs;

    debug_parenthesis(seed->topo->str_tag.str, seed->topo->num_parens,
                      seed->topo->parens_bidx, seed->topo->parens_eidx);

    do_mylm(num_pars, seed->x, &lmdata);

    seed->w = dvector(0, num_pars - 3);  /* 0-indexed, excluding F & C */
    normalize_psam(num_pars, seed->x, seed->w);  /* only optimized positions */
}

void lm_optimize_expt_psam(struct_seed * seed, struct_data * tdat, struct_seqArr * seqs)
{
    long num_pars;
    double W0 = sqrt(WEPS);
    struct_lm lmdata;

    log_msg("Optimizing:");

    check_opt_allpos(seed->topo);

    num_pars = NUM_BASE4 * seed->topo->Xcount + 2;
    seed->x = dvector(1, num_pars);

    init_lmseed_x0(num_pars, W0, seed->fitpars.a, seed->fitpars.b, seed->w, seed->x);

    lmdata.nrow = tdat->nrow;
    lmdata.num = seed->fitpars.num;
    lmdata.strand = seed->strand;
    lmdata.topo = seed->topo;
    lmdata.yval = tdat->resid[seed->col_idx];
    lmdata.seqs = seqs;

    do_mylm(num_pars, seed->x, &lmdata);

    normalize_psam(num_pars, seed->x, seed->w);
}
