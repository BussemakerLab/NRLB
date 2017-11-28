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

long is_motif_to_skip(long idx, double *sumX, char *cklist, long dyn_stnd,
                      long strand, long *rc_lookup)
{
    if (sumX[idx] < Gvars.misc.MIN_COUNTS)  /* too few matched motifs */
        return TRUE;

    if (cklist && !cklist[idx])  /* not in the selected list */
        return TRUE;

    if (Gvars.PALINDROME && idx != rc_lookup[idx])
        return TRUE;

    /* When auto-determining strand-ness, +1 and -1 have 1-to-1
     * correspondence, thus the -1 direction needs not to be checked */
    if (dyn_stnd && strand != 1 && idx >= rc_lookup[idx])
        return TRUE;  /* to avoid palindrome ATATAT */

    return FALSE;
}

long is_nonbit_motif_to_skip(long sumX, char *motif, long dyn_stnd, long strand)
{
    if (sumX < Gvars.misc.MIN_COUNTS)
        return TRUE;

    if (Gvars.PALINDROME && !is_rc_palindrome(motif))
        return TRUE;

    if (dyn_stnd && strand != 1 && is_rc_palindrome(motif))
        return TRUE;

    return FALSE;
}

/* check the ls-fitted parameters of each motif, and keep only the best one */
static void get_mfit_per_topo_expt(long topo_idx, long nc, long num, long msize,
                                   long *rc_lookup, double *sumX, double *sumX2,
                                   double *sumXY, double sY, double sY2, long strand,
                                   char *cklist, long dyn_stnd,
                                   struct_ntopArr * ntop_seeds)
{
    long i, tnum, lnum, bonf = 0, ntop = ntop_seeds->ntop;
    double b_tval;  /* t-value of the slope */

    struct_ntop *nte;  /* to simplify code: n-top-entry */

    tnum = ntop + 1;  /* total number of entries: [0, ntop] */
    lnum = ntop - 1;  /* index of last used entry */

    for (i = 0; i < msize; i++) {
        if (is_motif_to_skip(i, sumX, cklist, dyn_stnd, strand, rc_lookup))
            continue;

        bonf++;  /* number of tests performed */
        b_tval = get_slope_abs_tvalue(num, sumX[i], sumX2[i], sumXY[i], sY, sY2);

        /* Add NR_EPS to avoid round-off led ambiguity for palindrome
         * sequence, e.g., ACGCGT with -topo=X6 in the Spellman dataset
         * with -strand=0 setting */
        if (b_tval > ntop_seeds->ntops[lnum].b_tval + NR_EPS) {
            nte = &ntop_seeds->ntops[ntop];
            nte->col_idx = nc;
            nte->topo_idx = topo_idx;
            nte->index = i;  /* decimal ==> to deduce motif sequence */
            nte->b_tval = b_tval;  /* absolute value */
            nte->strand = strand;  /* the strand-ness for dynamic option */
            qsort(ntop_seeds->ntops, tnum, sizeof(struct_ntop), ntop_compare);
        }
    }

    ntop_seeds->bonf += bonf;  /* only ONE bonf counts */
}

void populate_mcount(struct_topo * t, struct_seq * seq, long *mcount)
{
    char *p;
    long i, idx, m, nb, nseq;

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {
        nb = seq->num_nb[nseq];  /* # of bases in this fragment */
        p = seq->seq_fragments[nseq];  /* pointer to the sequence [0123] */

        for (m = 0; m <= nb - t->Tcount; m++) {  /* # of windows along the sequence */
            if (t->num_parens && !with_matched_parens_bin(t, p + m))
                continue;

            idx = 0;
            for (i = 0; i < t->Xcount; i++) {  /* 0-based index */
                idx <<= 2;
                idx += p[m + t->offset[i]];
            }
            mcount[idx]++;
        }
    }
}

/* get all possible motif counts of pattern 't' in a sequence */
void get_mcount_per_seq(struct_topo * t, double dval, struct_seq * seq, long msize,
                        long *mcount, long *rc_lookup, long strand, double *sumX,
                        double *sumX2, double *sumXY)
{
    long i, idx, m = 0;

    for (i = 0; i < msize; i++)
        mcount[i] = 0;

    populate_mcount(t, seq, mcount);

    if (strand == 1) {  /* leading strand */
        for (i = 0; i < msize; i++) {  /* NB: '<', not '<=' */
            m = mcount[i];  /* count of motif represented by decimal 'i' */
            if (m) {  /* the sequence has motif 'i' */
                sumX[i] += m;
                sumX2[i] += m * m;
                sumXY[i] += m * dval;
            }
        }

    } else if (strand == -1) {  /* reverse complementary strand */
        for (i = 0; i < msize; i++) {
            m = mcount[rc_lookup[i]];
            if (m) {
                sumX[i] += m;
                sumX2[i] += m * m;
                sumXY[i] += m * dval;
            }
        }

    } else if (strand == 2) {  /* both strand -- check to eliminate degeneracy */
        for (i = 0; i < msize; i++) {
            idx = rc_lookup[i];
            if (i > idx)  /* idx == i for palindrome sequence, e.g. ATATAT */
                continue;

            m = mcount[i] + mcount[idx];
            if (m) {
                sumX[i] += m;
                sumX2[i] += m * m;
                sumXY[i] += m * dval;
            }
        }

    } else
        fatal("wrong strand value <%ld>: [-1, 1, 2]", strand);
}

void init_fitpars(struct_ls5fit * fitpars)
{
    fitpars->num = 0;
    fitpars->sumX = 0.0;
    fitpars->SSY = 0.0;
    fitpars->SSE = 0.0;
    fitpars->SSR = 0.0;
    fitpars->r2 = 0.0;

    fitpars->F = 0.0;
    fitpars->a = 0.0;
    fitpars->a_tval = 0.0;
    fitpars->a_pval = 0.0;
    fitpars->b = 0.0;
    fitpars->b_tval = 0.0;
    fitpars->b_pval = 0.0;
}

void init_seed_motif(struct_seed * seed, long strand, double p_value)
{
    init_fitpars(&seed->fitpars);

    seed->x = NULL;
    seed->w = NULL;
    seed->psam = NULL;

    seed->num_mids = 0;
    seed->bonf = 0;  /* for multiple test correction */
    seed->col_idx = 0;  /* experiment column in the measurement profile */
    seed->col_name = NULL;

    seed->topo_idx = DUMMY;  /* index of the topological pattern */
    seed->topo = NULL;

    seed->index = -1;  /* motif decimal index: to decode for sequence */
    seed->motif = NULL;
    seed->full = NULL;
    seed->optimal = NULL;

    seed->strand = strand;
    seed->stnd_msg = NULL;
    seed->p_value = p_value;
}

long *populate_rc_lookup(long msize, long Xcount)
{
    long i, *rc_lookup;

    rc_lookup = lvector(0, msize);

    for (i = 0; i < msize; i++)
        rc_lookup[i] = fast_cidx(i, Xcount);

    return rc_lookup;
}

static void get_nogap_topo(long k, struct_topo * t)
{
    char pattern[BUF512];
    long i;

    repeat_char_ntimes_string(XCHR, k, pattern);

    init_strtag_with_strings(pattern, NULL, &t->str_tag);

    t->Tcount = k;
    t->Xcount = k;
    t->num_parens = 0;  /* w/o () pairs */

    t->offset = lvector(0, k);
    for (i = 0; i < k; i++)
        t->offset[i] = i;

    if (k > MAXNTS_DIC)
        t->msize = -1;
    else
        t->msize = 1 << 2 * k;  /* 4^Xcount; ONE element more */
}

static void get_motif_fitpars(struct_data * tdat, struct_seqArr * seqs, char *motif,
                              struct_seed * seed, struct_ls5fit * fitpars, long debug)
{
    double *y, *counts;

    counts = dvector(1, seqs->num_seqs);
    y = tdat->resid[seed->col_idx];

    get_motif_counts(motif, seed->strand, seqs, counts);

    get_y_vs_x_lsfitpars(tdat->nrow, counts, y, fitpars);

    seed->num_mids = get_num_mids(tdat->nrow, counts, y);

    free_dvector(counts, 1, DUMMY);

    if (debug) {
        char str[BUF512];

        sprintf(str, "Motif <%s>", motif);
        print_fitpars(str, fitpars);
    }
}

static void populate_seed_from_motif(struct_data * tdat, struct_seqArr * seqs,
                                     struct_seed * seed, struct_dic * mlist)
{
    char msg[BUF512], *motif;
    long k;

    k = seed->index;  /* link to struct_dic arrays */
    motif = mlist->motifs[k].str;

    sprintf(msg, "\t\t...%s\t%ld", motif, k);
    log_msg_vchk(msg);

    get_motif_fitpars(tdat, seqs, motif, seed, &seed->fitpars, FALSE);

    expanded_motif_to_seedtopo(motif, seed);

    seed->motif = my_strdup(motif);
    seed->full = my_strdup(motif);

    seed->col_name = my_strdup(tdat->col_names[seed->col_idx]);
    seed->stnd_msg = strand2msg(seed->strand);
}

void populate_seed_fitpars(struct_data * tdat, struct_seqArr * seqs,
                           struct_topo * topos, struct_seed * seed)
{
    long sidx, nc, nr, msize, strand;
    long *mcount, *rc_lookup, m = 0, num_mids = 0;
    double dval, sX, sX2, sXY, sY, sY2;

    struct_topo *t = &topos[seed->topo_idx];

    nc = seed->col_idx;

    msize = t->msize;
    mcount = lvector(0, msize);

    rc_lookup = populate_rc_lookup(msize, t->Xcount);

    sidx = seed->index;
    strand = seed->strand;

    sX = sX2 = sXY = sY = sY2 = 0.0;

    for (nr = 1; nr <= tdat->nrow; nr++) {  /* loop over all measurement values */
        dval = tdat->resid[nc][nr];

        if (dval < XBIG_CUTOFF) {
            sY += dval;
            sY2 += dval * dval;

            populate_mcount(t, &seqs->sequences[nr], mcount);

            if (strand == 1)
                m = mcount[sidx];
            else if (strand == -1)
                m = mcount[rc_lookup[sidx]];
            else if (strand == 2)
                m = mcount[sidx] + mcount[rc_lookup[sidx]];
            else
                fatal("wrong strand value <%ld> [-1, 1, 2]", strand);

            if (m) {
                sX += m;
                sX2 += m * m;
                sXY += m * dval;
                num_mids++;

                mcount[sidx] = mcount[rc_lookup[sidx]] = 0;
            }
        }
    }

    free_lvector(rc_lookup, 0, DUMMY);
    free_lvector(mcount, 0, DUMMY);

    seed->num_mids = num_mids;
    lsfit_parameters(tdat->oknum[nc], sX, sX2, sXY, sY, sY2, &seed->fitpars);
}

void check_motif_for_iupac_degeneracy(args_reduce * args, struct_seed * seed,
                                      struct_data * tdat, struct_seqArr * seqs)
{
    char str[BUF512], motif[BUF512], motif_max[BUF512];
    double tval, tval_max = fabs(seed->fitpars.b_tval);
    double *y = tdat->resid[seed->col_idx], *counts;
    long num = 0, num0 = strlen(seed->motif), numa = strlen(seed->full);
    long num_sym = strlen(args->iupac_sym);
    long i, j, k, pos_max = 0, idx[BUF512] = { FALSE };  /* position index */

    if (args->iupac_pos == 0)
        return;

    if (args->iupac_pos > num0) {
        sprintf(str, "\nInput -iupac_pos=%ld over non-gapped motif length [%ld];"
                " reset it to all positions.", args->iupac_pos, num0);
        log_msg(str);
        args->iupac_pos = num0;
    }

    sprintf(str, "Checking IUPAC degeneracy:");
    log_msg(str);

    for (i = 0; i < numa; i++)
        if (strchr("-()", seed->full[i]))
            idx[i] = TRUE;  /* fixed positions */

    counts = dvector(1, seqs->num_seqs);
    while (num < args->iupac_pos) {
        k = FALSE;
        for (i = 0; i < numa; i++) {
            if (idx[i])
                continue;
            for (j = 0; j < num_sym; j++) {
                strcpy(motif, seed->full);  /* one mutation per iteration */
                motif[i] = args->iupac_sym[j];

                get_motif_counts(motif, seed->strand, seqs, counts);

                tval = get_y_vs_x_slope_abs_tvalue(tdat->nrow, counts, y);
                if (tval > tval_max) {
                    tval_max = tval;
                    k = TRUE;
                    pos_max = i;
                    strcpy(motif_max, motif);
                }
            }
        }
        if (!k)
            break;
        num++;
        idx[pos_max] = TRUE;
        strcpy(seed->full, motif_max);
        sprintf(str, "\titeration=%ld\tmotif=%s\t|t-value|=%g", num, motif_max, tval_max);
        log_msg(str);
    }

    get_motif_counts(seed->full, seed->strand, seqs, counts);

    get_y_vs_x_lsfitpars(tdat->nrow, counts, y, &seed->fitpars);
    seed->num_mids = get_num_mids(tdat->nrow, counts, y);
    revert_full_to_motif(seed);

    free_dvector(counts, 1, DUMMY);
}

void deduce_motif_full(struct_topo * t, long idx, char *motif, char *full)
{
    long i, j;

    convert_dec2seq(idx, t->Xcount, motif);

    for (i = 0; i < t->Tcount; i++)
        full[i] = DASH;  /* as default */
    full[i] = '\0';

    for (i = 0; i < t->Xcount; i++) {
        j = t->offset[i];
        full[j] = motif[i];
    }

    for (i = 1; i <= t->num_parens; i++) {
        full[t->parens_bidx[i]] = '(';
        full[t->parens_eidx[i]] = ')';
    }
}

void revert_full_to_motif(struct_seed * seed)
{
    long i, j;
    struct_topo *t = seed->topo;

    for (i = 0; i < t->Xcount; i++) {
        j = t->offset[i];
        seed->motif[i] = seed->full[j];
    }
}

void fulfill_smotif(struct_topo * topos, char **col_names, struct_seed * seed)
{
    struct_topo *t;  /* shorthand to simplify code */

    if (seed->index == -1)
        fatal("\tno seed motif found  -- motif length too long/specific?\n");

    t = &topos[seed->topo_idx];

    seed->topo = allocate_memory_for_one_topo();
    duplicate_topo(t, seed->topo);

    seed->motif = cvector(0, t->Xcount);
    seed->full = cvector(0, t->Tcount);

    deduce_motif_full(t, seed->index, seed->motif, seed->full);

    seed->col_name = my_strdup(col_names[seed->col_idx]);
    seed->stnd_msg = strand2msg(seed->strand);
}

void output_linefit_pars(FILE * fp, struct_ls5fit * fitpars)
{
    fprintf(fp, "\tintercept: coef=%+g\tt-value=%+g\tp-value=%g\n", fitpars->a,
            fitpars->a_tval, fitpars->a_pval);

    fprintf(fp, "\tslope:     coef=%+g\tt-value=%+g\tp-value=%g\n", fitpars->b,
            fitpars->b_tval, fitpars->b_pval);

    fprintf(fp, "\tr2=%g\tSSY=%g\tSSE=%g\tSSR=%g\n", fitpars->r2,
            fitpars->SSY, fitpars->SSE, fitpars->SSR);
}

void print_smotif(FILE * fp, struct_seed * seed, char *stnd)
{
    fprintf(fp, "Best seed motif:\n");

    output_linefit_pars(fp, &seed->fitpars);

    fprintf(fp, "\tmatches[matched_ids/all_ids]: %g[%ld/%ld]\tmotif: %s\n",
            seed->fitpars.sumX, seed->num_mids, seed->fitpars.num, seed->full);

    fprintf(fp, "\tbased on experiment coloum %ld; %s\n", seed->col_idx, seed->col_name);

    display_strand_msg(fp, "\t    and of sequence on", stnd, seed->strand);
}

void print_seed_expt(FILE * fp, struct_seed * seed, char *stnd)
{
    fprintf(fp, "Best seed experiment:\n");
    fprintf(fp, "\tnumber of tested candidate experiments: %ld\n", seed->bonf);

    output_linefit_pars(fp, &seed->fitpars);

    fprintf(fp, "\tmatches[matched-ids/total-ids]: %g[%ld/%ld]\texperiment: %s [%ld]\n",
            seed->fitpars.sumX, seed->num_mids, seed->fitpars.num, seed->col_name,
            seed->col_idx);

    display_strand_msg(fp, "\t    and of sequence on", stnd, seed->strand);
}

static void copy_ntop_seed_entry(struct_ntop * ntop, struct_seed * seed)
{
    seed->col_idx = ntop->col_idx;
    seed->topo_idx = ntop->topo_idx;
    seed->index = ntop->index;
    seed->strand = ntop->strand;
}

void setup_seed_from_topo_motif(struct_data * tdat, struct_seqArr * seqs,
                                struct_seed * seed, struct_ntop * nte,
                                struct_topo * topos, struct_dic * mlist)
{
    copy_ntop_seed_entry(nte, seed);

    if (seed->topo_idx > 0) {  /* normal topo-entries */
        populate_seed_fitpars(tdat, seqs, topos, seed);
        fulfill_smotif(topos, tdat->col_names, seed);

    } else  /* based directly on motif */
        populate_seed_from_motif(tdat, seqs, seed, mlist);
}

void print_ntop_entry(FILE * fp, long idx, struct_seed * seed, long maxlen)
{
    long is_rcp;

    is_rcp = is_rc_palindrome(seed->full);

    fprintf(fp, "\t%2ld\t%-*s", idx, (int) maxlen, seed->full);
    if (Gvars.seed_criterion == SEED_BY_SLOPE_ITSELF)
        fprintf(fp, "\t%+.3f(%+.3f)", seed->fitpars.b_tval, seed->fitpars.b);
    else
        fprintf(fp, "\t%+8.3f", seed->fitpars.b_tval);
    fprintf(fp, "\t%8.5f", seed->fitpars.r2);
    fprintf(fp, "\t%g[%ld/%ld]", seed->fitpars.sumX, seed->num_mids, seed->fitpars.num);
    fprintf(fp, "\t%-10s", is_rcp ? "palindrome" : seed->stnd_msg);
    fprintf(fp, "\t%s[%ld]\n", seed->col_name, seed->col_idx);
}

static long display_ntop(struct_data * tdat, struct_seqArr * seqs, struct_topo * topos,
                         struct_ntopArr * ntop_seeds, char *stnd, struct_dic * mlist)
{
    char msg[BUF512];
    long i, k = 0, nok = 0, ntop = ntop_seeds->ntop, nsel = ntop_seeds->nsel;
    long selected = 1;  /* default to the top seed */

    struct_seed *smotifs, *seed;
    struct_ntop *nte;  /* to simplify code: n-top-entry */

    for (i = 0; i < ntop; i++) {
        if (ntop_seeds->ntops[i].index == -1)
            break;
        nok++;
    }

    if (!nok)
        fatal("\tno seed motif found -- motif length too long/specific?\n");

    sprintf(msg, "List of top %ld best matched seeds%s:\n    Note: column 'counts'"
            " stand for number of matches[matched_ids/all_ids]", nok, stnd);
    log_msg(msg);

#if LUXDBG == TRUE
    for (i = 0; i < nok; i++) {
        nte = &ntop_seeds->ntops[i];
        sprintf(msg, "\t%ld\t%+g\t%ld\t%ld\t%ld\t%+ld", i + 1, nte->b_tval,
                nte->topo_idx, nte->index, nte->col_idx, nte->strand);
        log_msg_vchk(msg);
    }
    fprintf(Gvars.PRGLOG, "----------- above is a list of raw data -----------\n\n");
#endif

    smotifs = allocate_memory_for_smotifs(nok);

    for (i = 0; i < nok; i++) {
        nte = &ntop_seeds->ntops[i];

        seed = &smotifs[i];  /* as a shorthand notation */
        init_seed_motif(seed, nte->strand, 1.0);  /* p_value does not matter here */

        setup_seed_from_topo_motif(tdat, seqs, seed, nte, topos, mlist);
        k = lval_max(k, strlen(seed->full));
    }

    sprintf(msg, "\trank\tmotif\t");
    if (k >= 8)
        strcat(msg, "\t");
    strcat(msg, " t-value\tr-squared\tcounts\t\tdirection\texperiment");
    log_msg(msg);

    for (i = 0; i < nok; i++) {
        seed = &smotifs[i];  /* as a shorthand notation */
        print_ntop_entry(Gvars.RUNLOG, i + 1, seed, k);
        print_ntop_entry(Gvars.PRGLOG, i + 1, seed, k);
    }

    if (nsel) {  /* manually picked up from the list */
        fprintf(stdout, "    Please pick up your seed [1..%ld]: ", nok);
        selected = readline_cvt2long(stdin);
        if (!lval_in_range(selected, 1, nok)) {
            sprintf(msg, "    Picked motif not in range [1..%ld]: use the top one", nok);
            log_msg(msg);
            selected = 1;
        }
    }

    free_smotifs(0, nok - 1, smotifs);  /* NOTE: 0 ---> nok - 1 */

    return selected;
}

static void handle_topo(long topo_idx, struct_topoArr * topo, struct_data * tdat,
                        struct_seqArr * seqs, struct_ntopArr * ntop_seeds, long strand,
                        char *cklist, long dyn_stnd)
{
    char msg[BUF512];
    long i, nc, nr, msize, debug, dbg2;
    long *mcount, *rc_lookup;
    double ratio = 0.0, dval, sY, sY2, *sumX, *sumX2, *sumXY;

    struct_topo *t = &topo->topos[topo_idx];

    msize = t->msize;  /* to simplify the code a little bit */

    mcount = lvector(0, msize);  /* 1 position more */
    sumX = dvector(0, msize);
    sumX2 = dvector(0, msize);
    sumXY = dvector(0, msize);

    rc_lookup = populate_rc_lookup(msize, t->Xcount);

    debug = t->Xcount >= DEBUG_NCOUNT;
    dbg2 = debug && ((Gvars.RUNLOG == stderr) || (Gvars.RUNLOG == stdout));

    for (nc = 1; nc <= tdat->ncol; nc++) {  /* loop over each experiment */
        sY = sY2 = 0.0;  /* measurement value */

        if (debug) {
            sprintf(msg, "\t\ttopo: %ld; experiment: %ld ", topo_idx, nc);
            fprintf(Gvars.RUNLOG, "%s", msg);
            fprintf(Gvars.PRGLOG, "%s", msg);
            ratio = 100.0 / (double) tdat->nrow;
        }

        for (nr = 1; nr <= tdat->nrow; nr++) {  /* loop over all measurement values */
            dval = tdat->resid[nc][nr];

            if (dbg2)
                fprintf(Gvars.RUNLOG, "\r\t\t\t\t\t\t%6.2f", nr * ratio);

            if (dval < XBIG_CUTOFF) {
                sY += dval;
                sY2 += dval * dval;
                get_mcount_per_seq(t, dval, &seqs->sequences[nr], msize, mcount,
                                   rc_lookup, strand, sumX, sumX2, sumXY);
            }
        }

        if (debug)
            log_msg("");

        get_mfit_per_topo_expt(topo_idx, nc, tdat->oknum[nc], msize, rc_lookup,
                               sumX, sumX2, sumXY, sY, sY2, strand, cklist,
                               dyn_stnd, ntop_seeds);

        for (i = 0; i < msize; i++)
            sumX[i] = sumX2[i] = sumXY[i] = 0.0;
    }

    free_lvector(mcount, 0, DUMMY);
    free_lvector(rc_lookup, 0, DUMMY);
    free_dvector(sumX, 0, DUMMY);
    free_dvector(sumX2, 0, DUMMY);
    free_dvector(sumXY, 0, DUMMY);
}

void check_ntops_populate_seed(struct_ntopArr * ntop_seeds, struct_topo * topos,
                               long strand0, struct_data * tdat, struct_seqArr * seqs,
                               struct_seed * seed, struct_dic * mlist)
{
    char stnd[BUF512];
    long selected;

    strcpy(stnd, (strand0 == 0) ? " [auto-determined]" : "");

    selected = display_ntop(tdat, seqs, topos, ntop_seeds, stnd, mlist);

    setup_seed_from_topo_motif(tdat, seqs, seed, &ntop_seeds->ntops[selected - 1], topos, mlist);  /* Note - 1 to account for 0-index  */
    seed->bonf = ntop_seeds->bonf;

    print_smotif(Gvars.RUNLOG, seed, stnd);
    print_smotif(Gvars.PRGLOG, seed, stnd);

    free_ntop_entries(ntop_seeds->ntops);
}

void initialize_struct_ntopArr(long ntop, long nsel, struct_ntopArr * ntop_seeds)
{
    ntop_seeds->ntop = ntop;
    ntop_seeds->nsel = nsel;
    ntop_seeds->bonf = 0;  /* to sum-up */
    ntop_seeds->ntops = allocate_memory_for_ntop_entries(ntop);
}

void find_seed_by_topos(struct_data * tdat, struct_seqArr * seqs, struct_topoArr * topo,
                        struct_seed * seed, long ntop, long nsel)
{
    long topo_idx, strand0 = seed->strand;  /* initial setting for strandness */

    struct_ntopArr ntop_seeds;

    initialize_struct_ntopArr(ntop, nsel, &ntop_seeds);

    for (topo_idx = 1; topo_idx <= topo->num_topo; topo_idx++) {
        if (strand0 == 0) {  /* dynamically check forward and both strands */
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, 1, NULL, TRUE);
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, 2, NULL, TRUE);

        } else  /* follow user's command-line setting */
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, strand0, NULL, FALSE);
    }

    check_ntops_populate_seed(&ntop_seeds, topo->topos, strand0, tdat, seqs, seed, NULL);
}

static void populate_check_list(struct_dic * mlist, long topo_idx, char *cklist)
{
    long i, ib, ie, idx;

    ib = mlist->seidx[topo_idx][1];
    ie = mlist->seidx[topo_idx][2];

    for (i = ib; i <= ie; i++) {  /* serial index in mlist */
        idx = mlist->index[i];
        if (idx != -1)  /* for ACGT-only motif, with length <= MAXNTS_DIC */
            cklist[idx] = '\1';
    }
}

static void handle_nonbit_motifs(struct_data * tdat, struct_seqArr * seqs,
                                 struct_ntopArr * ntop_seeds, struct_dic * mlist,
                                 long strand, long dyn_stnd)
{
    char msg[BUF512], *motif;
    long i, k, nc, lnum, tnum, bonf = 0, ntop = ntop_seeds->ntop;
    double b_tval, *counts;

    struct_ls5fit fitpars;
    struct_ntop *nte;  /* to simplify code: n-top-entry */

    tnum = ntop + 1;  /* total number of entries: [0, ntop] */
    lnum = ntop - 1;  /* index of last used entry */

    counts = dvector(1, seqs->num_seqs);

    for (i = 1; i <= mlist->num_nonbit; i++) {  /* loop over all motifs */
        k = mlist->idx_nonbit[i];
        motif = mlist->motifs[k].str;  /* motif string */

        if (Gvars.VERBOSITY == VDEBUG) {
            sprintf(msg, "\t%ld\t%s\t%ld", i, motif, k);
            log_msg(msg);
        }

        get_motif_counts(motif, strand, seqs, counts);

        for (nc = 1; nc <= tdat->ncol; nc++) {  /* loop over each experiment */
            get_y_vs_x_lsfitpars(tdat->nrow, counts, tdat->resid[nc], &fitpars);

            if (is_nonbit_motif_to_skip(fitpars.sumX, motif, dyn_stnd, strand))
                continue;

            bonf++;

            b_tval = fabs(fitpars.b_tval);  /* absolute value */
            if (b_tval > ntop_seeds->ntops[lnum].b_tval + NR_EPS) {
                nte = &ntop_seeds->ntops[ntop];
                nte->col_idx = nc;
                nte->topo_idx = DUMMY;  /* non-topo entry */
                nte->index = k;  /* link to mlist arrays */
                nte->b_tval = b_tval;  /* absolute value */
                nte->strand = strand;  /* the strand-ness for dynamic option */
                qsort(ntop_seeds->ntops, tnum, sizeof(struct_ntop), ntop_compare);
            }
        }

        init_dvector(counts, 1, seqs->num_seqs, 0.0);  /* zero-out for next motif */
    }

    ntop_seeds->bonf += bonf;

    free_dvector(counts, 1, DUMMY);
}

void find_seed_by_dictfile(struct_data * tdat, struct_seqArr * seqs,
                           struct_topoArr * topo, struct_dic * mlist, struct_seed * seed,
                           long ntop, long nsel)
{
    char msg[BUF512], *cklist;
    long topo_idx, k, strand0 = seed->strand;  /* initial setting for strandness */

    struct_ntopArr ntop_seeds;

    initialize_struct_ntopArr(ntop, nsel, &ntop_seeds);

    for (topo_idx = 1; topo_idx <= topo->num_topo; topo_idx++) {  /* ideal cases */
        k = topo->topos[topo_idx].msize;

        cklist = cvector(0, k);  /* initialized with '\0' */
        populate_check_list(mlist, topo_idx, cklist);

        if (strand0 == 0) {  /* dynamically check forward and both strands */
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, 1, cklist, TRUE);
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, 2, cklist, TRUE);

        } else  /* follow user's command-line setting */
            handle_topo(topo_idx, topo, tdat, seqs, &ntop_seeds, strand0, cklist, FALSE);

        free_cvector(cklist, 0, DUMMY);
    }

    if (mlist->num_nonbit) {
        sprintf(msg, "\tcheck <%ld> more motif%s", mlist->num_nonbit,
                (mlist->num_nonbit == 1) ? "" : "s");
        log_msg_vchk(msg);

        if (strand0 == 0) {
            handle_nonbit_motifs(tdat, seqs, &ntop_seeds, mlist, 1, TRUE);
            handle_nonbit_motifs(tdat, seqs, &ntop_seeds, mlist, 2, TRUE);

        } else
            handle_nonbit_motifs(tdat, seqs, &ntop_seeds, mlist, strand0, FALSE);
    }

    check_ntops_populate_seed(&ntop_seeds, topo->topos, strand0, tdat, seqs, seed, mlist);
}

static void handle_topo_psam(struct_data * tdat, struct_seqArr * seqs, struct_seed * seed,
                             long strand)
{
    long nc;
    double tval, *counts;

    struct_ls5fit fitpars;

    counts = dvector(1, seqs->num_seqs);
    get_psam_counts(seqs, seed->topo, seed->strand, seed->w, counts);

    for (nc = 1; nc <= tdat->ncol; nc++) {
        get_y_vs_x_lsfitpars(tdat->nrow, counts, tdat->resid[nc], &fitpars);
        tval = fabs(fitpars.b_tval);

        if (tval > seed->fitpars.b_tval + NR_EPS) {
            seed->col_idx = nc;
            seed->fitpars.b_tval = tval;  /* absolute value */
            seed->strand = strand;  /* the strand-ness for dynamic option */
        }
    }

    seed->bonf += tdat->ncol;

    free_dvector(counts, 1, DUMMY);
}

void find_seed_experiment(struct_data * tdat, struct_seqArr * seqs, struct_seed * seed)
{
    char stnd[BUF512], msg[BUF512];
    long strand0 = seed->strand;  /* initial setting for strandness */

    if (strand0 == 0) {  /* dynamically check leading and both strands */
        handle_topo_psam(tdat, seqs, seed, 1);
        handle_topo_psam(tdat, seqs, seed, 2);

    } else  /* follow user's command-line setting */
        handle_topo_psam(tdat, seqs, seed, strand0);

    if (seed->col_idx == 0) {
        sprintf(msg, "no proper experiment match for seed PSAM/Motif\n");
        log_msg_exit(msg);
    }

    get_psam_fitpars(tdat, seqs, seed, &seed->fitpars, FALSE);

    seed->col_name = my_strdup(tdat->col_names[seed->col_idx]);
    seed->stnd_msg = strand2msg(seed->strand);

    strcpy(stnd, (strand0 == 0) ? " [auto-determined]" : "");

    print_seed_expt(Gvars.RUNLOG, seed, stnd);
    print_seed_expt(Gvars.PRGLOG, seed, stnd);
}

long check_dictfile_repeat(long num, struct_tag * str_tags)
{
    char msg[BUF512];
    long i, nok, k = 0;

    for (i = 1; i < num; i++) {
        if (is_equal_case_string(str_tags[i].str, str_tags[i + 1].str)) {
            sprintf(msg, "\t%4ld\t%s", i, str_tags[i].str);
            log_msg(msg);
            free_cvector(str_tags[i].str, 0, DUMMY);
            str_tags[i].str = my_strdup(REPEATED_ID);  /* mask-off previous one */
            k++;
        }
    }

    if (k) {
        sprintf(msg, "\tThere are %ld motif repeats", k);
        log_msg(msg);
    }

    nok = num - k;
    if (!nok)
        fatal("\tno valid motifs in dict file\n");

    return nok;
}

void validate_expanded_pattern(char *str)
{
    char msg[BUF512], *str0;

    str0 = my_strdup(str);

    upperstr(str);
    unify_parens_chars(str);
    cvtstr_c1toc2(str, 'X', 'N');
    cvtstr_c1toc2(str, '.', '-');

    if (!has_matched_parens(str)) {
        sprintf(msg, "\tpattern <%s> contains unmatched ()", str0);
        log_msg(msg);
        strcpy(str, "");
    }

    if (!string_contains_only_those_characters(str, IUPAC_PATTERN)) {
        sprintf(msg, "\tpattern <%s> contains invalid char", str0);
        log_msg(msg);
        strcpy(str, "");
    }

    free_cvector(str0, 0, DUMMY);
}

/* tidy up dictfile: expanding abbreviations, checking for validness */
void tidy_dictfile(char *dictfile, char *cln_dictfile)
{
    char *line, *p0, str[BUF512], fullstr[BUF512];
    long num = 0;
    FILE *fp1, *fp2;

    fp1 = open_file(dictfile, "r");
    fp2 = open_file(cln_dictfile, "w");

    while ((p0 = my_getline(fp1)) != NULL) {
        line = strtok(ltrim(p0), WSPACES);  /* keep the original value of p0 */
        if (line && !is_skip_line(line)) {
            strcpy(str, line);
            expand_char_num_to_full(str, fullstr);
            validate_expanded_pattern(fullstr);
            if (!is_empty_string(fullstr)) {
                cvtstr_c1toc2(fullstr, 'U', 'T');  /* for consistency */
                fprintf(fp2, "%s\n", fullstr);
                num++;
            }
        }

        free(p0);
    }

    close_file(fp1);
    close_file(fp2);

    sprintf(str, "\t %ld valid motifs from dictionary file [%s]", num, dictfile);
    log_msg(str);
}

void read_dictfile(char *dictfile, struct_dic * mlist)
{
    char msg[BUF512], *cln_dictfile = "tidy_dictfile.dat";
    long num, nok;
    struct_tag *str_tags;

    sprintf(msg, "reading motifs from dictionary file [%s]", dictfile);
    log_msg(msg);

    tidy_dictfile(dictfile, cln_dictfile);  /* change to UPPER CASE */

    str_tags = fillup_str_tags(TRUE, cln_dictfile, &num);
    qsort(str_tags + 1, num, sizeof(struct_tag), strtags_compare);

    nok = check_dictfile_repeat(num, str_tags);

    if (nok < num) {
        char *p;
        long i, k = 0;
        struct_tag *ok_str_tags;

        ok_str_tags = allocate_memory_for_strtags(nok);
        for (i = 1; i <= num; i++) {
            p = str_tags[i].str;
            if (!is_equal_string(p, REPEATED_ID)) {
                k++;
                init_strtag(&str_tags[i], &ok_str_tags[k]);
            }
        }
        free_strtags(1, num, str_tags);

        mlist->num = nok;
        mlist->motifs = ok_str_tags;

    } else {
        mlist->num = num;
        mlist->motifs = str_tags;
    }
}

static long **get_motifs_seidx(struct_dic * mlist)
{
    long i, ib, n, num = mlist->num;
    long *temp, **seidx;

    temp = lvector(1, num);

    for (i = 1; i < num; i++)  /* by string length */
        if (strlen(mlist->motifs[i].str) != strlen(mlist->motifs[i + 1].str))
            temp[i] = 1;
    temp[num] = 1;

    n = 0;  /* get number of different length */
    for (i = 1; i <= num; i++)
        if (temp[i])
            ++n;

    seidx = lmatrix(1, n, 0, 3);  /* allocate spaces */
    seidx[1][1] = 1;

    n = 0;
    for (i = 1; i <= num; i++)
        if (temp[i])
            seidx[++n][2] = i;

    free_lvector(temp, 1, DUMMY);

    for (i = 2; i <= n; i++)
        seidx[i][1] = seidx[i - 1][2] + 1;

    for (i = 1; i <= n; i++) {
        ib = seidx[i][1];
        seidx[i][0] = strlen(mlist->motifs[ib].str);
        seidx[i][3] = seidx[i][2] - seidx[i][1] + 1;
    }

    mlist->num_seidx = n;

    return seidx;
}

void get_topo_from_topology(char *topology, struct_topo * t)
{
    long i, Xcount = 0;

    init_strtag_with_strings(topology, NULL, &t->str_tag);

    t->Tcount = strlen(topology);
    t->offset = lvector(0, t->Tcount);  /* 0-index */

    for (i = 0; i < t->Tcount; i++)
        if (topology[i] == XCHR) {
            t->offset[Xcount] = i;
            Xcount++;
        }

    t->Xcount = Xcount;
    t->num_parens = fillup_parens(topology, BUF16, t->parens_bidx, t->parens_eidx);
    t->msize = -1;
}

static void deduce_acgt_nogap_topos(struct_dic * mlist, struct_topoArr * topo)
{
    long i, k, num_topo = 0;

    struct_topo *topos;

    for (i = 1; i <= mlist->num_seidx; i++) {  /* already sorted by length */
        if (mlist->seidx[i][0] > MAXNTS_DIC)
            break;
        num_topo++;
    }

    topos = allocate_memory_for_topos(num_topo);

    for (i = 1; i <= num_topo; i++) {
        k = mlist->seidx[i][0];  /* motif length */
        get_nogap_topo(k, &topos[i]);
    }

    topo->num_topo = num_topo;
    topo->topos = topos;
}

static void write_sorted_mlist(struct_dic * mlist)
{
    char *p;
    long i, j;
    FILE *fp;

    fp = open_file("sort_dictlist.dat", "w");

    fprintf(fp, "# ------------------------------------------------------------\n");
    for (i = 1; i <= mlist->num_seidx; i++) {
        fprintf(fp, "#%ld", i);
        for (j = 0; j <= 3; j++)
            fprintf(fp, "\t%ld", mlist->seidx[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "# ------------------------------------------------------------\n");

    for (i = 1; i <= mlist->num_nonbit; i++) {
        j = mlist->idx_nonbit[i];
        p = mlist->motifs[j].str;
        fprintf(fp, "#%ld\t%s\t%ld\t%ld\n", i, p, j, (long) strlen(p));
    }
    fprintf(fp, "# ------------------------------------------------------------\n");

    for (i = 1; i <= mlist->num; i++) {
        p = mlist->motifs[i].str;
        fprintf(fp, "%s\t#%ld\t%ld\t%ld\n", p, i, (long) strlen(p), mlist->index[i]);
    }

    close_file(fp);
}

static void classify_mlist_to_2parts(struct_dic * mlist)
{
    char *p;
    long i, k = 0, num = mlist->num;

    mlist->index = lvector(1, num);

    for (i = 1; i <= num; i++) {
        p = mlist->motifs[i].str;

        if (!string_contains_only_those_characters(p, BASES) || strlen(p) > MAXNTS_DIC) {
            /* non-ACGT or long motif */
            mlist->index[i] = -1;
            k++;

        } else  /* normal ACGT motifs */
            mlist->index[i] = fast_str_q2d(p);  /* decimal representation */
    }

    mlist->num_nonbit = k;
    mlist->idx_nonbit = lvector(1, k + 1);  /* + 1 to avoid error message with k = 0 */

    k = 0;
    for (i = 1; i <= num; i++)
        if (mlist->index[i] == -1)
            mlist->idx_nonbit[++k] = i;  /* map to full list */

    write_sorted_mlist(mlist);
}

void categorize_mlist(struct_topoArr * topo, struct_dic * mlist)
{
    mlist->seidx = get_motifs_seidx(mlist);

    classify_mlist_to_2parts(mlist);

    /* simplest cases: ACGT only, no gap, within 1..MAXNTS_DIC (10) */
    deduce_acgt_nogap_topos(mlist, topo);
}
