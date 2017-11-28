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

struct_win *allocate_memory_for_windows(long num)
{
    struct_win *windows;

    windows = (struct_win *) malloc((num + 1) * sizeof(struct_win));

    if (windows == NULL)
        fatal("malloc failure for allocate_memory_for_windows()\n");

    return windows;
}

void free_windows(struct_winArr * wins)
{
    long i;
    struct_win *w;

    for (i = 1; i <= wins->num_wins; i++) {
        w = &wins->windows[i];

        if (w->seqid)
            free_cvector(w->seqid, 0, DUMMY);

        if (w->winid)
            free_cvector(w->winid, 0, DUMMY);
    }

    free(wins->windows);
}

struct_id *allocate_memory_for_idstr(long num)
{
    struct_id *idstrs;

    idstrs = (struct_id *) malloc((num + 1) * sizeof(struct_id));

    if (idstrs == NULL)
        fatal("malloc failure for allocate_memory_for_idstr()\n");

    return idstrs;
}

void free_idstr(long ib, long ie, struct_id * idstrs)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (idstrs[i].idstr)
            free_cvector(idstrs[i].idstr, 0, DUMMY);

    free(idstrs);
}

struct_tag *allocate_memory_for_strtags(long num)
{
    struct_tag *str_tags;

    str_tags = (struct_tag *) malloc((num + 1) * sizeof(struct_tag));

    if (str_tags == NULL)
        fatal("malloc failure for allocate_memory_for_strtags()\n");

    return str_tags;
}

void free_strtags(long ib, long ie, struct_tag * str_tags)
{
    free_extra_strtags(ib, ie, str_tags);
    free(str_tags);
}

void free_p_strtag(struct_tag * str_tag)
{
    free_cvector(str_tag->str, 0, DUMMY);
    free_cvector(str_tag->tag, 0, DUMMY);
}

void free_extra_strtags(long ib, long ie, struct_tag * str_tags)
{
    long i;

    for (i = ib; i <= ie; i++)
        free_p_strtag(&str_tags[i]);
}

char **allocate_char_marray(long num)
{
    long i;
    char **pp;

    pp = (char **) malloc((num + 1) * sizeof(char *));

    if (pp == NULL)
        fatal("malloc failure for allocate_char_marray()\n");

    for (i = 0; i <= num; i++)
        pp[i] = NULL;

    return pp;
}

void free_char_marray(long ib, long ie, char **pp)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (pp[i])
            free_cvector(pp[i], 0, DUMMY);

    free(pp);
}

long **allocate_long_marray(long num)
{
    long i, **pp;

    pp = (long **) malloc((num + 1) * sizeof(long *));

    if (pp == NULL)
        fatal("malloc failure for allocate_long_marray()\n");

    for (i = 0; i <= num; i++)
        pp[i] = NULL;

    return pp;
}

void free_long_marray(long ib, long ie, long **pp)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (pp[i])
            free_lvector(pp[i], 1, DUMMY);

    free(pp);
}

double **allocate_double_marray(long num)
{
    double **pp;
    long i;

    pp = (double **) malloc((num + 1) * sizeof(double *));

    if (pp == NULL)
        fatal("malloc failure for allocate_double_marray()\n");

    for (i = 0; i <= num; i++)
        pp[i] = NULL;

    return pp;
}

void free_double_marray(long ib, long ie, double **pp)
{
    long i;

    for (i = ib; i <= ie; i++)
        if (pp[i])
            free_dvector(pp[i], 1, DUMMY);

    free(pp);
}

struct_ntop *allocate_memory_for_ntop_entries(long num)
{
    long i;
    struct_ntop *ntops;

    ntops = (struct_ntop *) malloc((num + 1) * sizeof(struct_ntop));

    if (ntops == NULL)
        fatal("malloc failure for allocate_memory_for_ntop_entries()\n");

    for (i = 0; i <= num; i++) {  /* note: 0..num; 0 for the selected */
        ntops[i].col_idx = 0;
        ntops[i].topo_idx = DUMMY;
        ntops[i].index = -1;
        ntops[i].strand = UNSET_LVAL;
        ntops[i].b_tval = 0.0;
    }

    return ntops;
}

void free_ntop_entries(struct_ntop * ntops)
{
    free(ntops);
}

struct_seed *allocate_memory_for_smotifs(long num)
{
    struct_seed *smotifs;

    smotifs = (struct_seed *) malloc((num + 1) * sizeof(struct_seed));

    if (smotifs == NULL)
        fatal("malloc failure for allocate_memory_for_smotifs()\n");

    return smotifs;
}

void free_smotifs(long ib, long ie, struct_seed * smotifs)
{
    long i;

    for (i = ib; i <= ie; i++)
        free_seed(&smotifs[i]);

    free(smotifs);
}

/* The smotifs[] is allocated via allocate_memory_for_smotifs(). Its
 * component is of type struct_seed, not a pointer to it: no free() */
void free_seed(struct_seed * seed)
{
    if (seed) {
        if (seed->col_name)
            free_cvector(seed->col_name, 0, DUMMY);

        if (seed->motif)
            free_cvector(seed->motif, 0, DUMMY);
        if (seed->full)
            free_cvector(seed->full, 0, DUMMY);
        if (seed->optimal)
            free_cvector(seed->optimal, 0, DUMMY);
        if (seed->stnd_msg)
            free_cvector(seed->stnd_msg, 0, DUMMY);

        if (seed->x)
            free_dvector(seed->x, 1, DUMMY);
        if (seed->w)
            free_dvector(seed->w, 0, DUMMY);
        if (seed->psam)
            free_dvector(seed->psam, 0, DUMMY);

        if (seed->topo)
            free_one_topo(seed->topo);
    }
}

struct_fitpars *allocate_memory_for_lfpars(long num, long vsize)
{
    long i;
    struct_fitpars *lfpars;

    lfpars = (struct_fitpars *) malloc((num + 1) * sizeof(struct_fitpars));

    if (lfpars == NULL)
        fatal("malloc failure for allocate_memory_for_lfpars()\n");

    for (i = 1; i <= num; i++) {
        lfpars[i].fval = dvector(0, vsize);  /* 0 for intercept */
        lfpars[i].tval = dvector(0, vsize);
        lfpars[i].pval = dvector(0, vsize);
    }

    return lfpars;
}

void free_lfpars(long ib, long ie, struct_fitpars * lfpars)
{
    long i;

    for (i = ib; i <= ie; i++) {
        free_dvector(lfpars[i].fval, 0, DUMMY);
        free_dvector(lfpars[i].tval, 0, DUMMY);
        free_dvector(lfpars[i].pval, 0, DUMMY);
    }

    free(lfpars);
}

struct_ls5fit *allocate_memory_for_fitpars(long num)
{
    struct_ls5fit *fitpars;

    fitpars = (struct_ls5fit *) malloc((num + 1) * sizeof(struct_ls5fit));

    if (fitpars == NULL)
        fatal("malloc failure for allocate_memory_for_fitpars()\n");

    return fitpars;
}

void free_fitpars(struct_ls5fit * fitpars)
{
    free(fitpars);
}

void init_mlist(struct_dic * mlist)
{
    mlist->num = 0;
    mlist->index = NULL;

    mlist->num_nonbit = 0;
    mlist->idx_nonbit = NULL;

    mlist->num_seidx = 0;
    mlist->seidx = NULL;

    mlist->motifs = NULL;
}

void free_mlist(struct_dic * mlist)
{
    if (mlist->index)
        free_lvector(mlist->index, 1, DUMMY);

    if (mlist->idx_nonbit)
        free_lvector(mlist->idx_nonbit, 1, DUMMY);

    if (mlist->seidx)
        free_lmatrix(mlist->seidx, 1, DUMMY, 0, DUMMY);

    if (mlist->motifs)
        free_strtags(1, mlist->num, mlist->motifs);
}

void free_resources(long psam_num, double **psam_counts, struct_seed * smotifs,
                    struct_data * tdat, struct_seqArr * seqs, struct_topoArr * topo,
                    struct_fitpars * lfpars, struct_dic * mlist)
{
    free_double_marray(1, psam_num, psam_counts);
    free_smotifs(1, psam_num, smotifs);
    free_tdt_data(tdat);
    free_sequences(NULL, seqs);
    free_topos(topo);
    free_lfpars(1, tdat->ncol, lfpars);
    free_mlist(mlist);
}

void initialize_psam_struct(struct_psam * psam_data)
{
    strcpy(psam_data->source, "");
    strcpy(psam_data->topology, "");
    strcpy(psam_data->seed_motif, "");
    strcpy(psam_data->measurement_file, "");
    strcpy(psam_data->experiment_name, "");
    psam_data->experiment_column = -1;
    psam_data->bonferroni = 1;
    psam_data->p_value = 1.0;

    psam_data->strand = 0;
    psam_data->psam_length = -1;  /* required field */
    strcpy(psam_data->optimal_sequence, "");

    psam_data->psam = NULL;
}

void free_psam_data(struct_psam * psam_data)
{
    if (psam_data->psam)
        free_dvector(psam_data->psam, 0, DUMMY);
}

struct_logo_info *allocate_memory_for_logoData(long num)
{
    struct_logo_info *logoData;

    logoData = (struct_logo_info *) malloc((num + 1) * sizeof(struct_logo_info));

    if (logoData == NULL)
        fatal("malloc failure for allocate_memory_for_logoData()\n");

    return logoData;
}
