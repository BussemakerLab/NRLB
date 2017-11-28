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

void get_full_psam(struct_seed * seed)
/* expanding seed->w[] to include non-optimized positions (N): 1.0 */
{
    long i, j, k, nW, idx = 0;
    struct_topo *t = seed->topo;  /* shorthand to simplify code */

    nW = NUM_BASE4 * t->Tcount - 1;  /* 0-indexed */

    seed->psam = dvector(0, nW);
    init_dvector(seed->psam, 0, nW, 1.0);

    for (i = 0; i < t->Xcount; i++) {
        k = t->offset[i] * NUM_BASE4;
        for (j = 0; j < NUM_BASE4; j++)
            seed->psam[k + j] = seed->w[idx++];
    }
}

char posvec_iupac(double *v0)
{
    char msg[BUF512], iupac = 'x';
    double dval, cutoff = Gvars.misc.IUPAC_CUTOFF - WEPS, v[4];
    long i, A, C, G, T;

    copy_dvector(v0, 0, TOP_BIDX3, v);  /* use a copy to keep original unchanged */

    dval = max_dvector(v, 0, TOP_BIDX3);

    if (dval <= 0) {
        sprintf(msg, "position with maximum=%g < 0: no IUPAC code, use '%c' instead",
                dval, iupac);
        log_msg(msg);

        return iupac;
    }

    if (dval != 1.0)  /* make the maximum at each position to be 1.0 */
        for (i = 0; i <= TOP_BIDX3; i++)
            v[i] /= dval;

    A = v[0] > cutoff;
    C = v[1] > cutoff;
    G = v[2] > cutoff;
    T = v[3] > cutoff;

    if (A && !C && !G && !T)
        iupac = 'A';

    else if (C && !A && !G && !T)
        iupac = 'C';

    else if (G && !A && !C && !T)
        iupac = 'G';

    else if (T && !A && !C && !G)
        iupac = (Gvars.RNA) ? 'U' : 'T';

    else if (A && T && !C && !G)  /* W --> [A, T] (weak) */
        iupac = 'W';

    else if (C && G && !A && !T)  /* S --> [C, G] (strong) */
        iupac = 'S';

    else if (A && G && !C && !T)  /* R --> [A, G] (puRine) */
        iupac = 'R';

    else if (C && T && !A && !G)  /* Y -- [C, T] (pYrimidine) */
        iupac = 'Y';

    else if (G && T && !A && !C)  /* K --> [G, T] (Keto) */
        iupac = 'K';

    else if (A && C && !G && !T)  /* M --> [A, C] (aMino) */
        iupac = 'M';

    else if (!A && C && G && T)  /* B --> [C, G, T] (not A) */
        iupac = 'B';

    else if (!C && A && G && T)  /* D --> [A, G, T] (not C) */
        iupac = 'D';

    else if (!G && A && C && T)  /* H --> [A, C, T] (not G) */
        iupac = 'H';

    else if (!T && A && C && G)  /* V --> [A, C, G] (not U) */
        iupac = 'V';

    else if (A && C && G && T)  /* N --> [A, C, G, T] (aNy) */
        iupac = 'N';

    else
        log_msg("this should not happy!");

    return iupac;
}

void get_psam_optimal_seq(double *psam, long num, char *optimal)
{
    long i, k = 0;

    for (i = 0; i < num; i++) {
        optimal[i] = posvec_iupac(&psam[k]);
        k += NUM_BASE4;
    }

    optimal[i] = '\0';
}

void write_psam_data(FILE * fp, long Tcount, char *opt_seq, char *full, double *psam)
{
    int c, X;
    long i, j, idx = 0;

    fprintf(fp, "\n<psam>\n");
    fprintf(fp, "# A            C            G            T             # no. opt  \n");
    fprintf(fp, "# +============+============+============+============ # ==+===+==\n");

    for (i = 0; i < Tcount; i++) {
        fprintf(fp, "  ");
        for (j = 0; j < NUM_BASE4; j++)
            fprintf(fp, "%-12.6g ", psam[idx++]);

        if (opt_seq && full) {
            X = toupper((int) opt_seq[i]);
            fprintf(fp, " # %3ld  %2c", i + 1, X);
            c = toupper((int) full[i]);
            fprintf(fp, " %c", (X == c || strchr("-()", c)) ? ' ' : 'x');
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "</psam>\n");
}

void write_meta(FILE * fp, struct_seed * seed, char *measfile, double p_value)
{
    char str[BUF512];

    get_currentTimeString(str);

    fprintf(fp, "<meta>\n");
    fprintf(fp, "    <source>%s</source>\n", VERSION_NUMBER);
    fprintf(fp, "    <comment>%s</comment>\n", "xl2134@columbia.edu");
    fprintf(fp, "    <date>%s</date>\n", str);

    fprintf(fp, "    <topology>%s</topology>\n", seed->topo->str_tag.str);
    fprintf(fp, "    <seed_motif>%s</seed_motif>\n", seed->full);

    fprintf(fp, "    <measurement_file>%s</measurement_file>\n", basename(measfile));
    fprintf(fp, "    <experiment_name>%s</experiment_name>\n", seed->col_name);
    fprintf(fp, "    <experiment_column>%ld</experiment_column>\n", seed->col_idx);

    fprintf(fp, "    <bonferroni>%ld</bonferroni>\n", seed->bonf);
    fprintf(fp, "    <p_value>%g</p_value>\n", p_value);
    fprintf(fp, "</meta>\n");
}

void print_psam(struct_seed * seed, char *measfile, char *outdir, double p_value,
                long motif_num)
/* direct output from MatrixREDUCE/MotifREDUCE/OptimizePSAM */
{
    char psam_file[BUF512];
    long Tcount = seed->topo->Tcount;
    FILE *fp;

    get_full_psam(seed);

    sprintf(psam_file, "%s/psam_%3.3ld.xml", outdir, motif_num);
    fp = open_file(psam_file, "w");

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(fp, "<matrix_reduce>\n");
    write_meta(fp, seed, measfile, p_value);

    fprintf(fp, "\n<directionality>%s</directionality>\n", seed->stnd_msg);
    fprintf(fp, "<psam_length>%ld</psam_length>\n", Tcount);

    if (seed->optimal)  /* e.g., initialized in recover_seed_from_psam() */
        free_cvector(seed->optimal, 0, DUMMY);
    seed->optimal = my_strdup(seed->full);
    get_psam_optimal_seq(seed->psam, Tcount, seed->optimal);
    fprintf(fp, "<optimal_sequence>%s</optimal_sequence>\n", seed->optimal);

    write_psam_data(fp, Tcount, seed->optimal, seed->full, seed->psam);
    fprintf(fp, "</matrix_reduce>\n");

    close_file(fp);
}

void get_tag_string_pair(char *tag, char *btag, char *etag)
{
    sprintf(btag, "<%s>", tag);
    sprintf(etag, "</%s>", tag);
}

void extract_xml_line_string(FILE * fp, char *tag, char *strval, long required)
{
    char *p, *bp, *ep, btag[BUF512], etag[BUF512], temp[BUF512];

    get_tag_string_pair(tag, btag, etag);
    strcpy(strval, "");  /* default to empty string */

    rewind(fp);

    while ((p = my_getline(fp)) != NULL) {
        bp = strstr(p, btag);
        ep = strstr(p, etag);
        if (bp && ep) {
            *ep = '\0';
            strcpy(temp, bp + strlen(btag));
            strcpy(strval, trim(temp));
            free(p);
            return;
        }
        free(p);
    }

    if (required)
        fatal("required tag pair <%s> not found\n", tag);
}

long extract_xml_line_long(FILE * fp, char *tag, long dft_value, long required)
{
    char strval[BUF512] = "";

    extract_xml_line_string(fp, tag, strval, required);

    return is_empty_string(strval) ? dft_value : cvt2long(strval);
}

double extract_xml_line_double(FILE * fp, char *tag, double dft_value, long required)
{
    char strval[BUF512] = "";

    extract_xml_line_string(fp, tag, strval, required);

    return is_empty_string(strval) ? dft_value : cvt2double(strval);
}

void extract_xml_psam(FILE * fp, long npos, double *psam)
{
    char *p0, *p1, *line, btag[BUF512], etag[BUF512], *items[BUF512];
    long i, nitem, idx = 0, num = 0;

    get_tag_string_pair("psam", btag, etag);

    rewind(fp);

    while ((p0 = my_getline(fp)) != NULL) {
        if (strstr(p0, btag)) {
            free(p0);
            while ((p1 = my_getline(fp)) != NULL) {
                if (strstr(p1, etag)) {
                    free(p1);
                    assert_equal_longs("Check PSAM length", num, npos);
                    return;
                }

                line = trim(p1);
                if (!is_skip_line(line)) {
                    num++;
                    if (num > npos)
                        fatal("PSAM tag position size > PSAM_LENGTH (%ld)\n", npos);

                    nullify_line_comment(line);

                    line = trim(p1);
                    nitem = item_list(line, items, NUM_BASE4, WSPACES);
                    if (nitem != NUM_BASE4)
                        fatal("wrong number of Ws [%ld <> %ld] for position #%ld\n",
                              nitem, NUM_BASE4, num);
                    for (i = 1; i <= NUM_BASE4; i++)
                        psam[idx++] = cvt2double(items[i]);
                }
                free(p1);
            }
        } else
            free(p0);
    }

    fatal("required tag <psam> not found\n");
}

/* matched ()s and gap '-' as is; other IUPAC as X positions: */
void take_expanded_motif_as_topology(char *motif, char *topology)
{
    long i;

    validate_expanded_pattern(motif);
    if (is_empty_string(motif))
        fatal("Invalid motif pattern: %s\n", motif);

    strcpy(topology, motif);
    for (i = 0; i < (long) strlen(motif); i++)
        if (!strchr("-()", motif[i]))
            topology[i] = XCHR;
}

/* to all XCHR, i.e., all positions to be optimized */
void set_default_topology(char *topology, long num)
{
    long i;

    for (i = 0; i < num; i++)
        topology[i] = XCHR;

    topology[i] = '\0';
}

void check_psam_value(long num, double *psam)
{
    char msg[BUF512];
    long i, num_bad = 0, idx = 0;
    double Wmin, Wmax;

    for (i = 0; i < num; i++) {
        Wmax = get_position_Wmax(&psam[idx]);

        if (Wmax != 1.0) {
            sprintf(msg, "\t\t>>position #%ld, maximum W = %g<<", i + 1, Wmax);
            log_msg(msg);
            num_bad++;
        }

        idx += NUM_BASE4;
    }

    if (num_bad) {
        sprintf(msg, "\t%ld position%s with Wmax <> 1.0: not a compliant PSAM",
                num_bad, (num_bad == 1) ? "" : "s");
        log_msg(msg);
    }

    Wmin = min_dvector(psam, 0, num * NUM_BASE4 - 1);
    if (Wmin < 0.0) {
        sprintf(msg, "\tWmin=%g < 0.0: not a compliant PSAM", Wmin);
        log_msg(msg);
    }
}

void check_psam_topology(long num, char *topology)
{
    char msg[BUF512], *valid_topology_chars = "X-()";
    long topology_len, isOk = TRUE;

    if (is_empty_string(topology)) {
        set_default_topology(topology, num);  /* all Xs */
        return;
    }

    topology_len = strlen(topology);
    if (topology_len != num) {
        sprintf(msg, "topology [%s] length (%ld) does NOT match length of PSAM"
                " (%ld)", topology, topology_len, num);
        log_msg(msg);
        isOk = FALSE;
    }

    if (!string_contains_only_those_characters(topology, valid_topology_chars)) {
        sprintf(msg, "topology [%s] contains non [%s] characters", topology,
                valid_topology_chars);
        log_msg(msg);
        isOk = FALSE;
    }

    if (!isOk) {
        sprintf(msg, "\treset topology [%s] to the default: all Xs", topology);
        log_msg(msg);
        set_default_topology(topology, num);
        sprintf(msg, "\t\t[%s] (%ld)", topology, num);
        log_msg(msg);
    }
}

void check_psam_optimal_seq(long num, double *psam, char *optseq)
{
    char msg[BUF512], derived_optseq[BUF512], valid_chars[BUF512] = "NACGT";
    long i, j, k, m = -1, isN, idx = 0;
    double Wmax;

    for (i = 0; i < num; i++) {
        Wmax = -XBIG;
        isN = TRUE;

        for (j = 0; j < NUM_BASE4; j++) {
            k = idx + j;

            if (Wmax < psam[k]) {  /* A preferred over C > G > T with equal W */
                Wmax = psam[k];
                m = j;
            }

            if (psam[k] != psam[idx])  /* any W is different */
                isN = FALSE;
        }

        derived_optseq[i] = isN ? NCHR : BASES[m];

        idx += NUM_BASE4;  /* move to the next position */
    }

    derived_optseq[i] = '\0';

    if (Gvars.RNA) {
        convert_T_to_U(optseq);
        convert_T_to_U(derived_optseq);
        strcpy(valid_chars, "NACGU");
    }

    if (is_empty_string(optseq))
        strcpy(optseq, derived_optseq);

    else if (!is_equal_string(optseq, derived_optseq) &&  /* be cautious re this message */
             string_contains_only_those_characters(optseq, valid_chars)) {
        sprintf(msg, "\t*** OPTIMAL_SEQUENCE [%s] does not match derived one [%s]",
                optseq, derived_optseq);
        log_msg(msg);
    }
}

void write_pwm_entry(char *pwm_file, struct_psam * pd)
/* 'standard' PWM entry: >id followed by an N-by-4 matrix */
{
    double *mtx = pd->psam;
    long i, j, npos = pd->psam_length, idx = 0;
    FILE *fp;

    fp = open_file(pwm_file, "w");

    fprintf(fp, ">%s\n", pd->source);
    for (i = 0; i < npos; i++) {
        for (j = 0; j < NUM_BASE4; j++)
            fprintf(fp, " %11.6g", mtx[idx++]);
        fprintf(fp, "\n");
    }
    close_file(fp);
}

void write_psam_struct(char *psam_file, struct_psam * pd)
/* more generic PSAM struct, e.g., converted from other sources */
{
    char *stnd_msg, str[BUF512];
    FILE *fp;

    fp = open_file(psam_file, "w");

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(fp, "<matrix_reduce>\n");
    get_currentTimeString(str);

    fprintf(fp, "<meta>\n");
    fprintf(fp, "    <source>%s</source>\n", pd->source);
    fprintf(fp, "    <comment>%s</comment>\n", "xl2134@columbia.edu");
    fprintf(fp, "    <date>%s</date>\n", str);

    fprintf(fp, "    <topology>%s</topology>\n", pd->topology);
    fprintf(fp, "    <seed_motif>%s</seed_motif>\n", pd->seed_motif);

    if (!is_empty_string(pd->measurement_file))
        fprintf(fp, "    <measurement_file>%s</measurement_file>\n",
                pd->measurement_file);
    if (!is_empty_string(pd->experiment_name))
        fprintf(fp, "    <experiment_name>%s</experiment_name>\n", pd->experiment_name);
    if (pd->experiment_column != -1)
        fprintf(fp, "    <experiment_column>%ld</experiment_column>\n",
                pd->experiment_column);

    if (pd->bonferroni != 1)
        fprintf(fp, "    <bonferroni>%ld</bonferroni>\n", pd->bonferroni);

    if (pd->p_value < 1.0)
        fprintf(fp, "    <p_value>%g</p_value>\n", pd->p_value);
    fprintf(fp, "</meta>\n");

    stnd_msg = strand2msg(pd->strand);
    fprintf(fp, "\n<directionality>%s</directionality>\n", stnd_msg);
    fprintf(fp, "<psam_length>%ld</psam_length>\n", pd->psam_length);
    fprintf(fp, "<optimal_sequence>%s</optimal_sequence>\n", pd->optimal_sequence);

    write_psam_data(fp, pd->psam_length, pd->optimal_sequence, pd->seed_motif, pd->psam);
    fprintf(fp, "</matrix_reduce>\n");

    close_file(fp);

    free_cvector(stnd_msg, 0, DUMMY);
}

void read_psam_struct(char *psam_file, struct_psam * pd)
{
    char str[BUF512];
    FILE *fp;

    initialize_psam_struct(pd);

    fp = open_file(psam_file, "r");

    extract_xml_line_string(fp, "source", pd->source, FALSE);
    extract_xml_line_string(fp, "topology", pd->topology, FALSE);
    extract_xml_line_string(fp, "seed_motif", pd->seed_motif, FALSE);
    extract_xml_line_string(fp, "measurement_file", pd->measurement_file, FALSE);
    extract_xml_line_string(fp, "experiment_name", pd->experiment_name, FALSE);
    pd->experiment_column = extract_xml_line_long(fp, "experiment_column", -1, FALSE);
    pd->bonferroni = extract_xml_line_long(fp, "bonferroni", 1, FALSE);
    pd->p_value = extract_xml_line_double(fp, "p_value", 1.0, FALSE);

    extract_xml_line_string(fp, "directionality", str, FALSE);
    pd->strand = msg2strand(str);

    pd->psam_length = extract_xml_line_long(fp, "psam_length", -1, TRUE);
    if (pd->psam_length == -1)
        fatal("\trequired tag <psam_length> not specified");
    extract_xml_line_string(fp, "optimal_sequence", pd->optimal_sequence, FALSE);

    pd->psam = dvector(0, NUM_BASE4 * pd->psam_length - 1);
    extract_xml_psam(fp, pd->psam_length, pd->psam);

    check_psam_value(pd->psam_length, pd->psam);
    check_psam_topology(pd->psam_length, pd->topology);
    check_psam_optimal_seq(pd->psam_length, pd->psam, pd->optimal_sequence);

    if (is_empty_string(pd->seed_motif))  /* no <seed_motif>...</seed_motif> */
        strcpy(pd->seed_motif, pd->optimal_sequence);

    close_file(fp);
}

void write_psam_struct_lite(char *psam_file, struct_psam * pd)
{
    char str[BUF512];
    FILE *fp;

    fp = open_file(psam_file, "w");

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(fp, "<matrix_reduce>\n");
    get_currentTimeString(str);

    fprintf(fp, "<meta>\n");
    fprintf(fp, "    <comment>%s</comment>\n", "xl2134@columbia.edu");
    fprintf(fp, "    <date>%s</date>\n", str);
    fprintf(fp, "</meta>\n");

    fprintf(fp, "\n<psam_length>%ld</psam_length>\n", pd->psam_length);

    write_psam_data(fp, pd->psam_length, pd->optimal_sequence, pd->seed_motif, pd->psam);
    fprintf(fp, "</matrix_reduce>\n");

    close_file(fp);
}

void read_psam_struct_lite(char *psam_file, struct_psam * pd)
{
    FILE *fp;

    initialize_psam_struct(pd);

    fp = open_file(psam_file, "r");

    pd->psam_length = extract_xml_line_long(fp, "psam_length", -1, TRUE);
    if (pd->psam_length == -1)
        fatal("\trequired tag <psam_length> not specified");

    pd->psam = dvector(0, NUM_BASE4 * pd->psam_length - 1);
    extract_xml_psam(fp, pd->psam_length, pd->psam);

    close_file(fp);
}

/* the w[] contains only Xcount * 4 Ws, as in seed motif */
void populate_w_from_psam(double *psam, struct_seed * seed)
{
    long i, idx, j, k = 0;
    struct_topo *t = seed->topo;

    for (i = 0; i < t->Xcount; i++) {
        idx = NUM_BASE4 * t->offset[i];
        for (j = 0; j < NUM_BASE4; j++)
            seed->w[k + j] = psam[idx + j];
        k += NUM_BASE4;
    }
}

void expanded_motif_to_seedtopo(char *motif, struct_seed * seed)
{
    char topology[BUF512];

    take_expanded_motif_as_topology(motif, topology);

    if (seed->topo) {
        char msg[BUF512];

        sprintf(msg, "\tmotif '%s' seed already has topo struct.", motif);
        log_msg(msg);

        free_one_topo(seed->topo);
    }

    seed->topo = allocate_memory_for_one_topo();
    get_topo_from_topology(topology, seed->topo);
}

void expanded_motif_to_seed(char *motif, long strand, struct_seed * seed)
{
    double *psam;
    long nw, num_nt = strlen(motif);

    seed->strand = strand;
    init_seed_motif(seed, strand, 1.0);

    expanded_motif_to_seedtopo(motif, seed);

    psam = dvector(0, NUM_BASE4 * num_nt - 1);
    iupac2psam(motif, psam);

    nw = NUM_BASE4 * seed->topo->Xcount - 1;
    if (nw > 0) {  /* with Xcount */
        seed->w = dvector(0, nw);
        populate_w_from_psam(psam, seed);
    }

    free_dvector(psam, 0, DUMMY);
}

void get_motif_counts(char *motif, long strand, struct_seqArr * seqs, double *counts)
{
    struct_seed sx;

    expanded_motif_to_seed(motif, strand, &sx);

    get_psam_counts(seqs, sx.topo, sx.strand, sx.w, counts);

    free_seed(&sx);
}

void motif2psam_struct(char *motif, long strand, double p_value, struct_psam * psam_data)
{
    char str0[BUF512];
    long num_nt;

    strcpy(str0, motif);
    expand_char_num_to_full(str0, motif);
    validate_expanded_pattern(motif);

    if (is_empty_string(motif))
        fatal("invalid motif: %s [%s]\n", str0, motif);

    num_nt = strlen(motif);  /* based on extended form */

    initialize_psam_struct(psam_data);
    sprintf(psam_data->source, "Seeded with motif: %s", motif);

    take_expanded_motif_as_topology(motif, psam_data->topology);
    strcpy(psam_data->seed_motif, motif);

    psam_data->p_value = p_value;

    psam_data->strand = strand;
    psam_data->psam_length = num_nt;
    strcpy(psam_data->optimal_sequence, motif);

    psam_data->psam = dvector(0, NUM_BASE4 * num_nt - 1);
    iupac2psam(motif, psam_data->psam);
}

void recover_seed_from_psam(struct_psam * pd, struct_seed * seed)
{
    long num;

    init_seed_motif(seed, pd->strand, pd->p_value);

    seed->topo = allocate_memory_for_one_topo();
    seed->topo_idx = DUMMY;
    get_topo_from_topology(pd->topology, seed->topo);

    seed->motif = my_strdup(pd->seed_motif);  /* full, i.e., with gaps */
    seed->full = my_strdup(pd->seed_motif);
    seed->optimal = my_strdup(pd->optimal_sequence);

    num = NUM_BASE4 * seed->topo->Xcount - 1;
    if (num > 0) {  /* with Xcount */
        seed->w = dvector(0, num);
        populate_w_from_psam(pd->psam, seed);
    }
}

void populate_motif_psam(struct_seed * seed)
{
    long num_pars;

    if (seed->topo->Xcount == 0)
        return;

    num_pars = NUM_BASE4 * seed->topo->Xcount + 2;  /* plus C & F */
    seed->x = dvector(1, num_pars);  /* 1-index, per NRC convention */

    init_lmseed(num_pars, 0.0, seed);

    seed->w = dvector(0, num_pars - 3);  /* only Ws, excluding C & F */
    normalize_psam(num_pars, seed->x, seed->w);
}

void rc_psam(struct_psam * psam_data)
/* get the reverse complementary PSAM representation in place */
{
    long i, j, k, idx = 0, nW = NUM_BASE4 * psam_data->psam_length;
    double *temp;

    temp = dvector(0, nW - 1);

    for (i = 0; i < psam_data->psam_length; i++) {  /* get complementary first */
        for (j = 0; j < NUM_BASE4; j++)
            temp[idx + j] = psam_data->psam[idx + TOP_BIDX3 - j];
        idx += NUM_BASE4;
    }

    idx = 0;
    for (i = 0; i < psam_data->psam_length; i++) {  /* then reverse it */
        k = (psam_data->psam_length - 1 - i) * NUM_BASE4;
        for (j = 0; j < NUM_BASE4; j++)
            psam_data->psam[k + j] = temp[idx + j];
        idx += NUM_BASE4;
    }

    free_dvector(temp, 0, DUMMY);

    reverse_cmpl(psam_data->optimal_sequence, FALSE);
    reverse_cmpl(psam_data->seed_motif, FALSE);
    reverse_cmpl(psam_data->topology, FALSE);
}

void iupac2psam(char *motif, double *psam)
/* IUPAC to PSAM conversion: assuming psam[] initialized to 0s */
{
    size_t i, k = 0;
    enum { A = 0, C = 1, G = 2, T = 3 };

    for (i = 0; i < strlen(motif); i++) {
        switch (toupper((int) motif[i])) {
        case 'A':
            psam[k + A] = 1.0;
            break;
        case 'C':
            psam[k + C] = 1.0;
            break;
        case 'G':
            psam[k + G] = 1.0;
            break;
        case 'T':
        case 'U':
            psam[k + T] = 1.0;
            break;
        case 'W':  /* W --> [A, T] (weak) */
            psam[k + A] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'S':  /* S --> [C, G] (strong) */
            psam[k + C] = 1.0;
            psam[k + G] = 1.0;
            break;
        case 'R':  /* R --> [A, G] (puRine) */
            psam[k + A] = 1.0;
            psam[k + G] = 1.0;
            break;
        case 'Y':  /* Y -- [C, T] (pYrimidine) */
            psam[k + C] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'K':  /* K --> [G, T] (Keto) */
            psam[k + G] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'M':  /* M --> [A, C] (aMino) */
            psam[k + A] = 1.0;
            psam[k + C] = 1.0;
            break;
        case 'B':  /* B --> [C, G, T] (not A) */
            psam[k + C] = 1.0;
            psam[k + G] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'D':  /* D --> [A, G, T] (not C) */
            psam[k + A] = 1.0;
            psam[k + G] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'H':  /* H --> [A, C, T] (not G) */
            psam[k + A] = 1.0;
            psam[k + C] = 1.0;
            psam[k + T] = 1.0;
            break;
        case 'V':  /* V --> [A, C, G] (not U) */
            psam[k + A] = 1.0;
            psam[k + C] = 1.0;
            psam[k + G] = 1.0;
            break;
        case 'N':  /* N --> [A, C, G, T] (aNy) */
        case '(':
        case ')':
        case '-':
            psam[k + A] = 1.0;
            psam[k + C] = 1.0;
            psam[k + G] = 1.0;
            psam[k + T] = 1.0;
            break;
        default:
            fatal("illegal base %c -- not in IUPAC set %s\n", motif[i], IUPAC_BASES);
            break;
        }

        k += NUM_BASE4;
    }
}

/* check if motif contains only valid IUPAC characters */
long is_valid_iupac_motif(char *motif, long debug)
{
    char msg[BUF512];
    size_t i;

    for (i = 0; i < strlen(motif); i++) {
        if (strchr(IUPAC_BASES, toupper((int) motif[i])) == NULL) {
            if (debug) {
                sprintf(msg, "\t<%s> contains non-IUPAC base character '%c'",
                        motif, motif[i]);
                log_msg(msg);
            }
            return FALSE;
        }
    }

    return TRUE;
}

void normalize_psam(long ma, double *a, double *psam)
/* psam[] is 0-indexed, following 0123/ACGT code */
{
    long i, nW = ma - 2;  /* number of Ws */

    for (i = 1; i <= nW; i++)  /* a[] is 1-indexed */
        psam[i - 1] = a[i] * a[i];

    /* commented out following Harmen's suggestion */
#if LUXDBG == TRUE
    {
        char msg[BUF512];

        sprintf(msg, "\tC=%g\tF=%g\n", a[ma - 1], a[ma]);
        log_msg_vchk(msg);
    }
#endif

    normalize_each_psam_position(nW / NUM_BASE4, psam, TRUE);
}

void normalize_each_psam_position(long num, double *psam, long debug)
{
    char msg[BUF512];
    long i, j, k, idx = 0;
    double Wmax;

    for (i = 0; i < num; i++) {
        Wmax = get_position_Wmax(&psam[idx]);

        if (Wmax != 1.0) {
            if (debug) {
                sprintf(msg, "\t\t>>position #%ld, maximum W = %g<<", i + 1, Wmax);
                log_msg(msg);
            }

            for (j = 0; j < NUM_BASE4; j++) {
                k = idx + j;
                psam[k] /= Wmax;
            }
        }

        idx += NUM_BASE4;  /* move to the next position */
    }
}
