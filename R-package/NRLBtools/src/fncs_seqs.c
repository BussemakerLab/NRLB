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

void read_sequences(char *seqfile, struct_seqArr * seqs, long raw_seq, long sort_seq)
/* [1] header line starts with '>' at the first column
 * [2] empty lines or lines starting with '#' or ';' are ignored
 * [3] each sequence line is terminated by a white-space */
{
    char msg[BUF512], *line, *p0, *pb = NULL;
    long nb = 0, num, maxline;
    FILE *fp;
    struct_seq *sequences;  /* 1-index */

    sprintf(msg, "reading in sequences from file [%s]", seqfile);
    log_msg(msg);

    num = get_number_of_fasta_entries(seqfile);
    sequences = allocate_memory_for_sequences(num);

    fp = open_file(seqfile, "r");
    num = 0;  /* re-initialized */

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original address value of p0 */
        if (is_fasta_skip_line(line)) {
            free(p0);
            continue;  /* skip empty and commented lines [; #] */
        }

        if (line[0] == '>') {  /* header line */
            assign_sequence(num, pb, nb, sequences, raw_seq);  /* for previous sequence */

            num++;
            extract_id_header(num, line, &sequences[num]);

            maxline = BUF1024;  /* reset starting size */
            pb = (char *) malloc(maxline * sizeof(char));  /* initial size */
            if (pb == NULL)
                fatal("malloc failure in reading sequence\n");
            nb = 0;  /* reset count for next sequence entry */

        } else if (num) {  /* now comes the sequence itself */
            while (*line && !isspace((int) *line)) {  /* any char upto the 1st \w */
                if (nb >= maxline - 1)
                    pb = enlarge_cline(&maxline, pb);
                pb[nb++] = *line++;
            }
        }

        free(p0);  /* release memory */
    }

    assign_sequence(num, pb, nb, sequences, raw_seq);  /* for the last sequence */
    close_file(fp);

    if (!num)
        fatal("no valid sequence in file [%s]\n", seqfile);

    if (sort_seq) {
        qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);  /* 1-based index */
        check_redundancy(seqfile, num, sequences);  /* mask repeat id with REPEATED_ID */
        qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);  /* 1-based index */
    }

    seqs->num_seqs = num;
    seqs->sequences = sequences;
}

void extract_id_header(long num, char *line, struct_seq * seq)
{
    char *p, *p2, str[BUF512];
    long n = strlen(line);

    /* make a copy of the original title line w/o leading '>' */
    seq->title = my_strdup(line + 1);

    p = strtok(line, WSPACES);  /* get sequence ID */
    if (is_equal_string(p, ">")) {
        if (n == 1) {  /* all empty, i.e., '>' */
            sprintf(str, "id_%ld", num);
            seq->id = my_strdup(str);
        } else {  /* in case like: '> BAHG_VITSP' */
            p2 = trim(p + 2);
            p = strtok(p2, WSPACES);
            seq->id = my_strdup(p);
        }

    } else  /* normal case: e.g. '>S13421' */
        seq->id = my_strdup(p + 1);  /* excluding the leading '>' */
}

void assign_sequence(long num, char *pb, long nb, struct_seq * sequences, long raw_seq)
{
    char *p0, **seq_fragments;
    long i, k, num_fragments = 0;

    if (!num)  /* no header line read in yet */
        return;

    pb[nb] = '\0';
    upperstr(pb);
    sequences[num].nb = nb;  /* # of total bases in sequence */
    sequences[num].seq = my_strdup(pb);
    free(pb);  /* use only as many/few memory space as needed */

    if (!raw_seq) {  /* break seq[] into only ACGT containing fragments */
        /* set any non-"ACGT" char to '\0', except for U which is converted to T */
        p0 = sequences[num].seq;  /* set pointer p0 to the sequence */
        convert_U_to_T(p0);  /* to allow for handling of RNA */

        while (*p0) {
            if (strchr(CASE_BASES, *p0) == NULL)  /* case-insensitive */
                *p0 = '\0';
            p0++;
        }
    }

    /* count number of fragments */
    p0 = sequences[num].seq;  /* reset pointer p0 */
    i = 0;
    while (i < nb)
        if (p0[i] != '\0') {  /* fragment breaking point */
            num_fragments++;
            k = strlen(p0 + i);  /* current fragment length */
            i += k;
        } else
            i++;

    /* now filling out the arrays */
    if (!num_fragments)  /* ID w/o corresponding sequence */
        sequences[num].num_nb = NULL;
    else
        sequences[num].num_nb = lvector(1, num_fragments);  /* 1-index */

    seq_fragments = allocate_seq_fragments(num_fragments);

    num_fragments = 0;
    i = 0;
    while (i < nb)
        if (p0[i] != '\0') {  /* "ACGT" */
            num_fragments++;
            k = strlen(p0 + i);  /* current fragment length */
            sequences[num].num_nb[num_fragments] = k;
            seq_fragments[num_fragments] = p0 + i;  /* pointer to seq array */
            i += k;
        } else
            i++;

    sequences[num].num_fragments = num_fragments;
    sequences[num].seq_fragments = seq_fragments;
}

void write_oneline_sequences(char *midx, struct_seqArr * seqs, char *filename)
{
    long i;
    FILE *fp;
    struct_seq *seq;

    fp = open_file(filename, "w");

    for (i = 1; i <= seqs->num_seqs; i++)
        if (!midx || (midx && midx[i])) {
            seq = &seqs->sequences[i];
            fprintf(fp, ">%s\t%ld\n", seq->id, seq->nb);
            fprintf(fp, "%s\n", seq->seq);
        }

    close_file(fp);
}

void write_hjb_sequences(char *midx, struct_seqArr * seqs, char *filename)
{
    long i;
    FILE *fp;
    struct_seq *seq;

    fp = open_file(filename, "w");

    for (i = 1; i <= seqs->num_seqs; i++)
        if (!midx || (midx && midx[i])) {
            seq = &seqs->sequences[i];
            fprintf(fp, "%s:%ld:-%ld:%s:\n", seq->id, seq->nb, seq->nb, seq->seq);
        }

    close_file(fp);
}

void format_fasta_sequence(FILE * fp, long nb, char *seq)
{
    char *str;
    long k, j, nline;

    str = cvector(0, BASE_PER_LINE);

    nline = nb / BASE_PER_LINE;
    if (nb % BASE_PER_LINE)
        nline++;

    k = 0;  /* offset */
    for (j = 1; j <= nline; j++) {
        strncpy(str, seq + k, BASE_PER_LINE);
        k += BASE_PER_LINE;
        str[BASE_PER_LINE] = '\0';
        fprintf(fp, "%s\n", str);
    }

    free_cvector(str, 0, DUMMY);
}

void write_one_seq(FILE * fp, struct_seq * seq)
{
    long nseq, nb_sum;

    nb_sum = 0;  /* total number of 'ACGT' bases */
    for (nseq = 1; nseq <= seq->num_fragments; nseq++)
        nb_sum += seq->num_nb[nseq];
    fprintf(fp, ">%s %ld %ld %ld\n", seq->id, seq->nb, nb_sum, seq->num_fragments);

    for (nseq = 1; nseq <= seq->num_fragments; nseq++) {
        if (nseq > 1) {  /* as a clear separator */
            fprintf(fp, "x #");  /* x as a separator */
            print_sep(fp, DASH, BASE_PER_LINE - 3);
        }

        format_fasta_sequence(fp, seq->num_nb[nseq], seq->seq_fragments[nseq]);
    }
}

void write_sequences(char *midx, struct_seqArr * seqs, char *filename)
{
    long i;
    FILE *fp;

    fp = open_file(filename, "w");

    for (i = 1; i <= seqs->num_seqs; i++)
        if (!midx || (midx && midx[i]))
            write_one_seq(fp, &seqs->sequences[i]);

    close_file(fp);
}

void write_sequences_in_tdat_order(char *filename, struct_seqArr * seqs,
                                   struct_data * tdat)
{
    char *str;
    long i, idx;
    FILE *fp;

    fp = open_file(filename, "w");

    str = cvector(0, BASE_PER_LINE);
    for (i = 1; i <= tdat->nrow; i++) {
        idx = tdat->seqidx[i];
        if (idx)
            write_one_seq(fp, &seqs->sequences[idx]);
    }
    free_cvector(str, 0, DUMMY);

    close_file(fp);
}

void convert_sequences_acgt_to_0123(struct_seqArr * seqs)
{
    long i, j;
    struct_seq *seq;

    for (i = 1; i <= seqs->num_seqs; i++) {
        seq = &seqs->sequences[i];
        for (j = 1; j <= seq->num_fragments; j++)
            seq2number(seq->seq_fragments[j]);
    }
}

void reverse_cmpl_sequences(struct_seqArr * seqs)
{
    long i, j;
    struct_seq *seq;

    for (i = 1; i <= seqs->num_seqs; i++) {
        seq = &seqs->sequences[i];
        for (j = 1; j <= seq->num_fragments; j++)
            reverse_cmpl(seq->seq_fragments[j], TRUE);
    }
}

void reverse_sequences(struct_seqArr * seqs)
{
    long i, j;
    struct_seq *seq;

    for (i = 1; i <= seqs->num_seqs; i++) {
        seq = &seqs->sequences[i];
        for (j = 1; j <= seq->num_fragments; j++)
            reverse_string(seq->seq_fragments[j]);  /* in-place reverse */
    }
}

long get_number_of_fasta_entries(char *seqfile)
/* count the number of '>' starting at the beginning of each line */
{
    char *p0, *line;
    long num = 0;
    FILE *fp;

    fp = open_file(seqfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = ltrim(p0);  /* allowing for leading white spaces before '>' */
        if (line[0] == '>')
            num++;
        free(p0);
    }

    close_file(fp);

    return num;
}

void free_seq_entry(struct_seq * seq)
{
    free_cvector(seq->title, 0, DUMMY);
    free_cvector(seq->id, 0, DUMMY);
    free_lvector(seq->num_nb, 1, DUMMY);
    free(seq->seq_fragments);  /* sequence fragments */
    free_cvector(seq->seq, 0, DUMMY);  /* whole sequence */
}

void free_sequences(char *fidx, struct_seqArr * seqs)
{
    long i;

    for (i = 1; i <= seqs->num_seqs; i++)
        if (!fidx || (fidx && fidx[i]))
            free_seq_entry(&seqs->sequences[i]);

    free(seqs->sequences);
}

char **allocate_seq_fragments(long num)
{
    char **str;

    str = (char **) malloc((num + 1) * sizeof(char *));

    if (str == NULL)
        fatal("malloc failure for allocate_seq_fragments()\n");

    return str;
}

struct_seq *allocate_memory_for_sequences(long num)
{
    struct_seq *sequences;

    sequences = (struct_seq *) malloc((num + 1) * sizeof(struct_seq));

    if (sequences == NULL)
        fatal("malloc failure forallocate_memory_for_sequences()\n");

    return sequences;
}

long base_offset(int base)
{
    long offset = -1;

    if (base == 'A')
        offset = 0;
    else if (base == 'C')
        offset = 1;
    else if (base == 'G')
        offset = 2;
    else if (base == 'T')
        offset = 3;
    else
        fatal("illegal character '%c' [non-%s] in base_offset()\n", base, BASES);

    return offset;
}

/* convert sequence "ACGT" to "0123" */
void seq2number(char *p)
{
    char c;

    while (*p) {
        c = (char) toupper((int) *p);

        if (c == 'A')
            c = '\0';
        else if (c == 'C')
            c = '\1';
        else if (c == 'G')
            c = '\2';
        else if (c == 'T')
            c = '\3';
        else  /* this should not happen! */
            fatal("illegal character '%c' [non-%s] in seq2number()\n", c, BASES);

        *p++ = c;
    }
}

/* convert "0123" to sequence "ACGT" */
void number2seq(char *p, long len)
{
    long i;

    for (i = 0; i < len; i++)
        if (p[i] == 0)
            p[i] = 'A';
        else if (p[i] == 1)
            p[i] = 'C';
        else if (p[i] == 2)
            p[i] = 'G';
        else if (p[i] == 3)
            p[i] = 'T';
        else
            fatal("illegal number '%c' [non-0123] in number2seq()\n", p[i]);
}

char cvt_base_U2T(char base)
{
    if (base == RNA_U)
        base = DNA_T;
    else if (base == RNA_u)
        base = DNA_t;

    return base;
}

char cvt_base_T2U(char base)
{
    if (base == DNA_T)
        base = RNA_U;
    else if (base == DNA_t)
        base = RNA_u;

    return base;
}

void convert_U_to_T(char *seq)
/* convert each U to T in seq[] for handling of RNA */
{
    while (*seq) {
        *seq = cvt_base_U2T(*seq);
        seq++;
    }
}

void convert_T_to_U(char *seq)
/* convert each T to U in seq[] for LogoGenerator */
{
    while (*seq) {
        *seq = cvt_base_T2U(*seq);
        seq++;
    }
}

int seqid_compare(const void *v1, const void *v2)
/* comparison function based on 'id' from FASTA sequence */
{
    const struct_seq *p1, *p2;

    p1 = (const struct_seq *) v1;
    p2 = (const struct_seq *) v2;

    return case_strcmp(p1->id, p2->id);  /* case-insensitive */
}

void check_redundancy(char *seqfile, long num, struct_seq * sequences)
/* after sorting sequence by ID (case-insensitive), mark redundant
 * sequence IDs with REPEATED_ID: keep only the last */
{
    char msg[BUF512];
    long i, num_duplicate = 0;

    for (i = 1; i < num; i++)
        if (is_equal_case_string(sequences[i].id, sequences[i + 1].id)) {
            sprintf(msg, "\t%s", sequences[i].id);
            log_msg(msg);
            free_cvector(sequences[i].id, 0, DUMMY);
            sequences[i].id = my_strdup(REPEATED_ID);  /* mark out previous one */
        }

    for (i = 1; i <= num; i++)
        if (is_equal_string(sequences[i].id, REPEATED_ID))
            num_duplicate++;

    if (num_duplicate) {
        sprintf(msg, "*** %ld repeated IDs in sequence file <%s> ***",
                num_duplicate, seqfile);
        log_msg(msg);
        log_msg("*** only SEQUENCE for the last ID is used! ***");
    }
}

void match_ids_seqs(long num, struct_tag * str_tags, struct_seqArr * seqs, char *midx)
{
    char msg[BUF512];
    long i, idx, num_bad = 0;

    struct_seq id_key, *np;

    for (i = 1; i <= num; i++) {
        id_key.id = str_tags[i].str;  /* matched by ID; just a copy of pointer */
        np = (struct_seq *) bsearch(&id_key, seqs->sequences + 1, seqs->num_seqs, sizeof(struct_seq), seqid_compare);  /* case-insensitive */

        if (np == NULL) {  /* w/o corresponding sequence */
            num_bad++;
            sprintf(msg, "\t\t%ld ID [%s (%s)] has NO matched sequence", num_bad,
                    str_tags[i].str, str_tags[i].tag);
            log_msg(msg);

        } else {
            idx = np - seqs->sequences;
            midx[idx] = '\1';
        }
    }

    if (num_bad) {
        sprintf(msg, "\ttotal number of missing IDs: %ld", num_bad);
        log_msg(msg);
    }
}

long get_tdat_ID_match_number(struct_data * tdat)
{
    long i, num = 0;

    for (i = 1; i <= tdat->nrow; i++)
        if (tdat->seqidx[i])
            num++;

    return num;
}

void get_matched_seq_tdat(struct_seqArr * seqs, struct_data * tdat)
{
    char msg[BUF512], *fidx, **cln_ids;
    long i, idx, j, num;
    double **cln_data;

    struct_seq *cln_sequences;

    num = get_tdat_ID_match_number(tdat);
    sprintf(msg, "   >>>> number of matched IDs: %ld <<<<", num);
    log_msg(msg);

    cln_sequences = allocate_memory_for_sequences(num);
    cln_ids = allocate_char_marray(num);
    cln_data = dmatrix(1, tdat->ncol, 1, num);

    fidx = cvector(1, seqs->num_seqs);  /* free-idx */
    init_cvector(fidx, 1, seqs->num_seqs, '\1');  /* default to free all */

    num = 0;
    for (i = 1; i <= tdat->nrow; i++) {
        idx = tdat->seqidx[i];
        if (idx) {
            num++;
            cln_sequences[num] = seqs->sequences[idx];
            fidx[idx] = '\0';  /* keep the matched ones -- 0 means not to free */

            cln_ids[num] = my_strdup(tdat->ids[i]);
            for (j = 1; j <= tdat->ncol; j++)
                cln_data[j][num] = tdat->data[j][i];
        }
    }

    free_sequences(fidx, seqs);
    free_cvector(fidx, 1, DUMMY);

    seqs->num_seqs = num;
    seqs->sequences = cln_sequences;

    free_char_marray(1, tdat->nrow, tdat->ids);
    free_dmatrix(tdat->data, 1, DUMMY, 1, DUMMY);

    tdat->ids = cln_ids;
    tdat->nrow = num;
    tdat->data = cln_data;

    free_lvector(tdat->seqidx, 1, DUMMY);

    tdat->seqidx = lvector(1, num);
    set_parallel_idx(1, num, tdat->seqidx);
}

void get_seq_tdat(char *outdir, char *seqfile, struct_seqArr * seqs, char *measfile,
                  struct_data * tdat)
{
    char msg[BUF512], filename[BUF512];
    long raw_seq = Gvars.RAWSEQ, sort_seq = TRUE;

    read_sequences(seqfile, seqs, raw_seq, sort_seq);
    sprintf(msg, "\tnumber of initial sequences: %ld", seqs->num_seqs);
    log_msg(msg);

    read_tdt_data(measfile, tdat);
    sprintf(msg, "\tnumber of initial measurement data: %ld", tdat->nrow);
    log_msg(msg);

    link_expr2seq(tdat, seqs, outdir, measfile);

    get_matched_seq_tdat(seqs, tdat);

    tdat->resid = dmatrix(1, tdat->ncol, 1, tdat->nrow);
    copy_dmatrix(tdat->data, 1, tdat->ncol, 1, tdat->nrow, tdat->resid);

    tdat->oknum = lvector(1, tdat->ncol);
    get_tdat_oknum(tdat);

    sprintf(filename, "%s/clean_seq.fasta", outdir);
    write_sequences_in_tdat_order(filename, seqs, tdat);

    sprintf(filename, "%s/clean_data.tsv", outdir);
    write_tdat(tdat, filename, tdat->data, NULL);

    if (!raw_seq)
        convert_sequences_acgt_to_0123(seqs);
}

void change_seq_case(struct_seqArr * seqs, long base_case)
{
    long i, nseq;
    struct_seq *seq;

    if (base_case == UNSET_LVAL)  /* keep as is */
        return;

    for (i = 1; i <= seqs->num_seqs; i++) {
        seq = &seqs->sequences[i];
        for (nseq = 1; nseq <= seq->num_fragments; nseq++) {
            if (base_case) {
                upperstr(seq->seq_fragments[nseq]);
                upperstr(seq->id);

            } else {
                lowerstr(seq->seq_fragments[nseq]);
                lowerstr(seq->id);
            }
        }
    }
}

void read_windows(char *winfile, struct_winArr * wins, long width)
/* struct_winArr file format -- 5 fields per line, as follows:
 * [1] seqid  -- sequence ID from which to extract a fragment
 * [2] winid  -- window ID to specify the extracted sequence fragment
 * [3] strand -- +1 for forward direction; -1 for reverse complementary
 * [4] ibeg   -- starting position
 * [5] iend   -- ending position
 * [6] width   -- fragment length: iend - ibeg + 1
 * Note: either [5] or [6] (with -width option), but not both */
{
    long num;
    struct_win *windows;

    num = get_line_number(winfile, TRUE);  /* could contain invalid entries */
    windows = allocate_memory_for_windows(num);

    populate_windows(winfile, windows, width);

    wins->num_wins = num;
    wins->windows = windows;
}

void initialize_struct_win(struct_win * w)
{
    w->seqid = NULL;
    w->winid = NULL;
    w->strand = 0;
    w->ibeg = -1;
    w->iend = -1;
    w->width = -1;
}

void populate_windows(char *winfile, struct_win * windows, long width)
{
    char msg[BUF512], *p0, *line, *items[BUF512], str0[BUF512];
    long nitem, NUM_REQUIRED_ITEM = 5, num = 0, snum = 0;
    long strand, ib, ie;
    FILE *fp;

    struct_win *w;

    fp = open_file(winfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        snum++;
        line = trim(p0);

        if (is_skip_line(line)) {
            free(p0);
            continue;
        }

        num++;
        strcpy(str0, line);
        nitem = item_list(line, items, NUM_REQUIRED_ITEM, WSPACES);

        w = &windows[num];

        if (nitem == NUM_REQUIRED_ITEM) {
            strand = cvt2long(items[3]);
            ib = cvt2long(items[4]);
            ie = cvt2long(items[5]);

            if ((strand != 1 && strand != -1) || (!width && ib > ie) ||
                (ib == LONG_MAX || ib <= 0) || (ie == LONG_MAX || ie <= 0)) {
                sprintf(msg, "\tinvalid entry <line #%ld: %s>", snum, str0);
                log_msg(msg);
                initialize_struct_win(w);

            } else {
                w->seqid = my_strdup(items[1]);
                w->winid = my_strdup(items[2]);
                w->strand = strand;
                w->ibeg = ib;
                if (!width) {
                    w->iend = ie;
                    w->width = ie - ib + 1;
                } else {
                    w->width = ie;
                    w->iend = ib + ie - 1;
                }
            }

        } else {
            sprintf(msg, "\tinvalid entry <line #%ld: %s>", snum, str0);
            log_msg(msg);
        }

        free(p0);
    }

    close_file(fp);
}

void write_windows(char *outdir, char *filename, struct_seqArr * seqs,
                   struct_winArr * wins)
{
    char msg[BUF512], str[BUF512];
    long i, nok = 0;
    FILE *fp;

    struct_win *w;
    struct_seq id_key, *np;

    sprintf(str, "%s/%s", outdir, filename);
    fp = open_file(str, "w");

    for (i = 1; i <= wins->num_wins; i++) {
        w = &wins->windows[i];
        if (!w->seqid)  /* wins->window could contain invalid entries */
            continue;

        id_key.id = w->seqid;  /* matched by seq ID; just a copy of pointer */
        np = (struct_seq *) bsearch(&id_key, seqs->sequences + 1, seqs->num_seqs, sizeof(struct_seq), seqid_compare);  /* case-insensitive */

        if (np == NULL) {
            sprintf(msg, "\tno match for entry [%s %s]", w->seqid, w->winid);
            log_msg(msg);

        } else if (w->iend > np->nb) {
            sprintf(msg, "\tentry [%s %s] end position too large [%ld > %ld]",
                    w->seqid, w->winid, w->iend, np->nb);
            log_msg(msg);

        } else {
            nok++;
            fprintf(fp, ">%s [%s %+ld %ld %ld %ld]\n", w->winid, w->seqid, w->strand,
                    w->ibeg, w->iend, w->width);
            print_struct_win(fp, w, np);
        }
    }

    close_file(fp);

    sprintf(msg, "a total of %ld %s written", nok, (nok == 1) ? "entry" : "entries");
    log_msg(msg);
}

void print_struct_win(FILE * fp, struct_win * w, struct_seq * np)
{
    char *p;

    p = cvector(0, w->width);

    strncpy(p, np->seq + w->ibeg - 1, w->width);
    p[w->width] = '\0';

    if (w->strand == -1)
        reverse_cmpl(p, FALSE);

    format_fasta_sequence(fp, w->width, p);

    free_cvector(p, 0, DUMMY);
}

/* -------------------- processing logo fasta sequence -------------------- */
void free_logo_sequences(long num_seqs, struct_seq * sequences)
/* free the logo sequence data structure */
{
    long i;

    for (i = 1; i <= num_seqs; i++) {
        free_cvector(sequences[i].id, 0, DUMMY);
        free_cvector(sequences[i].seq, 0, DUMMY);
    }

    free(sequences);
}

void assign_logo_sequence(long num, char *pb, long nb, struct_seq * sequences)
/* c.f. assign_sequence() */
{
    if (!num)  /* no header line read in yet */
        return;

    pb[nb] = '\0';

    sequences[num].nb = nb;  /* # of total bases in sequence */
    sequences[num].seq = my_strdup(pb);

    free(pb);  /* use only as many/few memory space as needed */
}

struct_seq *read_logo_fasta(char *seqfile, long *num_seqs, char *outdir)
/* c.f. read_sequences() */
{
    long nb = 0, num, maxline;
    char *line, *p0, *p, *pb = NULL;
    FILE *fp;

    struct_seq *sequences;  /* 1-indexing */

    num = get_number_of_fasta_entries(seqfile);
    sequences = allocate_memory_for_sequences(num);

    log_msg_vchk("Reading in the logo sequence in fasta format");

    fp = open_file(seqfile, "r");
    num = 0;  /* re-initialization */
    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (is_fasta_skip_line(line)) {
            free(p0);
            continue;  /* skip empty and commented lines [; #] */
        }

        if (line[0] == '>') {  /* header line */
            assign_logo_sequence(num, pb, nb, sequences);  /* for previous sequence */
            p = strtok(line, WSPACES);  /* get gene ID */
            num++;
            sequences[num].id = my_strdup(p + 1);  /* without leading '>' */

            maxline = BUF1024;  /* reset starting size */
            pb = (char *) malloc(maxline * sizeof(char));  /* initial size */
            if (pb == NULL)
                fatal("malloc failure in reading sequence\n");
            nb = 0;

        } else if (num) {
            while (*line && !isspace((int) *line)) {  /* any char upto the 1st \w */
                if (nb >= maxline - 1)
                    pb = enlarge_cline(&maxline, pb);
                pb[nb++] = *line++;
            }
        }

        free(p0);  /* release memory */
    }

    assign_logo_sequence(num, pb, nb, sequences);  /* for the last sequence */
    close_file(fp);

    if (!num)
        fatal("no valid sequences in file [%s]\n", seqfile);

    qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);
    check_redundancy(seqfile, num, sequences);
    qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);

    check_logo_sequences(num, sequences, seqfile, outdir);

    *num_seqs = num;

    return sequences;
}

void check_logo_sequences(long num_seqs, struct_seq * sequences, char *seqfile,
                          char *outdir)
/* check for consistency of sequence length, change no-ACGT base to DASH */
{
    char msg[BUF512], str[BUF512], *p;
    long i, j, k, nb, num_err = 0, nb1 = sequences[1].nb;
    FILE *fp;

    sprintf(str, "%s/logoflat_%s", outdir, basename(seqfile));
    fp = open_file(str, "w");

    fprintf(fp, "# ID\tsequence\t# length\tserial-#\tdelta_len1\n# ");
    print_sep(fp, DASH, 66);

    for (i = 1; i <= num_seqs; i++) {
        nb = sequences[i].nb;
        p = sequences[i].seq;
        convert_U_to_T(p);  /* for handling of RNA */

        if (nb != nb1) {
            num_err++;
            sprintf(msg, "%ld: %s [%ld vs %ld]", i, p, nb, nb1);
            log_msg(msg);
        }

        for (j = 0; j < nb; j++)
            if ((strchr(CASE_BASES, p[j]) == NULL))
                p[j] = DASH;  /* convert non-ACGT to DASH */

        k = nb - nb1;
        fprintf(fp, "%s\t%s\t# %ld\t%ld\t%ld%s\n", sequences[i].id, sequences[i].seq,
                nb, i, k, (k ? "*" : ""));
    }

    close_file(fp);

    if (num_err) {
        sprintf(msg, "sequences in file <%s> are of unequal length [%ld]\n",
                seqfile, num_err);
        log_msg(msg);

        sprintf(msg, "check file <%s> for non-zero entries at the last column\n", str);
        log_msg_exit(msg);
    }
}

/* -------------------- processing flat sequence -------------------- */
struct_seq *read_logo_flat(char *seqfile, long *num_seqs, char *outdir)
{
    long num, nitem, NUM_REQUIRED_ITEM = 2;
    char *items[BUF512], *line, *p0;
    FILE *fp;

    struct_seq *sequences;  /* 1-indexing */

    num = get_line_number(seqfile, TRUE);
    sequences = allocate_memory_for_sequences(num);

    log_msg_vchk("Reading in the logo sequence in flat format\n");

    fp = open_file(seqfile, "r");
    num = 0;  /* re-initialization */
    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (is_skip_line(line)) {
            free(p0);
            continue;
        }

        num++;
        nitem = item_list(line, items, NUM_REQUIRED_ITEM, WSPACES);
        if (nitem != NUM_REQUIRED_ITEM)
            fatal("wrong format <%s>: require 'ID\tsequence'\n", line);

        sequences[num].id = my_strdup(items[1]);
        sequences[num].seq = my_strdup(items[2]);
        sequences[num].nb = strlen(sequences[num].seq);

        free(p0);  /* release memory */
    }
    close_file(fp);

    if (!num)
        fatal("no valid sequences in file [%s]\n", seqfile);

    qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);
    check_redundancy(seqfile, num, sequences);
    qsort(sequences + 1, num, sizeof(sequences[1]), seqid_compare);

    check_logo_sequences(num, sequences, seqfile, outdir);

    *num_seqs = num;

    return sequences;
}
