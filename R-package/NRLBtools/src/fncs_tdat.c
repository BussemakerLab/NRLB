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

void read_tdt_data(char *tdtfile, struct_data * tdat)
/* read in tab-delimited-text data file */
{
    char msg[BUF512];
    long ncol, nrow, header_num, header_type;

    sprintf(msg, "reading in tab-delimited data file [%s]", tdtfile);
    log_msg(msg);

    get_tdat_col_info(tdtfile, &header_num, &header_type, &ncol);

    tdat->col_names = allocate_char_marray(ncol);
    get_tdat_col_names(tdtfile, header_num, header_type, ncol, tdat->col_names);

    get_tdat_nrow(tdtfile, header_num, ncol, &nrow);

    tdat->ids = allocate_char_marray(nrow);
    tdat->data = dmatrix(1, ncol, 1, nrow);
    read_tdat(tdtfile, header_num, tdat->ids, tdat->data);

    tdat->seqidx = lvector(1, nrow);  /* linkage to sequence */

    tdat->nrow = nrow;
    tdat->ncol = ncol;

    tdat->oknum = NULL;
    tdat->resid = NULL;
}

void free_tdt_data(struct_data * tdat)
{
    free_char_marray(1, tdat->ncol, tdat->col_names);
    free_char_marray(1, tdat->nrow, tdat->ids);
    free_dmatrix(tdat->data, 1, DUMMY, 1, DUMMY);
    free_lvector(tdat->seqidx, 1, DUMMY);

    if (tdat->oknum)
        free_lvector(tdat->oknum, 1, DUMMY);

    if (tdat->resid)
        free_dmatrix(tdat->resid, 1, DUMMY, 1, DUMMY);
}

void get_tdat_col_info(char *tdtfile, long *header_num, long *header_type, long *ncol)
/* the 1st valid line is taken as the header line if
 *    a) it has ONE field LESS than the following lines
 * or b) second and following fields are ALL non-numeric */
{
    char *p0, *p = NULL, *line, *items[BUFBIG + 1], *items2[BUFBIG + 1];
    long i, nitem = 0, nitem2 = 0, num = 0;
    FILE *fp;

    fp = open_file(tdtfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        num++;  /* total number of lines */
        line = ltrim(p0);  /* NOT trim to protect possible ending \t */

        if (!is_skip_line(line)) {  /* 1st non-empty or comment line */
            nitem = csplit(line, items, BUFBIG, TABCHAR);
            p = p0;  /* keep a copy of the 1st valid line in p */
            break;  /* can't free(p0) yet, because of items[] refer to it */
        }

        free(p0);
    }

    *header_num = num;  /* including all lines upto the 1st valid line */

    /* presumably we have got the header line (in items[]); now read in a
     * data line to check for the number of fields */
    while ((p0 = my_getline(fp)) != NULL) {
        line = ltrim(p0);

        if (!is_skip_line(line)) {  /* 2nd non-empty or comment line */
            nitem2 = csplit(line, items2, BUFBIG, TABCHAR);
            if (nitem2 < 2)
                fatal("<%s> must be tab-delimited with >= ONE data field column:"
                      " (%s)!\n", tdtfile, p0);
            free(p0);  /* items2[] no longer needed -- only nitem2 used later */
            break;
        }

        free(p0);
    }

    if (nitem2 == nitem + 1)
        *header_type = 1;  /* 1..nitem */

    else if (nitem2 == nitem) {
        num = 0;
        for (i = 2; i <= nitem; i++)
            num += is_numeric(items[i]);  /* pure numeric field */

        if (!num)  /* all fields are non-numeric: taken as headline */
            *header_type = 2;  /* 2..nitem */

        else {  /* at least one numeric field: taken as a data line */
            *header_type = -1;  /* no header line */
            (*header_num)--;
        }

    } else
        fatal("number of fields in the first TWO lines does NOT match: "
              "<%ld vs. %ld>\n", nitem, nitem2);

    free(p);  /* items[] no longer needed! */

    close_file(fp);

    *ncol = nitem2 - 1;
}

void get_tdat_col_names(char *tdtfile, long header_num, long header_type, long ncol,
                        char **col_names)
{
    char *p0, *line, *items[BUFBIG + 1], str[BUF512];
    long i, nitem, num = 0;
    FILE *fp;

    if (header_type == -1) {
        for (i = 1; i <= ncol; i++) {
            sprintf(str, "col_%3.3ld", i);
            col_names[i] = my_strdup(str);
        }
        return;
    }

    if (header_type != 1 && header_type != 2)
        fatal("illegal header type <%ld>: [-1, 1, 2]\n", header_type);

    fp = open_file(tdtfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        num++;
        if (num == header_num)
            break;
        free(p0);
    }

    line = ltrim(p0);  /* NOT trim to protect possible ending \t */
    nitem = csplit(line, items, BUFBIG, TABCHAR);

    if (header_type == 1) {
        for (i = 1; i <= nitem; i++) {
            if (Gvars.RAWID)
                get_alnum(items[i], '_');  /* convert non-alphanumeric chars to '_' */
            col_names[i] = my_strdup(items[i]);
        }

    } else {  /* header_type == 2 */
        for (i = 2; i <= nitem; i++) {
            if (Gvars.RAWID)
                get_alnum(items[i], '_');
            col_names[i - 1] = my_strdup(items[i]);
        }
    }

    free(p0);

    close_file(fp);
}

void get_tdat_nrow(char *tdtfile, long header_num, long ncol, long *nrow)
{
    char *p0, *line, *items[BUFBIG + 1];
    long line_number, nitem, nitem2 = ncol + 1, num = 0;
    FILE *fp;

    fp = open_file(tdtfile, "r");

    skip_lines(fp, header_num);

    line_number = header_num;  /* used for real line number info */
    num = 0;

    while ((p0 = my_getline(fp)) != NULL) {
        line_number++;
        line = ltrim(p0);

        if (!is_skip_line(line)) {
            nitem = csplit(line, items, BUFBIG, TABCHAR);
            if (nitem != nitem2)
                fatal("data line '#%ld [%s]' has unequal # of fields <%ld vs. %ld>\n",
                      line_number, p0, nitem, nitem2);
            num++;
        }

        free(p0);
    }

    close_file(fp);

    if (!num)
        fatal("NO valid data entries in tab-delimited file <%s>\n", tdtfile);

    *nrow = num;  /* total # of lines with equal # of data fields */
}

void read_tdat(char *tdtfile, long header_num, char **ids, double **data)
/* all fields must be tab-delimited; empty lines and lines starting with
 * '#' are ignored; NB: "strtok" is GREEDY -- it consums consecutive \t, so
 * use csplit() */
{
    char *p0, *line, *items[BUFBIG + 1];
    long i, nitem, num = 0;
    FILE *fp;

    fp = open_file(tdtfile, "r");

    skip_lines(fp, header_num);

    num = 0;
    while ((p0 = my_getline(fp)) != NULL) {
        line = ltrim(p0);

        if (is_skip_line(line)) {
            free(p0);
            continue;  /* skip empty and commented lines */
        }

        num++;
        nitem = csplit(line, items, BUFBIG, TABCHAR);
        ids[num] = my_strdup(items[1]);  /* first field for id */
        for (i = 2; i <= nitem; i++)  /* data fields */
            data[i - 1][num] = cvt2double(items[i]);

        free(p0);
    }

    close_file(fp);
}

static void write_tdat_header(FILE * fp, struct_data * tdat)
{
    double r;
    long i, k, ncol = tdat->ncol, nrow = tdat->nrow;

    fprintf(fp, "## %ld data column%s:\n", ncol, (ncol == 1) ? "" : "s");

    if (tdat->oknum == NULL) {
        tdat->oknum = lvector(1, tdat->ncol);
        get_tdat_oknum(tdat);
    }

    k = get_number_of_valid_rows(ncol, nrow, tdat->data);
    r = (double) k / nrow * 100.0;
    fprintf(fp, "## \tnumber of all-valid rows: %ld/%ld (%.1f)\n", k, nrow, r);

    for (i = 1; i <= ncol; i++) {
        r = (double) tdat->oknum[i] / nrow * 100;
        fprintf(fp, "## %3.3ld (%ld/%ld=%.1f) -- %s\n", i, tdat->oknum[i], nrow, r,
                tdat->col_names[i]);
    }

    for (i = 1; i <= ncol; i++)
        fprintf(fp, "\t%s", tdat->col_names[i]);
    fprintf(fp, "\n");
}

void write_dval(FILE * fp, double dval)
{
    if (dval < XBIG_CUTOFF)
        fprintf(fp, "\t%g", dval);
    else
        fprintf(fp, "\tNA");
}

void write_tdat_row(FILE * fp, long row_idx, struct_data * tdat, double **dmtx)
{
    long j;

    if (tdat->seqidx[row_idx]) {  /* with matching sequence */
        fprintf(fp, "%s", tdat->ids[row_idx]);

        for (j = 1; j <= tdat->ncol; j++)
            write_dval(fp, dmtx[j][row_idx]);

        fprintf(fp, "\n");
    }
}

void write_tdat_transpose(struct_data * tdat, char *filename, double **dmtx,
                          struct_id * my_ids)
{
    long i, idx, k, num_ok = 0;
    FILE *fp;

    fp = open_file(filename, "w");

    for (i = 1; i <= tdat->nrow; i++) {
        idx = my_ids ? my_ids[i].idx : i;
        if (tdat->seqidx[idx])
            num_ok++;
    }

    fprintf(fp, "### transpose: %ld column%s ###\n", num_ok, (num_ok == 1) ? "" : "s");

    for (i = 1; i <= tdat->nrow; i++) {
        idx = my_ids ? my_ids[i].idx : i;
        if (tdat->seqidx[idx])
            fprintf(fp, "\t%s", tdat->ids[idx]);
    }
    fprintf(fp, "\n");

    for (i = 1; i <= tdat->ncol; i++) {
        fprintf(fp, "%s", tdat->col_names[i]);

        for (k = 1; k <= tdat->nrow; k++) {
            idx = my_ids ? my_ids[k].idx : k;
            if (tdat->seqidx[idx])
                write_dval(fp, dmtx[i][idx]);
        }
        fprintf(fp, "\n");
    }

    close_file(fp);
}

void write_tdat(struct_data * tdat, char *filename, double **dmtx, struct_id * my_ids)
{
    long i, k;
    FILE *fp;

    fp = open_file(filename, "w");

    write_tdat_header(fp, tdat);

    for (i = 1; i <= tdat->nrow; i++) {
        k = my_ids ? my_ids[i].idx : i;
        write_tdat_row(fp, k, tdat, dmtx);
    }

    close_file(fp);
}

struct_id *get_sorted_ids(long num_ids, char **ids, char *filename)
{
    char msg[BUF512];
    long i, num_duplicate = 0;
    struct_id *my_ids;

    my_ids = allocate_memory_for_idstr(num_ids);

    for (i = 1; i <= num_ids; i++) {
        my_ids[i].idx = i;
        my_ids[i].idstr = my_strdup(ids[i]);
    }

    qsort(my_ids + 1, num_ids, sizeof(my_ids[1]), idstr_compare);

    for (i = 1; i < num_ids; i++)
        if (is_equal_string(my_ids[i].idstr, my_ids[i + 1].idstr)) {
            free_cvector(my_ids[i].idstr, 0, DUMMY);
            my_ids[i].idstr = my_strdup(REPEATED_ID);  /* mark out previous one */
        }

    for (i = 1; i <= num_ids; i++)
        if (is_equal_string(my_ids[i].idstr, REPEATED_ID))
            num_duplicate++;

    if (num_duplicate) {
        sprintf(msg, "*** %ld repeated IDs in file <%s> ***", num_duplicate, filename);
        log_msg(msg);
        log_msg("*** only the last entry for the ID is used ***");
    }

    qsort(my_ids + 1, num_ids, sizeof(my_ids[1]), idstr_compare);

    return my_ids;
}

void match_ids_tdat(long num, struct_tag * str_tags, long nrow, struct_id * my_ids,
                    long *midx)
{
    char msg[BUF512];
    long i, idx, num_ok = 0, num_bad = 0;

    struct_id id_key, *np;

    for (i = 1; i <= num; i++) {
        id_key.idstr = str_tags[i].str;  /* just a copy of pointer */
        np = (struct_id *) bsearch(&id_key, my_ids + 1, nrow, sizeof(struct_id),
                                   idstr_compare);

        if (np == NULL) {  /* w/o corresponding tdat entry */
            num_bad++;
            sprintf(msg, "\t\t%ld ID [%s (%s)] w/o matched tdat entry", num_bad,
                    str_tags[i].str, str_tags[i].tag);
            log_msg(msg);

        } else {
            num_ok++;
            idx = my_ids[np - my_ids].idx;  /* link to tdat row */
            midx[idx] = 1;
        }
    }

    if (num_bad) {
        sprintf(msg, "\ttotal number of missing IDs: %ld", num_bad);
        log_msg(msg);
    }

    if (!num_ok)
        fatal("no valid IDs in your list!\n");
}

void mask_off_na_rows(char *fitmtx_file, struct_data * tdat)
{
    char msg[BUF512];
    long i, j, num_na_rows = 0;

    for (i = 1; i <= tdat->nrow; i++) {
        for (j = 1; j <= tdat->ncol; j++)
            if (tdat->data[j][i] > XBIG_CUTOFF)
                break;

        if (j <= tdat->ncol) {
            num_na_rows++;
            free_cvector(tdat->ids[i], 0, DUMMY);
            tdat->ids[i] = my_strdup(REPEATED_ID);
        }
    }

    sprintf(msg, "%s: ignoring %ld NA-contaning row(s)", fitmtx_file, num_na_rows);
    log_msg(msg);
}

void write_problematic_ids(char *outdir, char *tdtfile, long nrow, char **ids,
                           long num_missing_seq, long num_duplicate, long *check_idx)
{
    char filename[BUF512];
    long i;
    FILE *fp;

    if (!num_missing_seq && !num_duplicate)  /* all are good */
        return;

    if (nrow == num_missing_seq)
        fatal("no matching sequence for data file <%s>: check ids!\n", tdtfile);

    sprintf(filename, "%s/check_%s", outdir, basename(tdtfile));

    fp = open_file(filename, "w");

    if (num_missing_seq) {
        fprintf(fp, "%s%ld\n", NUM_PROBE_NO_SEQ, num_missing_seq);
        for (i = 1; i <= nrow; i++)
            if (check_idx[i] == -1)
                fprintf(fp, "\t%s\n", ids[i]);
    }

    if (num_duplicate) {
        print_sep(fp, DASH, 76);
        fprintf(fp, "%s%ld\n", NUM_PROBE_REPEAT, num_duplicate);
        for (i = 1; i <= nrow; i++)
            if (check_idx[i] == 1)
                fprintf(fp, "\t%s\n", ids[i]);
    }

    close_file(fp);
}

void link_expr2seq(struct_data * tdat, struct_seqArr * seqs, char *outdir, char *tdtfile)
/* link measurement data with corresponding sequence by matching ids */
{
    long i, idx, k, nrow = tdat->nrow;
    long num_missing_seq = 0, num_duplicate = 0;
    long *check_idx, *seq2expr_idx;

    struct_seq id_key, *np;

    seq2expr_idx = lvector(1, seqs->num_seqs);
    check_idx = lvector(1, nrow);

    for (i = 1; i <= nrow; i++) {
        id_key.id = tdat->ids[i];  /* matched by ID; just a copy of pointer */

        np = (struct_seq *) bsearch(&id_key, seqs->sequences + 1, seqs->num_seqs, sizeof(struct_seq), seqid_compare);  /* case-insensitive */

        if (np == NULL) {  /* w/o corresponding sequence */
            check_idx[i] = -1;
            num_missing_seq++;

        } else {
            idx = np - seqs->sequences;
            k = seq2expr_idx[idx];  /* linking from sequence to measurement */
            if (k) {  /* repeated probe ID in measurement file: keep the last */
                tdat->seqidx[k] = 0;  /* mask the previous one */
                check_idx[k] = 1;
                num_duplicate++;
            }

            tdat->seqidx[i] = idx;  /* link from expr row to sequence info */
            seq2expr_idx[idx] = i;  /* link back from seq to measurement row */
        }
    }

    write_problematic_ids(outdir, tdtfile, nrow, tdat->ids, num_missing_seq,
                          num_duplicate, check_idx);

    free_lvector(seq2expr_idx, 1, DUMMY);
    free_lvector(check_idx, 1, DUMMY);
}

void change_id_case(struct_data * tdat, long base_case)
{
    long i;

    if (base_case == UNSET_LVAL)  /* keep as is */
        return;

    for (i = 1; i <= tdat->nrow; i++) {
        if (base_case)
            upperstr(tdat->ids[i]);
        else
            lowerstr(tdat->ids[i]);
    }
}

void log_transform_tdat(struct_data * tdat, char *log_base)
{
    long i, j;
    double dval, cvt;

    if (is_empty_string(log_base))
        return;

    cvt = (is_equal_string(log_base, "e")) ? 1.0 : log(cvt2double(log_base));

    for (i = 1; i <= tdat->nrow; i++) {
        for (j = 1; j <= tdat->ncol; j++) {
            dval = tdat->data[j][i];

            if (dval <= 0)
                fatal("can't take log of %g [id %s (%ld)]\n", dval, tdat->ids[i], i);

            if (dval < XBIG_CUTOFF)  /* valid number */
                tdat->data[j][i] = log(dval) / cvt;
        }
    }
}

void get_tdat_oknum(struct_data * tdat)
{
    long i, j;

    for (i = 1; i <= tdat->ncol; i++) {
        for (j = 1; j <= tdat->nrow; j++)
            if (tdat->data[i][j] < XBIG_CUTOFF)
                tdat->oknum[i]++;
    }
}

long get_number_of_valid_rows(long nvec, long nlen, double **Xmtx)
/* all entries in a row must be valid to be counted */
{
    long i, j, num = 0;

    for (i = 1; i <= nlen; i++) {
        for (j = 1; j <= nvec; j++)
            if (Xmtx[j][i] > XBIG_CUTOFF)
                break;

        if (j > nvec)
            num++;
    }

    return num;
}

void check_columns(struct_data * tdat, char *columns)
{
    char **ok_col_names;
    double **ok_data;
    long i, k, num, vnum[BUF512];  /* assuming maximum BUF512: enough command line */

    if (is_empty_string(columns))
        return;

    init_lvector(vnum, 0, BUF512 - 1, 0);
    num = extract_numlist(columns, LIST_SEP, 1, tdat->ncol, vnum);

    if (!num) {
        free_tdt_data(tdat);
        log_msg_exit("Your have picked 0 columns");
    }

    if (num == tdat->ncol) {  /* all is included! */
        log_msg("all columns included");
        return;
    }

    ok_col_names = allocate_char_marray(num);
    ok_data = dmatrix(1, num, 1, tdat->nrow);

    for (i = 1; i <= num; i++) {
        k = vnum[i];
        ok_col_names[i] = my_strdup(tdat->col_names[k]);
        copy_dvector(tdat->data[k], 1, tdat->nrow, ok_data[i]);
    }

    free_char_marray(1, tdat->ncol, tdat->col_names);
    free_dmatrix(tdat->data, 1, DUMMY, 1, DUMMY);

    /* now the selected data */
    tdat->ncol = num;
    tdat->col_names = ok_col_names;
    tdat->data = ok_data;
}
