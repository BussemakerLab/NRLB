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

void read_topos(char *topofile, long list, struct_topoArr * topo)
{
    char msg[BUF512];
    long num_topo;
    struct_tag *str_tags;
    struct_topo *topos;

    log_msg("reading topological patterns");

    str_tags = fillup_str_tags(list, topofile, &num_topo);
    expand_abbr_topos(&num_topo, str_tags);
    get_unique_strtag(&num_topo, str_tags);

    sprintf(msg, "\tnumber of topological patterns: %ld", num_topo);
    fprintf(Gvars.PRGLOG, "%s\n", msg);

    topos = allocate_memory_for_topos(num_topo);
    populate_topos(num_topo, topos, str_tags);

    topo->num_topo = num_topo;
    topo->topos = topos;

    free_strtags(1, num_topo, str_tags);
}

void unify_parens_chars(char *str)
{
    cvtstr_set1toc2(str, "[{<", '(');
    cvtstr_set1toc2(str, "]}>", ')');
}

/* The topological patterns can be supplied in short-hand form, e.g.,
 * X4-4X3 for XXXX----XXX; this function performs the transformation
 * and several quality control checkings. Note that topo-patterns are
 * just short hand forms which will be enumerated at each position to
 * A|C|G|T(U). X is for an unknown position [could be A|C|G|T(U)] to
 * be enumerated, - is for gap position to be skipped. N is converted
 * to - internally. Additionally, () pairs stand for matched canonical
 * base pairs, used for specifying double-helix stem region. */
void expand_abbr_topos(long *num, struct_tag * str_tags)
{
    char expanded_str[BUF512], msg[BUF512];
    char str[BUF512], tag[BUF512], *accept = "X-()0123456789";
    long i, Xcount, num_ok = 0;

    log_prg("    expanding & checking...");

    for (i = 1; i <= *num; i++) {
        strcpy(str, str_tags[i].str);
        strcpy(tag, str_tags[i].tag);

        /* cf. validate_expanded_pattern() -- with 'X' here */
        upperstr(str);  /* for consistence, all converted to upper case */
        unify_parens_chars(str);
        cvtstr_c1toc2(str, 'N', '-');  /* 'N' to '-' */

        expand_char_num_to_full(str, expanded_str);

        if (!string_contains_only_those_characters(str, accept)) {
            sprintf(msg, "\t***IGNORE pattern [%s %ld: %s] -- containing invalid"
                    " character(s) [%s]", tag, i, str, accept);
            log_prg(msg);
            continue;
        }

        if (isdigit((int) str[0])) {
            sprintf(msg, "\t***IGNORE pattern [%s %ld: %s] -- must start"
                    " with '-' or 'X|x', NOT a digit", tag, i, str);
            log_prg(msg);
            continue;
        }

        if (!has_matched_parens(expanded_str)) {
            sprintf(msg, "\t***IGNORE pattern [%s %ld: %s] -- has unmatched"
                    " () pairs", tag, i, str);
            log_prg(msg);
            continue;
        }

        Xcount = char_count_in_string(expanded_str, XCHR);

        if (Xcount > MAXNTS_TOP) {
            sprintf(msg, "\t***IGNORE pattern [%s %ld: %s(%s)] -- X characters "
                    "%ld > %ld", tag, i, str, expanded_str, Xcount, MAXNTS_TOP);
            log_prg(msg);

        } else {  /* tidy up version */
            num_ok++;
            free_p_strtag(&str_tags[num_ok]);
            init_strtag_with_strings(expanded_str, tag, &str_tags[num_ok]);
        }
    }

    free_extra_strtags(num_ok + 1, *num, str_tags);

    if (num_ok == 0)
        fatal("\tno valid topological patterns\n");

    *num = num_ok;
}

/* e.g., c4s2tt ---> ccccsstt */
long expand_char_num_to_full(char *src, char *dst)
{
    char c, str[BUF512], msg[BUF512];
    long ir, nr, idx = -1, j = 0, m = 0;
    long len = strlen(src);

    if (len == 0)
        return FALSE;

    if (isdigit((int) src[0])) {
        sprintf(msg, "\t<%s>: starts with digit", src);
        log_prg(msg);
        return FALSE;
    }

    while (j < len) {
        c = src[j];

        if (!isdigit((int) c)) {
            dst[m++] = c;
            j++;

        } else {  /* extract the digits */
            idx = 0;
            while (isdigit((int) src[j])) {
                str[idx++] = src[j];
                j++;
            }
            str[idx] = '\0';

            nr = cvt2long(str) - 1;  /* previous char already counted */
            c = dst[m - 1];  /* char to be repeated */
            for (ir = 0; ir < nr; ir++)
                dst[m + ir] = c;
            m += nr;
        }
    }

    dst[m] = '\0';

    if (idx >= 0) {
        sprintf(msg, "\texpanding %s ===> %s", src, dst);
        log_prg(msg);
    }

    return TRUE;  /* all fine */
}

long char_count_in_string(char *str, char c)
{
    long count = 0;

    while (*str)
        if (*str++ == c)
            count++;

    return count;
}

void populate_topos(long num_topo, struct_topo * topos, struct_tag * str_tags)
{
    char msg[BUF512], *p;  /* for simplicity */
    long i, j, Xcount;
    size_t nlen;
    struct_topo *t;  /* as a short hand */

    sprintf(msg, "summary of %ld topology pattern%s:", num_topo,
            (num_topo > 1) ? "s" : "");
    log_prg(msg);

    for (i = 1; i <= num_topo; i++) {
        t = &topos[i];
        p = str_tags[i].str;
        nlen = strlen(p);

        init_strtag(&str_tags[i], &t->str_tag);
        t->Tcount = nlen;

        t->offset = lvector(0, nlen);  /* 0-indexed: Tcount */

        Xcount = 0;
        for (j = 0; j < (long) nlen; j++)
            if (p[j] == XCHR) {
                t->offset[Xcount] = j;
                Xcount++;
            }

        t->Xcount = Xcount;
        t->msize = 1 << 2 * Xcount;  /* 4^Xcount; ONE element more */

        t->num_parens = fillup_parens(p, BUF16, t->parens_bidx, t->parens_eidx);

        sprintf(msg, "%s %4ld: Tcount=%2ld; Xcount=%ld; ()=%ld %s [%s]", t->str_tag.tag,
                i, t->Tcount, t->Xcount, t->num_parens, t->str_tag.str, t->str_tag.tag);
        log_prg(msg);

        sprintf(msg, "\t%2ld %s\t%s", i, t->str_tag.tag, t->str_tag.str);
        log_run(msg);
    }
}

void write_topos(struct_topoArr * topo, char *outdir, char *topo_file)
{
    char str[BUF512];
    long i;
    FILE *fp;

    struct_topo *t;

    sprintf(str, "%s/%s", outdir, topo_file);
    fp = open_file(str, "w");

    fprintf(fp, "# cleaned up and expanded topological patterns\n");
    for (i = 1; i <= topo->num_topo; i++) {
        t = &topo->topos[i];
        fprintf(fp, "%-30s\t%-20s\t#%ld\t%ld\t%ld\n", t->str_tag.str, t->str_tag.tag,
                i, t->Tcount, t->Xcount);
    }

    close_file(fp);
}

struct_topo *allocate_memory_for_topos(long num)
{
    struct_topo *topos;

    topos = (struct_topo *) malloc((num + 1) * sizeof(struct_topo));

    if (topos == NULL)
        fatal("malloc failure for allocate_memory_for_topos()\n");

    return topos;
}

void free_topos(struct_topoArr * topo)
{
    long i;

    for (i = 1; i <= topo->num_topo; i++) {
        free_p_strtag(&topo->topos[i].str_tag);
        free_lvector(topo->topos[i].offset, 0, DUMMY);
    }

    free(topo->topos);
}

struct_topo *allocate_memory_for_one_topo(void)
{
    struct_topo *topo;

    topo = (struct_topo *) malloc(sizeof(struct_topo));

    if (topo == NULL)
        fatal("malloc failure for allocate_memory_for_one_topo()\n");

    return topo;
}

void free_one_topo(struct_topo * topo)
{
    if (topo->str_tag.str)
        free_p_strtag(&topo->str_tag);

    if (topo->offset)
        free_lvector(topo->offset, 0, DUMMY);

    free(topo);
}

void duplicate_topo(struct_topo * src, struct_topo * dst)
{
    long i;

    init_strtag(&src->str_tag, &dst->str_tag);

    dst->Tcount = src->Tcount;
    dst->Xcount = src->Xcount;

    dst->offset = lvector(0, dst->Tcount);  /* more than enough */
    copy_lvector(src->offset, 0, dst->Xcount - 1, dst->offset);

    dst->msize = src->msize;

    dst->num_parens = src->num_parens;
    for (i = 1; i <= dst->num_parens; i++) {
        dst->parens_bidx[i] = src->parens_bidx[i];
        dst->parens_eidx[i] = src->parens_eidx[i];
    }
}

long has_matched_parens(char *str)
{
    long i, k = 0;

    for (i = 0; i < (long) strlen(str); i++) {
        if (str[i] == '(')
            k++;
        else if (str[i] == ')')
            k--;
        if (k < 0)
            break;
    }

    return k == 0;
}

long fillup_parens(char *str, long limit, long *parens_bidx, long *parens_eidx)
{
    long i, j, k = 0, m, num_parens = 0, len = strlen(str);
    long *matched, tidx[BUF16] = { 0 };

    if (!has_matched_parens(str))  /* just to be sure */
        fatal("pattern <%s> has unmatched () pairs\n", str);

    for (i = 0; i < len; i++) {
        if (str[i] == '(')
            parens_bidx[++num_parens] = i;  /* 1-index */
        if (str[i] == ')')
            tidx[++k] = i;
        if (num_parens >= limit)
            fatal("<%s> contains more then <%s> () pairs\n", str, limit);
    }

    matched = lvector(0, len);

    for (i = 1; i <= num_parens; i++) {  /* match to ')' from left to right */
        for (j = tidx[i] - 1; j >= 0; j--)  /* closest matched '(' to left */
            if (str[j] == '(' && !matched[j])
                for (m = num_parens; m >= 1; m--)  /* right to left for most common cases */
                    if (parens_bidx[m] == j) {
                        matched[j] = TRUE;
                        parens_eidx[m] = tidx[i];
                        goto NEXT_RIGHT_PAREN;
                    }
      NEXT_RIGHT_PAREN:
        do_nothing();
    }

    free_lvector(matched, 0, DUMMY);

    debug_parenthesis(str, num_parens, parens_bidx, parens_eidx);

    return num_parens;
}

/* str[] must be in upper case, containing ACGT, and no U */
void match_parens_str(long num_parens, long *parens_bidx, long *parens_eidx, char *str)
{
    char c;
    long i;

    /* no G-U Wobble check yet */
    for (i = 1; i <= num_parens; i++) {
        c = str[parens_bidx[i]];
        if (c == 'A')
            str[parens_eidx[i]] = 'T';
        else if (c == 'C')
            str[parens_eidx[i]] = 'G';
        else if (c == 'G')
            str[parens_eidx[i]] = 'C';
        else if (c == 'T')
            str[parens_eidx[i]] = 'A';
    }
}

/* ACGT -- 0123; thus ( + ) = 3; for Wobble, G + T = 5 */
long with_matched_parens_bin(struct_topo * t, char *b)
{
    long i, k;

    for (i = 1; i <= t->num_parens; i++) {
        k = b[t->parens_bidx[i]] + b[t->parens_eidx[i]];
        if (k != 3 && (!Gvars.WOBBLE || k != 5))
            return FALSE;  /* with any unmatched () */
#if 0  /* the following is more detailed */
        if (k == 3)  /* Waton-Crick pair: A-T or G-C */
            continue;
        if (Gvars.WOBBLE && k == 5)  /* Wobble G-T pair */
            continue;
        return FALSE;  /* with any unmatched () */
#endif
    }

    return TRUE;  /* all ()s matched */
}
