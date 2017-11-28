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

void output_html_spar_fitpars(long psam_num, struct_data * tdat,
                              struct_seqArr * seqs, struct_seed * seed)
/* output seed-psam univariate fit parameters for HTML output */
{
    char str[BUF512];
    long i, num_char = 95;
    double *y, *counts;

    struct_ls5fit fitpars;

    print_sep(Gvars.PRGLOG, '=', num_char);
    fprintf(Gvars.PRGLOG, "List of PSAM#%ld univariate fit parameters against each"
            " experiment\n", psam_num);
    fprintf(Gvars.PRGLOG, "## expt_info\tnum\tr2\tseed\tF-value\tt-value\tP-value\n");

    counts = dvector(1, seqs->num_seqs);

    for (i = 1; i <= tdat->ncol; i++) {
        sprintf(str, "%ld\t%s", i, tdat->col_names[i]);

        y = tdat->resid[i];  /* with residual */
        get_psam_counts(seqs, seed->topo, seed->strand, seed->w, counts);
        get_y_vs_x_lsfitpars(tdat->nrow, counts, y, &fitpars);

        fprintf(Gvars.PRGLOG, "%s\t%ld\t%g\t%s\t%+g\t%+g\t%g\n", str, fitpars.num,
                fitpars.r2, (i == seed->col_idx) ? "yes" : "no", fitpars.b,
                fitpars.b_tval, fitpars.b_pval);
    }

    print_sep(Gvars.PRGLOG, '=', num_char);

    free_dvector(counts, 1, DUMMY);
}

void read_spar_fitpars(long psam_num, long num_expts, char *logfile, struct_spar * sinfo)
{
    char *p0, *line, str[BUF512], dup[BUF512], *items[BUF512];
    long i, k, nitem, NUM_REQUIRED_ITEM = 8;

    FILE *fp;

    sprintf(str, "List of PSAM#%ld univariate fit parameters", psam_num);

    fp = open_file(logfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */

        if (str_pmatch(line, str)) {
            free(p0);
            p0 = my_getline(fp);  /* skip one line */
            free(p0);

            for (i = 1; i <= num_expts; i++) {
                p0 = my_getline(fp);
                line = trim(p0);
                strcpy(dup, line);

                nitem = item_list(line, items, NUM_REQUIRED_ITEM, WSPACES);
                if (nitem != NUM_REQUIRED_ITEM)
                    fatal("wrong # of entries in log file <%s: %s>\n", logfile, dup);

                k = cvt2long(items[1]);
                if (k == LONG_MAX)
                    fatal("wrong expt_col# in log file <%s: %s>\n", logfile, dup);

                sinfo[i].experiment_column = k;
                strcpy(sinfo[i].experiment_name, items[2]);

                k = cvt2long(items[3]);
                if (k == LONG_MAX)
                    fatal("wrong ok# in log file <%s: %s>\n", logfile, dup);

                sinfo[i].num = k;
                sinfo[i].r2 = cvt2double(items[4]);

                if (is_equal_string(items[5], "yes"))
                    sinfo[i].is_seed = 1;
                else if (is_equal_string(items[5], "no"))
                    sinfo[i].is_seed = 0;
                else
                    fatal("wrong yes/no info in log file <%s: %s>\n", logfile, dup);

                sinfo[i].F = cvt2double(items[6]);
                sinfo[i].t = cvt2double(items[7]);
                sinfo[i].P = cvt2double(items[8]);

                free(p0);
            }

            if (sinfo[num_expts].experiment_column != num_expts)
                fatal("wrong # of experiments in log file <%s: %s>\n", logfile, dup);

            break;
        }
        free(p0);
    }

    close_file(fp);
}

double get_valFtP(long idx, double F, double t, double P)
{
    if (idx == FVAL_IDX)
        return F;
    else if (idx == TVAL_IDX)
        return t;
    else
        return P;
}

void print_valFtP(FILE * fp, long idx, double F, double t, double P)
{
    if (idx == FVAL_IDX)
        fprintf(fp, "\t%+g", F);
    else if (idx == TVAL_IDX)
        fprintf(fp, "\t%+g", t);
    else
        fprintf(fp, "\t%g", P);  /* no + sign */
}

void write_FtP(char *type, char *outdir, char *expr_bname, struct_data * tdat,
               long psam_num, struct_fitpars * lfpars)
/* run "cluster -e 1 -g 1 -f model_coeff.tsv" for display with Java TreeView */
{
    char filename[BUF512];
    char *model_files[] = { MODEL_FILE_LIST };
    long i, j, k, num_files = sizeof(model_files) / sizeof(model_files[0]);

    FILE *fp;

    for (i = 0; i < num_files; i++) {
        sprintf(filename, "%s/model_%s", outdir, model_files[i]);
        fp = open_file(filename, "w");

        fprintf(fp, "### measurement_file=%s\n", expr_bname);
        fprintf(fp, "EXPID\tintercept");
        for (k = 1; k <= psam_num; k++)
            fprintf(fp, "\t%s%ld", type, k);
        fprintf(fp, "\n");

        for (j = 1; j <= tdat->ncol; j++) {  /* loop over experiments */
            fprintf(fp, "%s[%ld]", tdat->col_names[j], j);

            for (k = 0; k <= psam_num; k++)  /* 0 for intercept */
                print_valFtP(fp, i, lfpars[j].fval[k], lfpars[j].tval[k],
                             lfpars[j].pval[k]);
            fprintf(fp, "\n");
        }

        close_file(fp);
    }
}

void write2raw(char *filename, args_reduce * args, struct_data * tdat,
               long psam_num, struct_seed * smotifs, struct_fitpars * lfpars)
{
    long i, j, num;
    double dval, r, prob, z, tval, *org, *fit, *res;
    FILE *fp;

    fp = open_file(filename, "w");
    fprintf(fp, "sequence_file=%s\n", basename(args->seqfile));
    fprintf(fp, "measurement_file=%s\n", basename(args->measfile));

    fprintf(fp, "\nnumber_of_psams=%ld\n", psam_num);
    for (i = 1; i <= psam_num; i++)
        fprintf(fp, "\t%s|%ld\n", smotifs[i].col_name, smotifs[i].col_idx);

    fprintf(fp, "\nnumber_of_experiments=%ld\n", tdat->ncol);
    for (i = 1; i <= tdat->ncol; i++)
        fprintf(fp, "\t%s|%ld\n", tdat->col_names[i], i);

    fprintf(fp,
            "\n### The following is the multivariate linear ls-fit parameters\n"
            "### c0 means intercept, c1 for the 1st PSAM and so on\n"
            "### in the format of ci=F|t|P for coef, t-value, and P-value\n");
    for (i = 1; i <= tdat->ncol; i++) {
        fprintf(fp, "%s|%ld", tdat->col_names[i], i);
        for (j = 0; j <= psam_num; j++)
            fprintf(fp, "\tc%ld=%+g|%+g|%g", j, lfpars[i].fval[j], lfpars[i].tval[j],
                    lfpars[i].pval[j]);
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n### output for overall fit: raw_measurement_value = a + b * y_fit\n");
    org = dvector(1, tdat->nrow);
    fit = dvector(1, tdat->nrow);
    res = dvector(1, tdat->nrow);

    for (i = 1; i <= tdat->ncol; i++) {
        num = 0;
        for (j = 1; j <= tdat->nrow; j++) {
            dval = tdat->data[i][j];
            if (dval < XBIG_CUTOFF) {
                num++;
                org[num] = dval;
                res[num] = tdat->resid[i][j];
                fit[num] = dval - res[num];
            }
        }

        NRC_pearsn(fit, org, num, &r, &prob, &z, &tval);
        fprintf(fp, "%s|%ld\tnum=%ld\tr2=%g", tdat->col_names[i], i, num, r * r);
        fprintf(fp, "\tmean=%+g\tvar=%g", ave_dvector(org, 1, num),
                var_dvector(org, 1, num));
        fprintf(fp, "\tfit-var=%g", var_dvector(fit, 1, num));
        fprintf(fp, "\tres-var=%g\n", var_dvector(res, 1, num));
    }

    free_dvector(org, 1, DUMMY);
    free_dvector(fit, 1, DUMMY);
    free_dvector(res, 1, DUMMY);

    close_file(fp);
}

void write_psams_list(long motif_num, struct_seed * smotifs, char *outdir)
{
    char filename[BUF512];
    long i;
    FILE *fp;

    sprintf(filename, "%s/%s", outdir, "psams.list");
    fp = open_file(filename, "w");

    for (i = 1; i <= motif_num; i++)
        fprintf(fp, "psam_%3.3ld.xml\t%s|%ld\n", i, smotifs[i].col_name,
                smotifs[i].col_idx);

    close_file(fp);
}

void write_motifs_list(long motif_num, struct_seed * smotifs, char *outdir)
{
    char filename[BUF512];
    long i;
    FILE *fp;

    sprintf(filename, "%s/%s", outdir, "motifs.list");
    fp = open_file(filename, "w");

    for (i = 1; i <= motif_num; i++)
        fprintf(fp, "%-20s\t%s|%ld\n", smotifs[i].full, smotifs[i].col_name,
                smotifs[i].col_idx);

    close_file(fp);
}

void write_list_resid_predict(long psam_num, struct_seed * smotifs,
                              char *outdir, struct_data * tdat)
{
    char filename[BUF512];

    write_psams_list(psam_num, smotifs, outdir);
    write_motifs_list(psam_num, smotifs, outdir);

    sprintf(filename, "%s/%s", outdir, RESIDUAL_FILE);
    write_tdat(tdat, filename, tdat->resid, NULL);

    get_predicted(tdat);
    sprintf(filename, "%s/%s", outdir, PREDICT_FILE);
    write_tdat(tdat, filename, tdat->resid, NULL);
}

void output_mvar_fitpars(struct_data * tdat, long num, char *type, char *outdir,
                         struct_ls5fit * fitpars, struct_fitpars * lfpars)
{
    char filename[BUF512], *model_files[] = { MODEL_FILE_LIST };
    long num_files = sizeof(model_files) / sizeof(model_files[0]);
    long i, j, k;
    FILE *fp;

    for (i = 0; i < num_files; i++) {
        sprintf(filename, "%s/transfactivity_mvar_%s", outdir, model_files[i]);
        fp = open_file(filename, "w");

        fprintf(fp, "EXPID\tr2\tintercept");
        for (k = 1; k <= num; k++)
            fprintf(fp, "\t%s%ld", type, k);
        fprintf(fp, "\n");

        for (j = 1; j <= tdat->ncol; j++) {
            fprintf(fp, "%s[%ld]-(%ld)\t%g", tdat->col_names[j], j, fitpars[j].num,
                    fitpars[j].r2);

            for (k = 0; k <= num; k++)  /* 0 for intercept */
                print_valFtP(fp, i, lfpars[j].fval[k], lfpars[j].tval[k],
                             lfpars[j].pval[k]);
            fprintf(fp, "\n");
        }

        close_file(fp);
    }
}

static void write_FtP_HTML(FILE * fp, double F, double t, double P)
{
    char pn[BUF512];

    strcpy(pn, (F < 0.0) ? "negative" : "positive");

    fprintf(fp, "<td class=\"%s\">%+g</td>\n", pn, F);
    fprintf(fp, "<td class=\"%s tvalue\">%.2f</td>\n", pn, t);
    fprintf(fp, "<td>%g</td>\n", P);
}

void write_separate_FtP_HTML(FILE * fp, long intercept, long idx, double F, double t,
                             double P, double tc)
{
    char str[BUF512] = "tvalue";
    double dval = get_valFtP(idx, F, t, P);

    if (idx == PVAL_IDX)
        fprintf(fp, "<td>%.3g</td>\n", dval);

    else if (intercept)
        fprintf(fp, "<td>%+.3g</td>\n", dval);

    else if (idx == TVAL_IDX) {
        if (dval > tc)
            strcat(str, " positive");
        else if (dval < -tc)
            strcat(str, " negative");

        fprintf(fp, "<td class=\"%s\">%.2f</td>\n", str, dval);

    } else {
        strcpy(str, (dval < 0.0) ? "negative" : "positive");
        fprintf(fp, "<td class=\"%s\">%+g</td>\n", str, dval);
    }
}

static void write_unihead(FILE * fp)
{
    fprintf(fp, "<thead><tr>\n"
            "<th>&nbsp;</th><th>r<sup>2</sup></th>\n"
            "<th>intercept</th><th>t-value</th><th>P-value</th>\n"
            "<th class=\"markit\">slope</th><th class=\"markit\">t-value</th>"
            "<th>P-value</th>\n" "</tr></thead>\n");
}

static void write_r2head(FILE * fp)
{
    fprintf(fp, "<thead><tr>\n"
            "<th>&nbsp;</th><th>r<sup>2</sup></th>\n"
            "<th>df</th><th>t-value</th><th>P-value</th>\n" "</tr></thead>\n");
}

void open_table(FILE * fp, char *class)
{
    char str[BUF512];

    fprintf(fp, "<div class=\"table-wrapper\">\n");

    if (is_empty_string(class))
        strcpy(str, class);
    else
        sprintf(str, " %s", class);

    fprintf(fp, "<table cellspacing=\"0\" summary=\"\"%s>\n", str);
}

void close_table(FILE * fp)
{
    fprintf(fp, "</table></div> <!-- table-wrapper -->\n");
}

void html_recap_common(FILE * fp, char *options, char *seqfile, char *measfile,
                       struct_data * tdat)
{
    long i, ncol = tdat->ncol;

    fprintf(fp,
            "<li><a href=\"%s\">%s</a>, which contains all\n"
            "the necessary information so that the current job can be exactly\n"
            "reproduced. The info includes the version number of the %s package,\n"
            "the directory where %s was run, and all the command line options.</li>\n\n",
            options, options, MATRIX, Gvars.PROGNAME);

    fprintf(fp,
            "<li><a href=\"clean_seq.fasta\">%s</a>, the cleaned up sequence file\n"
            " in <a href=\"http://en.wikipedia.org/wiki/Fasta_format\">FASTA\n"
            "format</a>. The sequence ID is taken as the first word following\n"
            "the leading &gt; and is converted into lowercase. Each sequence is\n"
            "also converted into lowercase, and broken into fragments at\n"
            "non-'ACGT' bases. The cleaned up file is sorted into alphabetical\n"
            "order by the sequence ID.</li>\n", seqfile);

    fprintf(fp,
            "<li><a href=\"clean_data.tsv\">%s</a>, the measurement file, which"
            " contains %ld experiment%s\n", measfile, ncol, (ncol == 1) ? "" : "s");

    fprintf(fp, "<ol>\n");
    for (i = 1; i <= ncol; i++)
        fprintf(fp, "<li>%s</li>\n", basename(tdat->col_names[i]));
    fprintf(fp, "</ol></li>\n");
}

void set_transfactivity_univar_htmlHeader(FILE * fp, char *seqfile, struct_data * tdat,
                                          char *measfile)
{
    char options[BUF512];

    sprintf(options, "%s.opt", Gvars.PROGNAME);

    fprintf(fp, "<h2 id=\"recap\">Recap of used options</h2>\n");
    fprintf(fp, "<ul>\n");
    html_recap_common(fp, options, seqfile, measfile, tdat);
    fprintf(fp, "</ul>\n");
}

void set_transfactivity_htmlHeader(FILE * fp, long num, struct_tag * str_tags,
                                   char *SUBMENU, char *type, char *seqfile,
                                   struct_data * tdat, char *measfile, char *resid_file)
{
    char *id, options[BUF512];
    long i;

    sprintf(options, "%s.opt", Gvars.PROGNAME);

    fprintf(fp, "<h2 id=\"recap\">Recap of used options</h2>\n");
    fprintf(fp, "\n%s\n", SUBMENU);
    fprintf(fp, "<ul>\n");
    html_recap_common(fp, options, seqfile, measfile, tdat);

    fprintf(fp, "<li><a href=\"%s\">%s</a>, the residual</li>\n", resid_file, resid_file);
    fprintf(fp, "<li>In this %s run, you've used %ld %s%s:\n",
            Gvars.PROGNAME, num, type, (num == 1) ? "" : "s");

    fprintf(fp, "<ol>\n");
    for (i = 1; i <= num; i++) {
        id = str_tags[i].str;
        if (str_pmatch(type, "PSAM"))
            fprintf(fp, "<li><a href=\"%s\">%s</a></li>\n", id, basename(id));
        else
            fprintf(fp, "<li>%s</li>\n", basename(id));
    }

    fprintf(fp, "</ol></li>\n");
    fprintf(fp, "</ul>\n");
}

void set_overall_fit(FILE * fp, long num, char *SUBMENU, struct_data * tdat,
                     struct_ls5fit * fitpars, struct_fitpars * lfpars)
{
    long i;
    struct_ls5fit *fitp;

    fprintf(fp, "<h2 id=\"r2\">The overall fit r<sup>2</sup></h2>\n");
    fprintf(fp, "\n%s\n", SUBMENU);

    open_table(fp, "class=\"ftp\"");
    write_r2head(fp);

    fprintf(fp, "<tbody>\n");

    for (i = 1; i <= tdat->ncol; i++) {
        fitp = &fitpars[i];
        fprintf(fp, "<tr><td>%s[%ld]-(%ld)</td><td>%g</td>\n", tdat->col_names[i], i,
                lfpars[i].num, fitp->r2);
        fprintf(fp, "<td>%ld</td>", fitp->num - num - 1);  /* 1 for intercept */
        fprintf(fp, "<td class=\"%s\">", (fitp->b_tval < 0.0) ? "negative" : "positive");
        fprintf(fp, "%.2f</td><td>%.3g</td>", fitp->b_tval, fitp->b_pval);
        fprintf(fp, "</tr>\n");
    }

    fprintf(fp, "</tbody>\n");

    close_table(fp);
}

void mvar_FtP_transposed(FILE * fp, long FtP_idx, long ncol, long num, long intercept,
                         struct_fitpars * lfpars, struct_tag * str_tags, double tc)
{
    char *id;
    long i, j;

    fprintf(fp, "<thead><tr><th>&nbsp;</th>\n");
    for (j = 1; j <= ncol; j++)
        fprintf(fp, "<th>expt_%ld</th>", j);
    fprintf(fp, "</tr></thead>\n");

    fprintf(fp, "<tbody>\n");

    for (j = 0; j <= num; j++) {
        if (j == 0) {
            if (intercept) {
                fprintf(fp, "<tr>");
                fprintf(fp, "<td>intercept</td>");

                for (i = 1; i <= ncol; i++)
                    write_separate_FtP_HTML(fp, TRUE, FtP_idx, lfpars[i].fval[j],
                                            lfpars[i].tval[j], lfpars[i].pval[j], tc);
                fprintf(fp, "</tr>\n\n");
            }

        } else {
            fprintf(fp, "<tr>");
            id = str_tags[j].str;
            fprintf(fp, "<td>%s</td>", basename(id));

            for (i = 1; i <= ncol; i++)
                write_separate_FtP_HTML(fp, FALSE, FtP_idx, lfpars[i].fval[j],
                                        lfpars[i].tval[j], lfpars[i].pval[j], tc);
            fprintf(fp, "</tr>\n\n");
        }
    }

    fprintf(fp, "</tbody>\n");
}

void mvar_FtP_normal(FILE * fp, long FtP_idx, struct_data * tdat,
                     long num, struct_fitpars * lfpars, char *type, double tc)
{
    long i, j;

    fprintf(fp, "<thead><tr><th>&nbsp;</th><th>intercept</th>\n");
    for (j = 1; j <= num; j++)
        fprintf(fp, "<th>%s%ld</th>", type, j);
    fprintf(fp, "</tr></thead>\n");

    fprintf(fp, "<tbody>\n");

    for (i = 1; i <= tdat->ncol; i++) {
        fprintf(fp, "<tr>");
        fprintf(fp, "<td>%s[%ld]-(%ld)</td>", tdat->col_names[i], i, lfpars[i].num);

        for (j = 0; j <= num; j++)
            write_separate_FtP_HTML(fp, j == 0, FtP_idx, lfpars[i].fval[j],
                                    lfpars[i].tval[j], lfpars[i].pval[j], tc);
        fprintf(fp, "</tr>\n\n");
    }

    fprintf(fp, "</tbody>\n");
}

void output_mvar_fitparsHTML(struct_data * tdat, long num, char *type,
                             struct_ls5fit * fitpars, struct_fitpars * lfpars,
                             struct_tag * str_tags, args_transfactivity * args)
{
    char filename[BUF512];
    char *outdir = args->outdir;
    char *resid_file = args->resid_file;
    char *seqfile = basename(args->seqfile);
    char *measfile = basename(args->measfile);
    char *SUBMENU = "\n<p class=\"cmenu\">\
\n<a href=\"#wrap\">Top</a> | \
\n<a href=\"#recap\">Recap</a> | \
\n<a href=\"#r2\">r<sup>2</sup></a> | \
\n<a href=\"#ftable\">F-value</a> | \
\n<a href=\"#ttable\">t-value</a> | \
\n<a href=\"#ptable\">P-value</a>\
</p>";
    double tc = args->tval_cutoff;
    long idx, ncol = tdat->ncol;
    FILE *fp;

    sprintf(filename, "%s/%s", outdir, "transfactivity_mvar_results.html");
    fp = open_file(filename, "w");

    set_htmlFileHeader(fp, outdir, args->copy, "Transfactivity &ndash; multivariate fit");

    set_transfactivity_htmlHeader(fp, num, str_tags, SUBMENU, type, seqfile, tdat,
                                  measfile, resid_file);

    set_overall_fit(fp, num, SUBMENU, tdat, fitpars, lfpars);

    log_msg("\nperforming multivariate fit");
    for (idx = 0; idx < 3; idx++) {
        if (idx == FVAL_IDX)
            fprintf(fp, "<h2 id=\"ftable\">Summary table of F-values</h2>");

        else if (idx == TVAL_IDX)
            fprintf(fp, "<h2 id=\"ttable\">Summary table of t-values</h2>");

        else  /* idx == PVAL_IDX */
            fprintf(fp,
                    "<h2 id=\"ptable\">Summary table of P-values (<em>no</em> "
                    "multiple test bonferroni correction performed)</h2>");

        fprintf(fp, "\n%s\n", SUBMENU);
        open_table(fp, "class=\"ftp\"");

        if (args->transpose)
            mvar_FtP_transposed(fp, idx, ncol, num, args->intercept, lfpars,
                                str_tags, tc);
        else
            mvar_FtP_normal(fp, idx, tdat, num, lfpars, type, tc);

        close_table(fp);
    }

    set_htmlFileFooter(fp);

    close_file(fp);
}

void output_univar(struct_data * tdat, long num, char *type, double **Xmtx,
                   struct_tag * str_tags, char *outdir)
{
    char filename[BUF512], *model_files[] = { MODEL_FILE_LIST };
    long num_files = sizeof(model_files) / sizeof(model_files[0]);
    long i, j, m;
    FILE *fp;
    struct_ls5fit fitpars;

    for (i = 0; i < num_files; i++) {
        sprintf(filename, "%s/transfactivity_uvar_%s", outdir, model_files[i]);
        fp = open_file(filename, "w");

        fprintf(fp, "%s\tEXPID\tr2\tintercept\tslope\n", type);
        for (m = 1; m <= num; m++) {  /* loop over each motif/PSAM */
            for (j = 1; j <= tdat->ncol; j++) {  /* loop over experiments */
                fprintf(fp, "%s", str_tags[m].str);
                get_y_vs_x_lsfitpars(tdat->nrow, Xmtx[m], tdat->data[j], &fitpars);

                fprintf(fp, "\t%s[%ld]-(%ld)", tdat->col_names[j], j, fitpars.num);
                fprintf(fp, "\t%g", fitpars.r2);

                if (i == FVAL_IDX)
                    fprintf(fp, "\t%+g\t%+g\n", fitpars.a, fitpars.b);
                else if (i == TVAL_IDX)
                    fprintf(fp, "\t%+g\t%+g\n", fitpars.a_tval, fitpars.b_tval);
                else
                    fprintf(fp, "\t%g\t%g\n", fitpars.a_pval, fitpars.b_pval);
            }

            if (m < num)
                fprintf(fp, "\n");
        }

        close_file(fp);
    }
}

void output_univarHTML(struct_data * tdat, long num, char *type, double **Xmtx,
                       struct_tag * str_tags, args_transfactivity * args)
{
    char filename[BUF512], *outdir = args->outdir;
    char *seqfile = basename(args->seqfile);
    char *measfile = basename(args->measfile);
    long m;
    FILE *fp;

    sprintf(filename, "%s/%s", outdir, "transfactivity_uvar_results.html");
    fp = open_file(filename, "w");

    set_htmlFileHeader(fp, outdir, args->copy, "Transfactivity &ndash; univariate fit");

    set_transfactivity_univar_htmlHeader(fp, seqfile, tdat, measfile);

    fprintf(fp, "<h2>List of all %s%s</h2>\n", type, (num == 1) ? "" : "s");
    fprintf(fp, "<ol>\n");
    for (m = 1; m <= num; m++)
        fprintf(fp, "<li><a href=\"#m%ld\">%s</a></li>\n", m, str_tags[m].str);
    fprintf(fp, "</ol>\n");

    fprintf(fp, "<h2>Detailed univariate fit parameters (<em>no</em> multiple"
            " test correction performed on p-value)</h2>\n");
    log_msg("\nperforming univariate fit");

    if (tdat->ncol > 1)
        output_univarHTML_multiple_expt(fp, tdat, num, type, Xmtx, str_tags);
    else
        output_univarHTML_one_expt(fp, tdat, num, type, Xmtx, str_tags);

    set_htmlFileFooter(fp);
    close_file(fp);
}

void output_univarHTML_multiple_expt(FILE * fp, struct_data * tdat, long num, char *type,
                                     double **Xmtx, struct_tag * str_tags)
{
    char msg[BUF512], *p;
    long m, i;
    struct_ls5fit fitpars;

    fprintf(fp, "<dl>\n");
    for (m = 1; m <= num; m++) {  /* loop over each motif/PSAM */
        sprintf(msg, "\t%s%ld\t%s", type, m, str_tags[m].str);
        log_msg(msg);

        fprintf(fp, "<dt id=\"m%ld\">%ld. ", m, m);
        p = str_tags[m].str;
        if (str_pmatch(type, "PSAM"))
            fprintf(fp, "<a href=\"%s\">%s</a>", p, p);
        else
            fprintf(fp, "<strong>%s</strong>", p);

        fprintf(fp, " <a href=\"#wrap\">[Top]</a></dt>\n");
        fprintf(fp, "<dd>\n");
        open_table(fp, "class=\"ftp\"");
        write_unihead(fp);

        fprintf(fp, "<tbody>\n");
        for (i = 1; i <= tdat->ncol; i++) {
            get_y_vs_x_lsfitpars(tdat->nrow, Xmtx[m], tdat->data[i], &fitpars);
            fprintf(fp, "<tr>");
            fprintf(fp, "<td>%s[%ld]-(%ld)</td>", tdat->col_names[i], i, fitpars.num);
            fprintf(fp, "<td>%g</td\n>", fitpars.r2);
            write_FtP_HTML(fp, fitpars.a, fitpars.a_tval, fitpars.a_pval);
            write_FtP_HTML(fp, fitpars.b, fitpars.b_tval, fitpars.b_pval);
            fprintf(fp, "</tr>\n");
        }
        fprintf(fp, "</tbody>\n");

        close_table(fp);
        fprintf(fp, "</dd>\n");
    }
    fprintf(fp, "</dl>\n");
}

void output_univarHTML_one_expt(FILE * fp, struct_data * tdat, long num, char *type,
                                double **Xmtx, struct_tag * str_tags)
{
    char msg[BUF512], *p;
    long m, i = 1;
    struct_ls5fit fitpars;

    open_table(fp, "class=\"ftp\"");

    write_unihead(fp);
    fprintf(fp, "<tbody>\n");

    for (m = 1; m <= num; m++) {  /* loop over each motif/PSAM */
        sprintf(msg, "\t%s%ld\t%s", type, m, str_tags[m].str);
        log_msg(msg);

        fprintf(fp, "<tr>");

        fprintf(fp, "<td id=\"m%ld\">%ld. ", m, m);
        p = str_tags[m].str;
        if (str_pmatch(type, "PSAM"))
            fprintf(fp, "<a href=\"%s\">%s</a>", p, p);
        else
            fprintf(fp, "<strong>%s</strong>", p);
        fprintf(fp, " <a href=\"#wrap\">[Top]</a></td>\n");

        get_y_vs_x_lsfitpars(tdat->nrow, Xmtx[m], tdat->data[i], &fitpars);
        fprintf(fp, "<td>%g</td\n>", fitpars.r2);
        write_FtP_HTML(fp, fitpars.a, fitpars.a_tval, fitpars.a_pval);
        write_FtP_HTML(fp, fitpars.b, fitpars.b_tval, fitpars.b_pval);
        fprintf(fp, "</tr>\n");
    }

    fprintf(fp, "</tbody>\n");

    close_table(fp);

    if (Gvars.VERBOSITY == VDEBUG) {
        FILE *fpc;  /* for checking purpose */

        fpc = open_file("univar_fit.dat", "w");
        for (m = 1; m <= num; m++)
            fprintf(fpc, "\t%s", str_tags[m].str);
        fprintf(fpc, "\ty\n");
        for (i = 1; i <= tdat->nrow; i++) {
            fprintf(fpc, "%s", tdat->ids[i]);
            for (m = 1; m <= num; m++)
                fprintf(fpc, "\t%g", Xmtx[m][i]);
            fprintf(fpc, "\t%g\n", tdat->data[1][i]);  /* one-experiment */
        }
        close_file(fpc);
    }
}

void set_htmlFileHeader(FILE * fp, char *outdir, long copy, char *pgname)
{
    char css[BUF512], js[BUF512], fileLoc[BUF512], *bname = "HTMLSummary";

    get_sysdir(fileLoc, "html", "HTMLSummary.css");
    sprintf(css, "%s%s.css", fileLoc, bname);
    sprintf(js, "%s%s.js", fileLoc, bname);
    copy_gifs_css_js(copy, fileLoc, outdir, css, js);

    fprintf(fp, "<!DOCTYPE html>\n<html lang=\"en\">\n\n");
    fprintf(fp, "<head>\n<meta charset=\"UTF-8\">\n");
    fprintf(fp,
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");

    fprintf(fp, "  <title>Summary Page: %s</title>\n", pgname);
    fprintf(fp, "  <link rel=\"stylesheet\" type=\"text/css\" href=\"%s\">\n", css);
    fprintf(fp, "  <script type=\"text/javascript\" src=\"%s\">\n", js);
    fprintf(fp, "  </script>\n");
    fprintf(fp, "</head>\n\n");

    get_currentTimeString(css);  /* 'css' contains a time string */

    fprintf(fp,
            "<body>\n"
            "  <div id=\"wrap\">\n"
            "  <div id=\"header\">\n"
            "    <h1>Results Summary Page for %s</h1>\n"
            "   <p>%s</p>\n  </div>\n\n", pgname, css);

    fprintf(fp, "<div id=\"content\">\n");
    fprintf(fp, "  <div class=\"gutter\">\n");
}

void set_htmlFileFooter(FILE * fp)
{
    char year[BUF512];

    get_currentYear(year);

    fprintf(fp, "  </div>  <!-- end .gutter -->\n");
    fprintf(fp, "</div>  <!-- end #content -->\n");
    fprintf(fp,
            "\n  <div id=\"footer\">\n"
            "    <p>Copyright &copy; %s <a href=\"http://www.bussemakerlab.org/\">\n"
            "      Harmen Bussemaker Laboratory</a>, Columbia University.\n"
            "      All rights reserved.</p>\n"
            "  </div>  <!-- end #footer -->\n"
            "</div>  <!-- end #wrap -->\n  </body>\n</html>\n", year);
}

void copy_gifs_css_js(long copy, char *fileLoc, char *outdir, char *css, char *js)
{
    char *files[] = {
        "HTMLSummary.css",
        "HTMLSummary.js",
        "content-bg.gif",
        "footer-bg.gif",
        "header-bg.gif",
        "li-bullet.gif"
    };
    char src[BUF512], dst[BUF512];
    long i, num_files = sizeof(files) / sizeof(files[0]);

    if (!copy)
        return;

    for (i = 0; i < num_files; i++) {
        sprintf(src, "%s%s", fileLoc, files[i]);
        sprintf(dst, "%s/%s", outdir, files[i]);
        copy_file(src, dst);

        if (strstr(files[i], ".css")) {
            sprintf(dst, "./%s", files[i]);
            strcpy(css, dst);  /* point to the correct location */

        } else if (strstr(files[i], ".js")) {
            sprintf(dst, "./%s", files[i]);
            strcpy(js, dst);
        }
    }
}

long is_matrix_reduce_run(char *outdir)
{
    char filename[BUF512];

    sprintf(filename, "%s/%s", outdir, "MatrixREDUCE.opt");

    return exist_file(filename);
}

void get_raw_data_file(char *outdir, char *filename)
{
    if (is_matrix_reduce_run(outdir))
        sprintf(filename, "%s/%s", outdir, MATRIXREDUCE_RAW);
    else
        sprintf(filename, "%s/%s", outdir, MOTIFREDUCE_RAW);
}

void get_default_html_file(char *outdir, char *filename)
/* Note: 'outdir' not included in 'filename' */
{
    if (is_matrix_reduce_run(outdir))
        strcpy(filename, MATRIXREDUCE_HTML);
    else
        strcpy(filename, MOTIFREDUCE_HTML);
}
