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

void show_version()
{
    fprintf(stderr, "%s\nPlease join us at URL: http://reducesuite.bussemakerlab.org/\n",
            VERSION_NUMBER);
    exit(0);
}

void get_sysdir(char *dstdir, char *dirname, char *filename)
/* note: 'dstdir' has an ending slash (/) */
{
    char *temp, str[BUF512];

    if (exist_file(filename))  /* check current directory */
        strcpy(dstdir, "./");  /* current directory */

    else if ((temp = getenv("REDUCE_SUITE")) != NULL) {
        strcpy(str, temp);
        add_end_slash(str);
        sprintf(dstdir, "%s%s/", str, dirname);

    } else if ((temp = getenv("HOME")) != NULL ||  /* Unix */
               (temp = getenv("HOMEDRIVE")) != NULL) {  /* Windows: HOMEPATH? */
        strcpy(str, temp);
        add_end_slash(str);
        sprintf(dstdir, "%sREDUCE-Suite-v2.2/%s/", str, dirname);

    } else
        fatal("can't locate directory with file [%s]\n", filename);
}

void get_sysfile(char *fullfile, char *dirname, char *basefile)
{
    get_sysdir(fullfile, dirname, basefile);
    strcat(fullfile, basefile);

    if (!exist_file(fullfile))
        fatal("can't locate '%s' --\n\tDid you run REDUCE_Suite_setup?\n", fullfile);
}

void display_helpfile(char *helpfile)
{
    char fullfile[BUF512], *p0;
    long num_char = 82;
    FILE *fp;

    get_sysfile(fullfile, "doc", helpfile);
    fp = open_file(fullfile, "r");

    print_sep(stderr, '=', num_char);
    while ((p0 = my_getline(fp)) != NULL) {
        if (!is_comment_line(p0))
            fprintf(stderr, "%s\n", p0);
        free(p0);
    }

    fprintf(stderr, "\n"
            "********************************** RESOURCES **********************************\n"
            " Questions? bug reports? or suggestions? Indeed, if you have *any* REDUCE Suite\n"
            "    related issues, please join us at http://reducesuite.bussemakerlab.org/\n"
            "             We strive to give you a quick and concrete response.\n");

    print_sep(stderr, '=', num_char);

    close_file(fp);

    exit(0);
}

long read_pwm(char *pwmfile, double *w)
{
    char *p0, *line, *items[BUF512], str[BUF512], msg[BUF512];
    double dval;
    long i, nitem, NUM_REQUIRED_ITEM = 4, idx = 0, num_nt = 0;

    FILE *fp;

    fp = open_file(pwmfile, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original address value of p0 */

        if (!is_skip_line(line)) {
            strcpy(str, line);
            num_nt++;

            nitem = item_list(line, items, NUM_REQUIRED_ITEM, WSPACES);
            if (nitem != NUM_REQUIRED_ITEM) {
                sprintf(msg, "<%s> not in PWM format: [%s] contains %ld =/= %ld items\n",
                        pwmfile, str, nitem, NUM_REQUIRED_ITEM);
                log_msg_exit(msg);
            }

            for (i = 1; i <= NUM_REQUIRED_ITEM; i++) {
                dval = cvt2double(items[i]);
                if (dval > XBIG_CUTOFF) {
                    sprintf(msg, "<%s> not in PWM format: [%s] contains invalid W %s\n",
                            pwmfile, str, items[i]);
                    log_msg_exit(msg);
                }
                w[idx++] = dval;
            }
        }
        free(p0);
    }

    close_file(fp);

    return num_nt;
}

void populate_pd_struct_from_psam(long num_nt, double *psam, long strand0, double p_value,
                                  struct_psam * psam_data)
{
    char motif[BUFBIG];
    long strand, nW;

    get_psam_optimal_seq(psam, num_nt, motif);

    strand = get_motif_strand(strand0);
    motif2psam_struct(motif, strand, p_value, psam_data);

    nW = psam_data->psam_length * NUM_BASE4 - 1;  /* 0-indexed */
    copy_dvector(psam, 0, nW, psam_data->psam);
}

void debug_parenthesis(char *pat, long num_parens, long *parens_bidx, long *parens_eidx)
{
    if (Gvars.VERBOSITY == VDEBUG && num_parens) {
        long i;

        fprintf(Gvars.RUNLOG, "[debug] <%s> has %ld matched ()s:\n", pat, num_parens);
        for (i = 1; i <= num_parens; i++)
            fprintf(Gvars.RUNLOG, "\t%ld\t%ld-%ld\n", i, parens_bidx[i], parens_eidx[i]);
    }
}

void write2xml(char *filename, args_reduce * args, struct_data * tdat,
               long psam_num, struct_seed * smotifs, struct_fitpars * lfpars)
{
    FILE *fp;

    fp = open_file(filename, "w");

    fprintf(fp, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(fp, "<reduce_suite version=\"2.0\">\n");

    xml_parameters(fp, args, tdat->ncol, psam_num);
    xml_experiment(fp, tdat, psam_num, lfpars);
    xml_psam(fp, psam_num, smotifs);

    fprintf(fp, "</reduce_suite>\n");

    close_file(fp);
}

void xml_parameters(FILE * fp, args_reduce * args, long expt_num, long psam_num)
{
    char *stnd_msg, *topo_str;

    stnd_msg = strand2msg(args->strand);

    fprintf(fp, "   <parameters>\n");
    fprintf(fp, "      <measurement_file>%s</measurement_file>\n",
            basename(args->measfile));
    fprintf(fp, "      <number_of_experiments>%ld</number_of_experiments>\n", expt_num);
    fprintf(fp, "      <sequence_file>%s</sequence_file>\n", basename(args->seqfile));

    if (is_empty_string(args->dictfile))
        topo_str = is_empty_string(args->topo) ? args->topo_list : args->topo;
    else
        topo_str = args->dictfile;
    fprintf(fp, "      <topology>%s</topology>\n", basename(topo_str));

    fprintf(fp, "      <min_counts>%ld</min_counts>\n", Gvars.misc.MIN_COUNTS);
    fprintf(fp, "      <directionality>%s</directionality>\n", stnd_msg);

    fprintf(fp, "      <stop_criteria>\n");
    fprintf(fp, "         <max_motifs>%ld</max_motifs>\n", args->max_motif);
    fprintf(fp, "         <max_pvalue>%g</max_pvalue>\n", args->p_value);
    fprintf(fp, "      </stop_criteria>\n");

    fprintf(fp, "      <number_of_psams>%ld</number_of_psams>\n", psam_num);
    fprintf(fp, "   </parameters>\n");

    free_cvector(stnd_msg, 0, DUMMY);
}

void xml_experiment(FILE * fp, struct_data * tdat, long psam_num, struct_fitpars * lfpars)
{
    long i;

    fprintf(fp, "\n   <experiment_summary>\n");

    for (i = 1; i <= tdat->ncol; i++) {
        fprintf(fp, "      <experiment expt_id=\"%ld\">\n", i);
        xml_fitpars(fp, tdat, psam_num, i, lfpars);
        fprintf(fp, "      </experiment>\n");

        if (i < tdat->ncol)
            fprintf(fp, "\n");
    }

    fprintf(fp, "   </experiment_summary>\n");
}

void xml_fitpars(FILE * fp, struct_data * tdat, long psam_num, long icol,
                 struct_fitpars * lfpars)
{
    char *p = "         ";
    long i, k = 0, nrow = tdat->nrow;
    double r, prob, z, tval;
    double *okdat, *okfit, *okres;
    struct_fitpars *par = &lfpars[icol];

    okdat = dvector(1, nrow);
    okfit = dvector(1, nrow);
    okres = dvector(1, nrow);

    for (i = 1; i <= nrow; i++)
        if (tdat->data[icol][i] < XBIG_CUTOFF) {
            k++;
            okdat[k] = tdat->data[icol][i];
            okres[k] = tdat->resid[icol][i];
            okfit[k] = okdat[k] - okres[k];
        }

    NRC_pearsn(okfit, okdat, k, &r, &prob, &z, &tval);

    fprintf(fp, "%s<summary>\n", p);
    fprintf(fp, "%s   <description>%s</description>\n", p, tdat->col_names[icol]);
    fprintf(fp, "%s   <mean>%+g</mean>\n", p, ave_dvector(okdat, 1, k));
    fprintf(fp, "%s   <var>%g</var>\n", p, var_dvector(okdat, 1, k));
    fprintf(fp, "%s   <dimensionality>%ld</dimensionality>\n", p, k);
    fprintf(fp, "%s</summary>\n\n", p);

    fprintf(fp, "%s<multivariate_fit>\n", p);
    fprintf(fp, "%s   <rsquared>%g</rsquared>\n", p, r * r);
    fprintf(fp, "%s   <fit_var>%g</fit_var>\n", p, var_dvector(okfit, 1, k));
    fprintf(fp, "%s   <res_var>%g</res_var>\n", p, var_dvector(okres, 1, k));

    i = 0;  /* for intercept */
    fprintf(fp, "%s   <intercept>\n", p);
    xml_FtP_value(fp, par->fval[i], par->tval[i], par->pval[i]);
    fprintf(fp, "%s   </intercept>\n", p);

    for (i = 1; i <= psam_num; i++) {
        fprintf(fp, "%s   <slope psam_id=\"%ld\">\n", p, i);
        xml_FtP_value(fp, par->fval[i], par->tval[i], par->pval[i]);
        fprintf(fp, "%s   </slope>\n", p);
    }

    fprintf(fp, "%s</multivariate_fit>\n", p);

    free_dvector(okdat, 1, DUMMY);
    free_dvector(okfit, 1, DUMMY);
    free_dvector(okres, 1, DUMMY);
}

void xml_FtP_value(FILE * fp, double f, double t, double p)
{
    char *leading_space = "               ";

    fprintf(fp, "%s<coeff>%+g</coeff>\n", leading_space, f);
    fprintf(fp, "%s<t_value>%+g</t_value>\n", leading_space, t);
    fprintf(fp, "%s<p_value>%g</p_value>\n", leading_space, p);
}

void xml_psam(FILE * fp, long psam_num, struct_seed * smotifs)
{
    char base;
    long i, idx, j, k;
    struct_seed *seed;

    fprintf(fp, "\n   <psam_summary>\n");

    for (i = 1; i <= psam_num; i++) {
        seed = &smotifs[i];

        fprintf(fp, "      <psam psam_id=\"%ld\">\n", i);
        fprintf(fp, "         <seed_motif>%s</seed_motif>\n", seed->full);
        fprintf(fp, "         <optimal_sequence>%s</optimal_sequence>\n", seed->optimal);
        fprintf(fp, "         <directionality>%s</directionality>\n", seed->stnd_msg);
        fprintf(fp, "         <affinities>\n");

        idx = 0;
        for (j = 1; j <= seed->topo->Tcount; j++) {
            fprintf(fp, "            <base_pair pos=\"%ld\">\n", j);
            for (k = 0; k < NUM_BASE4; k++) {
                base = BASES[k];
                fprintf(fp, "               <%c>%+.6g</%c>\n", base, seed->psam[idx + k],
                        base);
            }

            idx += NUM_BASE4;
            fprintf(fp, "            </base_pair>\n");
        }

        fprintf(fp, "         </affinities>\n");
        fprintf(fp, "      </psam>\n");
        if (i < psam_num)
            fprintf(fp, "\n");
    }

    fprintf(fp, "   </psam_summary>\n");
}
