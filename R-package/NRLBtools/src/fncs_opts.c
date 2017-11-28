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

void extract_log_base(char *option, char *log_base)
{
    char str[BUF512];
    double dval;

    strcpy(log_base, "e");  /* default to natural log */

    if (!strchr(option, '='))
        return;

    get_strvalue(option, str, FALSE);
    lowerstr(str);
    if (str_pmatch(str, "n") || str_pmatch(str, "e"))  /* natural log */
        return;

    dval = cvt2double(str);
    if (dval > XBIG_CUTOFF)
        fatal("unrecognized -log=%s option [%s]\n", str, option);
    else if (dval <= 0)
        fatal("can't take log-base (%ld <= 0 [%s])\n", dval, option);
    else
        strcpy(log_base, str);
}

long extract_option_psam_motif(char *option, char *psam, char *psam_list, char *motif,
                               char *motif_list)
{
    if (str_pmatch(option, "-p=") || str_pmatch(option, "-psam="))
        get_strvalue(option, psam, TRUE);

    else if (str_pmatch(option, "-pl") || str_pmatch(option, "-psam_l"))
        get_strvalue(option, psam_list, TRUE);

    else if (str_pmatch(option, "-m=") || str_pmatch(option, "-motif="))
        get_strvalue(option, motif, FALSE);

    else if (str_pmatch(option, "-ml") || str_pmatch(option, "-motif_l"))
        get_strvalue(option, motif_list, TRUE);

    else
        return FALSE;

    return TRUE;
}

long extract_option_topo(char *option, char *topo, char *topo_list)
{
    if (str_pmatch(option, "-t=") || str_pmatch(option, "-topo="))
        get_strvalue(option, topo, FALSE);

    else if (str_pmatch(option, "-tl") || str_pmatch(option, "-topo_l"))
        get_strvalue(option, topo_list, TRUE);

    else
        return FALSE;

    return TRUE;
}

long extract_verbose_option(char *option)
{
    if (strchr(option, '=') != NULL)
        return get_lvalue(option, -1, BUF512);  /* -1 to repress any output */
    else
        return VERBOSE;  /* just -verb, use default value */
}

long extract_case_option(char *option)
{
    char msg[BUF512], str[BUF512];
    long is_upcase;

    get_strvalue(option, str, FALSE);
    lowerstr(str);

    if (str_pmatch(str, "l") || str_pmatch(str, "s") || str_pmatch(str, "0"))
        is_upcase = FALSE;  /* lower/small/0 */

    else if (str_pmatch(str, "u") || str_pmatch(str, "c") || str_pmatch(str, "1"))
        is_upcase = TRUE;  /* upper/captial/1 */

    else {
        sprintf(msg, "\tunrecognized option <%s>: taken as upper case", option);
        log_msg(msg);
        is_upcase = TRUE;
    }

    return is_upcase;
}

void write_option_case(FILE * fp, long case_setting)
{
    if (case_setting != UNSET_LVAL)
        fprintf(fp, "\t-case=%s \\\n", (case_setting) ? "upper" : "lower");
}

void write_option_outdir_runlog(FILE * fp, char *outdir, char *runlog)
{
    if (!is_equal_string(outdir, ".") && !is_equal_string(outdir, "./"))
        fprintf(fp, "\t-outdir=%s \\\n", outdir);

    if (!is_equal_string(runlog, "stderr"))
        fprintf(fp, "\t-runlog=%s \\\n", runlog);
}

void write_option_psam_motif(FILE * fp, char *psam, char *psam_list, char *motif,
                             char *motif_list)
{
    if (!is_empty_string(psam))
        fprintf(fp, "\t-psam=%s \\\n", psam);

    if (!is_empty_string(psam_list))
        fprintf(fp, "\t-psam_list=%s \\\n", psam_list);

    if (!is_empty_string(motif))
        fprintf(fp, "\t-motif=%s \\\n", motif);

    if (!is_empty_string(motif_list))
        fprintf(fp, "\t-motif_list=%s \\\n", motif_list);
}

void write_option_topo(FILE * fp, char *topo, char *topo_list)
{
    if (!is_empty_string(topo))
        fprintf(fp, "\t-topo='%s' \\\n", topo);

    if (!is_empty_string(topo_list))
        fprintf(fp, "\t-topo_list=%s \\\n", topo_list);
}

void write_option_expand(FILE * fp, char *expand, char *expand_list)
{
    if (!is_empty_string(expand))
        fprintf(fp, "\t-expand=%s \\\n", expand);

    if (!is_empty_string(expand_list))
        fprintf(fp, "\t-expand_list=%s \\\n", expand_list);
}

long extract_strand_option(char *option)
{
    char str[BUF512];
    long strand = 0;

    get_strvalue(option, str, FALSE);
    lowerstr(str);

    if (str_pmatch(str, "l") || str_pmatch(str, "f") ||  /* -leading/-forward */
        str_pmatch(str, "1") || str_pmatch(str, "+1"))
        strand = 1;

    else if (str_pmatch(str, "r") ||  /* -reverse */
             str_pmatch(str, "c") ||  /* -complementary */
             str_pmatch(str, "-1"))
        strand = -1;

    else if (str_pmatch(str, "b") ||  /* -both */
             str_pmatch(str, "2") || str_pmatch(str, "+2"))
        strand = 2;

    else if (str_pmatch(str, "a") || str_pmatch(str, "d") || str_pmatch(str, "0"))  /* -auto/-dynamic */
        strand = 0;

    else
        fatal("unrecognized -strand option <%s>\n", option);

    return strand;
}

long extract_strand_option_no0(char *option)
{
    long strand = extract_strand_option(option);

    if (strand == 0)
        strand = UNSET_LVAL;

    return strand;
}

void check_topo_option(char *topo_list, char *topo)
{
    char *dft_topo_file = "up_to_octamers";

    if (is_empty_string(topo_list) && is_empty_string(topo))  /* default */
        get_sysfile(topo_list, "data/topology", dft_topo_file);

    if (!is_empty_string(topo_list) && !is_empty_string(topo))
        fatal("Specify either -topo or -topo_list, but NOT both!\n");

    if (!is_empty_string(topo_list) && !exist_file(topo_list))
        fatal("-topo_list=%s file does NOT exist!\n", topo_list);
}

void set_my_globals(char *pgname)
{
    Gvars.RUNLOG = stderr;
    Gvars.PRGLOG = NULL;
    Gvars.VERBOSITY = NORMAL;
    Gvars.RNA = FALSE;
    Gvars.WOBBLE = FALSE;
    Gvars.RAWID = FALSE;
    Gvars.RAWSEQ = FALSE;
    Gvars.PALINDROME = FALSE;
    Gvars.PROGNAME = pgname;

    Gvars.misc.MIN_COUNTS = 5;
    Gvars.misc.ALLPOS = FALSE;
    Gvars.misc.FIT4 = FALSE;
    Gvars.misc.POS[BUF510] = 0;  /* BUF512 hold Tcount */
    Gvars.misc.IUPAC_CUTOFF = 0.90;
    Gvars.misc.SSY = 1.0;

    Gvars.seed_criterion = SEED_BY_SLOPE_TVAL;
}

void cleanup_my_globals()
{
    close_file(Gvars.RUNLOG);
    close_file(Gvars.PRGLOG);
}

long check_common_options(char *option, char *helpfile, char *runlog, char *outdir)
{
    char str[BUF512];
    long k;

    strcpy(str, option);
    lowerstr(str);

    if (is_equal_string(str, "-h") || is_equal_string(str, "--help"))
        display_helpfile(helpfile);

    else if (is_equal_string(str, "-v"))
        show_version();

    else if (str_pmatch(str, "-verb")) {
        k = extract_verbose_option(str);
        if (k < 0)
            Gvars.VERBOSITY = SILENT;
        else if (k == 0)
            Gvars.VERBOSITY = NORMAL;
        else if (k < 10)
            Gvars.VERBOSITY = VERBOSE;
        else
            Gvars.VERBOSITY = VDEBUG;

    } else if (str_pmatch(str, "-silent") || str_pmatch(str, "-quiet"))
        Gvars.VERBOSITY = SILENT;

    else if (str_pmatch(str, "-rna"))
        Gvars.RNA = set_switch_with_dft_true(str);

    else if (str_pmatch(str, "-wobble") || is_equal_string(str, "-gu"))
        Gvars.WOBBLE = set_switch_with_dft_true(str);

    else if (str_pmatch(str, "-rawid") || str_pmatch(str, "-raw_id"))
        Gvars.RAWID = set_switch_with_dft_true(str);

    else if (str_pmatch(str, "-rawseq") || str_pmatch(str, "-raw_seq"))
        Gvars.RAWSEQ = set_switch_with_dft_true(str);

    else if (str_pmatch(str, "-palin"))
        Gvars.PALINDROME = set_switch_with_dft_true(str);

    else if (str_pmatch(str, "-run"))
        get_strvalue(option, runlog, TRUE);

    else if (str_pmatch(str, "-o=") || str_pmatch(str, "-out"))
        get_stropt_wo_end_slash(option, outdir);

    else if (str_pmatch(str, "-min_count"))
        Gvars.misc.MIN_COUNTS = get_lvalue(str, 1, BUFBIG);

    else if (str_pmatch(str, "-allpos"))
        Gvars.misc.ALLPOS = TRUE;

    else if (str_pmatch(str, "-fit4"))
        Gvars.misc.FIT4 = TRUE;

    else if (str_pmatch(str, "-iupac_cutoff"))
        Gvars.misc.IUPAC_CUTOFF = get_dvalue(str, 0.5, 1.0);

    else if (str_pmatch(str, "-ssy"))
        Gvars.misc.SSY = get_dvalue(str, 0.1, 1.0e6);

    else if (str_pmatch(str, "-seed_criterion")) {
        if (strstr(str, "slope") && !strstr(str, "val")) {
            Gvars.seed_criterion = SEED_BY_SLOPE_ITSELF;
            fprintf(Gvars.RUNLOG, "Find seeds based on |slope| itself\n");
        } else
            fprintf(Gvars.RUNLOG, "Find seeds by the t-value of slope [default]\n");
    } else
        return FALSE;

    return TRUE;
}

void display_unrecognized_option(char *option)
{
    char str[BUF512], tmp[BUF512], *ptr;

    sprintf(str, "***** unrecognized option: <%s> *****\n", option);

    ptr = strchr(option, '=');
    if (ptr == NULL) {
        sprintf(tmp, "    Did you forget to add '=val'? (i.e., `%s=val')\n", option);
        strcat(str, tmp);
    }

    quit_with_help_reminder(str);
}

void set_logs_and_check_option(char *runlog, char *outdir, long argc, long end_opt_idx,
                               char *option)
{
    set_logs(runlog, outdir);
    check_option_format(argc, end_opt_idx, option);
}

void set_logs(char *runlog, char *outdir)
{
    char str[BUF512], logfile[BUF512];

    bname_ext(Gvars.PROGNAME, "log", str);
    sprintf(logfile, "%s/%s", outdir, str);

    Gvars.PRGLOG = open_file(logfile, "w");

    if (Gvars.VERBOSITY == SILENT)  /* write to a temporary file */
        Gvars.RUNLOG = open_tmpfile();
    else
        Gvars.RUNLOG = open_file(runlog, "w");
}

void check_option_format(long argc, long end_opt_idx, char *option)
{
    char str[BUF512];

    if (argc > end_opt_idx) {
        sprintf(str, "options must be set in -par[=val] format [%s]\n", option);
        quit_with_help_reminder(str);
    }
}

void check_either_A_or_B(char *optA, char *optB, char *msg)
{
    char str[BUF512];

    if (is_empty_string(optA) && is_empty_string(optB)) {
        sprintf(str, "\tNone of %s specified -- pick up one only\n", msg);
        quit_with_help_reminder(str);

    } else if (!is_empty_string(optA) && !is_empty_string(optB)) {
        sprintf(str, "\tBoth of %s specified -- pick up one only\n", msg);
        quit_with_help_reminder(str);
    }
}

void check_required_option(char *option, char *invalid_str, char *msg)
{
    char str[BUF512];

    if (is_equal_string(option, invalid_str)) {
        sprintf(str, "missing required option: %s\n", msg);
        quit_with_help_reminder(str);
    }
}

void check_required_file(char *filename, char *invalid_str, char *msg)
{
    check_required_option(filename, invalid_str, msg);

    if (!exist_file(filename))
        fatal("required file <%s> does NOT exist!\n", filename);
}

long set_switch_with_dft_true(char *option)
{
    char str[BUF512];

    if (!strchr(option, '='))  /* just a switch */
        return TRUE;

    get_strvalue(option, str, FALSE);
    lowerstr(str);

    if (str_pmatch(str, "of") || str_pmatch(str, "0") ||  /* off | 0 | no | false */
        str_pmatch(str, "n") || str_pmatch(str, "f"))
        return FALSE;

    return TRUE;
}

long set_switch_with_dft_false(char *option)
{
    char str[BUF512];

    if (!strchr(option, '='))  /* just a switch */
        return FALSE;

    get_strvalue(option, str, FALSE);
    lowerstr(str);

    if (str_pmatch(str, "on") || str_pmatch(str, "1") ||  /* on | 1 | yes | true */
        str_pmatch(str, "y") || str_pmatch(str, "t"))
        return TRUE;

    return FALSE;
}
