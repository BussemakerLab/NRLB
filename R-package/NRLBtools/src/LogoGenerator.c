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

int main(int argc, char **argv)
{
    args_logo_generator args;  /* command-line options */

    set_my_globals(argv[0]);
    logo_cmdline(argc, argv, &args);

    if (is_equal_string(args.type, "psam"))
        logo_psam(&args);

    else if (is_equal_string(args.type, "pwm"))
        logo_pwm(&args);

    else
        logo_seq(&args);

    get_image(&args);

    cleanup_my_globals();

    return 0;
}

void logo_psam(args_logo_generator * args)
{
    char str[BUF512], msg[BUF512];
    long i, j, k = 0, num_pos;
    double dmin, dsum;

    struct_psam psam_data;
    struct_logo_info *logoData, *pld;

    read_psam_struct_lite(args->file, &psam_data);

    num_pos = psam_data.psam_length;
    dmin = min_dvector(psam_data.psam, 0, num_pos * NUM_BASE4 - 1);

    if ((dmin < 0.0) && !is_equal_string(args->style, "raw")) {
        sprintf(msg, "\tMIN=%g < 0 for -style=%s", dmin, args->style);
        log_msg(msg);
    }

    if (args->rc) {  /* get reverse complementary */
        rc_psam(&psam_data);
        sprintf(str, "%s/%s", args->outdir, "reverse_comp.xml");
        write_psam_struct(str, &psam_data);
    }

    logoData = allocate_memory_for_logoData(num_pos);

    for (i = 1; i <= num_pos; i++) {
        pld = &logoData[i];  /* to simplify later on reference */
        init_struct_logo_info(pld);

        dsum = 0.0;
        for (j = 0; j < NUM_BASE4; j++) {
            pld->raw_value[j] = psam_data.psam[k + j];
            dsum += pld->raw_value[j];
        }
        k += NUM_BASE4;

        pld->num_valid = dsum;

        get_frequency(pld);  /* deduce frequency data */
        set_heightByStyle(args, pld);
        sort_height(pld);
    }

    free_psam_data(&psam_data);

    generate_eps(num_pos, logoData, 0, NULL, args);

    free(logoData);
}

void logo_pwm(args_logo_generator * args)
{
    double w[BUFBIG];
    long num_nt;

    struct_psam psam_data;

    num_nt = read_pwm(args->file, w);  /* initial PWM file name */
    populate_pd_struct_from_psam(num_nt, w, 1, 1.0, &psam_data);

    sprintf(args->file, "%s/%s", args->outdir, "_pwm2psam.xml");
    write_psam_struct_lite(args->file, &psam_data);

    free_psam_data(&psam_data);

    logo_psam(args);
}

void logo_seq(args_logo_generator * args)
{
    long num_seqs, num_pos;

    struct_seq *sequences;  /* base sequence info */
    struct_logo_info *logoData;

    if (is_equal_string(args->type, "fasta"))
        sequences = read_logo_fasta(args->file, &num_seqs, args->outdir);
    else  /* flat sequence */
        sequences = read_logo_flat(args->file, &num_seqs, args->outdir);

    num_pos = sequences[1].nb;  /* length of logo sequence */
    logoData = allocate_memory_for_logoData(num_pos);

    get_logo_data(num_pos, logoData, num_seqs, sequences, args);
    generate_eps(num_pos, logoData, num_seqs, sequences, args);

    free_logo_sequences(num_seqs, sequences);
    free(logoData);
}

void get_image(args_logo_generator * args)
{
    char str[BUF512], src[BUF512], dst[BUF512], cmd[BUFBIG];
    char gs[BUF512] = "", gs_opts[BUF512] = "", convert[BUF512] = "";
    long width, height;
    FILE *fp;

    sprintf(src, "%s/%s", args->outdir, LOGO_PS_TMP);
    sprintf(dst, "%s/%s", args->outdir, args->logo);

    if (is_equal_string(args->format, "eps")) {
        copy_file(src, dst);
        return;
    }

    get_cfg(gs, gs_opts, convert, args->format);

    width = lround(args->width * PPCM);
    height = lround(args->height * PPCM);

    if (is_equal_string(args->format, "pdf"))
        sprintf(cmd, "%s -sOutputFile=%s \\\n"
                "\t-sDEVICE=pdfwrite \\\n"
                "\t-dPDFSETTINGS=/printer \\\n"
                "\t-dEmbedAllFonts=true \\\n"
                "\t-dDEVICEWIDTHPOINTS=%ld \\\n"
                "\t-dDEVICEHEIGHTPOINTS=%ld \\\n"
                "\t%s \\\n\t%s\n", gs, dst, width, height, gs_opts, src);

    else if (is_equal_string(args->format, "jpeg"))
        sprintf(cmd, "%s -sOutputFile=%s \\\n"
                "\t-sDEVICE=jpeg \\\n"
                "\t-dDEVICEWIDTHPOINTS=%ld \\\n"
                "\t-dDEVICEHEIGHTPOINTS=%ld \\\n"
                "\t%s \\\n\t%s\n", gs, dst, width, height, gs_opts, src);

    else if (is_equal_string(args->format, "png"))
        sprintf(cmd, "%s -sOutputFile=%s \\\n"
                "\t-sDEVICE=png16m \\\n"
                "\t-dDEVICEWIDTHPOINTS=%ld \\\n"
                "\t-dDEVICEHEIGHTPOINTS=%ld \\\n"
                "\t%s \\\n\t%s\n", gs, dst, width, height, gs_opts, src);

    else  /* gif */
        sprintf(cmd, "%s -sOutputFile=- \\\n"
                "\t-sDEVICE=png16m \\\n"
                "\t-dDEVICEWIDTHPOINTS=%ld \\\n"
                "\t-dDEVICEHEIGHTPOINTS=%ld \\\n"
                "\t%s \\\n"
                "\t%s \\\n"
                "\t| %s png:- %s\n", gs, width, height, gs_opts, src, convert, dst);
    call_system(cmd);

    sprintf(str, "%s/%s", args->outdir, LOGO_OPTIONS);
    fp = open_file(str, "a");  /* note mode 'a' */
    fprintf(fp, "\n");
    fprintf(fp, "# covert from %s [eps] to %s [%s]\n\n", src, dst, args->format);
    fprintf(fp, "%s", cmd);
    close_file(fp);
}

void get_cfg(char *gs, char *gs_opts, char *convert, char *fmt)
{
    char msg[BUF512], parfile[BUF512];
    FILE *fp;

    get_sysfile(parfile, "config", PKG_CFG);

    fp = open_file(parfile, "r");

    extract_xml_line_string(fp, "gs", gs, TRUE);
    sprintf(msg, "gs: %s", gs);
    log_msg_vchk(msg);

    extract_xml_line_string(fp, "gs_opts", gs_opts, TRUE);
    sprintf(msg, "gs_opts: %s", gs_opts);
    log_msg_vchk(msg);

    if (is_equal_string(fmt, "gif")) {  /* 'convert' is needed only for gif images */
        extract_xml_line_string(fp, "convert", convert, TRUE);
        sprintf(msg, "convert: %s", convert);
        log_msg_vchk(msg);
    }

    close_file(fp);
}

void set_heightByStyle(args_logo_generator * args, struct_logo_info * pld)
{
    long i;

    if (is_equal_string(args->style, "raw"))
        for (i = 0; i < NUM_BASE4; i++)
            pld->height[i] = pld->raw_value[i];

    else if (is_equal_string(args->style, "frequency"))
        for (i = 0; i < NUM_BASE4; i++)
            pld->height[i] = pld->frequency[i];

    else if (is_equal_string(args->style, "bits_info")) {
        get_ssc(args->smallSampleCorrection, pld);
        get_H(pld);
        get_R(args->stretch, pld);
        get_height(pld);

    } else {  /* ddG */
        pld->num_valid = 0.0;
        get_ddG_height(pld, args);
    }
}

void get_y_range(args_logo_generator * args, long num_pos, struct_logo_info * logoData)
{
    char msg[BUF512];
    long i, j, ymin, ymax;
    double dneg_sum, dpos, dneg, dval;

    dpos = get_maxSum(num_pos, logoData);
    ymax = dbl2long(dpos);  /* rounded to the nearest integer for tick */

    dneg = get_minSum(num_pos, logoData);
    ymin = dbl2long(dneg);

    if ((fabs(args->ymax) < XBIG_CUTOFF) && args->ymax >= dpos)
        ymax = args->ymax;

    if ((fabs(args->ymin) < XBIG_CUTOFF) && args->ymin <= dneg)
        ymin = args->ymin;

    sprintf(msg, "yaxis [%g (%ld); %g (%ld)]", dpos, ymax, dneg, ymin);
    log_msg_vchk(msg);

    for (i = 1; i <= num_pos; i++) {
        dneg_sum = 0.0;  /* sum of negative components at each position */

        for (j = 0; j < NUM_BASE4; j++) {
            dval = logoData[i].height[j];
            if (dval < 0.0)
                dneg_sum += dval;
        }

        logoData[i].yoffset = dneg_sum - ymin;  /* >= 0.0 */
        logoData[i].ymin = ymin;
        logoData[i].ymax = ymax;
    }
}

double get_maxSum(long num_pos, struct_logo_info * logoData)
{
    long i, j;
    double dval, dsum, dmax = -XBIG;

    for (i = 1; i <= num_pos; i++) {
        dsum = 0.0;

        for (j = 0; j < NUM_BASE4; j++) {
            dval = logoData[i].height[j];
            if (dval > 0.0)
                dsum += dval;
        }

        dmax = dval_max(dmax, dsum);
    }

    return dmax;
}

double get_minSum(long num_pos, struct_logo_info * logoData)
{
    long i, j;
    double dval, dsum, dmin = XBIG;

    for (i = 1; i <= num_pos; i++) {
        dsum = 0.0;

        for (j = 0; j < NUM_BASE4; j++) {
            dval = logoData[i].height[j];
            if (dval < 0.0)
                dsum += dval;
        }

        dmin = dval_min(dmin, dsum);  /* dsum <= 0.0 */
    }

    return dmin;
}

void get_ddG_height(struct_logo_info * pld, args_logo_generator * args)
{
    long i;
    double dval, min_Ka = args->min_Ka, log_geoAve = 0.0;

    for (i = 0; i < NUM_BASE4; i++) {  /* using raw PSAM data for ddG */
        if (is_equal_string(args->type, "psam"))
            dval = pld->raw_value[i];

        else  /* using frequency data for ddG in other cases */
            dval = pld->frequency[i];

        pld->frequency[i] = (dval < min_Ka) ? min_Ka : dval;
        log_geoAve += log(pld->frequency[i]);
    }

    log_geoAve /= 4.0;

    for (i = 0; i < NUM_BASE4; i++)
        pld->height[i] = log(pld->frequency[i]) - log_geoAve;
}

void write_matrix_info(FILE * fp, char *consensus_seq, long num_pos, char *style,
                       struct_logo_info * logoData, char *fmt, char *fmt0)
{
    long i, j;

    fprintf(fp, "\noriginal weight matrix:\n");
    for (j = 0; j < NUM_BASE4; j++)
        fprintf(fp, "\t %c", BASES[j]);
    fprintf(fp, "\n");

    print_sep(fp, DASH, 39);
    for (i = 1; i <= num_pos; i++) {
        fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
        for (j = 0; j < NUM_BASE4; j++)
            print_dval(fp, logoData[i].raw_value[j], fmt, fmt0);
        fprintf(fp, "\n");
    }

    if (is_equal_string(style, "ddG")) {
        fprintf(fp, "\nweight matrix after applying min_Ka cut-off:\n");
        for (j = 0; j < NUM_BASE4; j++)
            fprintf(fp, "\t %c", BASES[j]);
        fprintf(fp, "\n");

        print_sep(fp, DASH, 39);
        for (i = 1; i <= num_pos; i++) {
            fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
            for (j = 0; j < NUM_BASE4; j++)
                print_dval(fp, logoData[i].frequency[j], fmt, fmt0);
            fprintf(fp, "\n");
        }

    } else {
        fprintf(fp, "\nbase frequency at each position:\n");
        fprintf(fp, "\t dsum");
        for (j = 0; j < NUM_BASE4; j++)
            fprintf(fp, "\t %c", BASES[j]);
        fprintf(fp, "\n");

        print_sep(fp, DASH, 47);
        for (i = 1; i <= num_pos; i++) {
            fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
            print_dval(fp, logoData[i].num_valid, fmt, fmt0);

            for (j = 0; j < NUM_BASE4; j++)
                print_dval(fp, logoData[i].frequency[j], fmt, fmt0);
            fprintf(fp, "\n");
        }
    }
}

void write_seq_info(FILE * fp, char *consensus_seq, long num_pos, long num_seqs,
                    struct_seq * sequences, struct_logo_info * logoData, char *fmt,
                    char *fmt0)
{
    long i, j;

    fprintf(fp, "%ld sequence information sorted by ID:\n", num_seqs);
    for (i = 1; i <= num_seqs; i++)
        fprintf(fp, "%4ld %s %s\n", i, sequences[i].seq, sequences[i].id);

    fprintf(fp, "\nbase count at each position:\n");
    for (j = 0; j < NUM_BASE4; j++)
        fprintf(fp, "\t%c", BASES[j]);
    fprintf(fp, "\n");

    print_sep(fp, DASH, 34);
    for (i = 1; i <= num_pos; i++) {
        fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
        for (j = 0; j < NUM_BASE4; j++)
            print_dval(fp, logoData[i].raw_value[j], "\t%-7.0f", "\t%-7.0f");
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nbase frequency at each position:\n");
    for (j = 0; j < NUM_BASE4; j++)
        fprintf(fp, "\t %c", BASES[j]);
    fprintf(fp, "\n");

    print_sep(fp, DASH, 39);
    for (i = 1; i <= num_pos; i++) {
        fprintf(fp, "%2ld %c %2.0f", i, consensus_seq[i - 1], logoData[i].num_valid);
        for (j = 0; j < NUM_BASE4; j++)
            print_dval(fp, logoData[i].frequency[j], fmt, fmt0);
        fprintf(fp, "\n");
    }
}

void verify_logoData(char *consensus_seq, long num_pos, struct_logo_info * logoData,
                     long num_seqs, struct_seq * sequences, args_logo_generator * args)
{
    char filename[BUF512], *fmt = "\t%-7.4f", *fmt0 = "\t%-7.0f";
    long i, j;
    FILE *fp;

    sprintf(filename, "%s/%s", args->outdir, LOGO_DATA);
    fp = open_file(filename, "w");

    if (is_equal_string(args->type, "psam") || is_equal_string(args->type, "pwm"))
        write_matrix_info(fp, consensus_seq, num_pos, args->style, logoData, fmt, fmt0);
    else  /* based on sequence */
        write_seq_info(fp, consensus_seq, num_pos, num_seqs, sequences, logoData,
                       fmt, fmt0);

    fprintf(fp, "\nbase height at each position:\n");
    fprintf(fp, "\t e\t H\t R");
    for (j = 0; j < NUM_BASE4; j++)
        fprintf(fp, "\t %c", BASES[j]);
    fprintf(fp, "\n");

    print_sep(fp, DASH, 63);
    for (i = 1; i <= num_pos; i++) {
        fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
        print_dval(fp, logoData[i].e, fmt, fmt0);
        print_dval(fp, logoData[i].H, fmt, fmt0);
        print_dval(fp, logoData[i].R, fmt, fmt0);
        for (j = 0; j < NUM_BASE4; j++)
            print_dval(fp, logoData[i].height[j], fmt, fmt0);
        fprintf(fp, "\n");
    }

    fprintf(fp, "\nbase sorted height at each position:\n");
    for (i = 1; i <= num_pos; i++) {
        fprintf(fp, "%2ld %c", i, consensus_seq[i - 1]);
        for (j = 0; j < NUM_BASE4; j++) {
            fprintf(fp, "\t%c:", logoData[i].height_sorted[j].base);
            print_dval(fp, logoData[i].height_sorted[j].height, fmt, fmt0);
        }
        fprintf(fp, "\n");
    }

    close_file(fp);
}

void print_dval(FILE * fp, double d, char *fmt, char *fmt0)
{
    if (fabs(d) > 1.0e-6)
        fprintf(fp, fmt, d);
    else
        fprintf(fp, fmt0, d);
}

/* ------------------------- get frequency/height etc ------------------------- */
void generate_eps(long num_pos, struct_logo_info * logoData, long num_seqs,
                  struct_seq * sequences, args_logo_generator * args)
{
    char *consensus_seq, epsfile[BUF512];
    long yscale, y_overall_offset;
    FILE *fp;

    get_y_range(args, num_pos, logoData);

    yscale = logoData[1].ymax - logoData[1].ymin;
    y_overall_offset = -logoData[1].ymin;  /* to make it >= 0.0 */

    if (args->sb == UNSET_LVAL)
        args->sb = args->start_num;
    if (args->se == UNSET_LVAL)
        args->se = num_pos + args->start_num - 1;

    consensus_seq = cvector(0, num_pos);
    get_consensus_seq(num_pos, logoData, consensus_seq);

    sprintf(epsfile, "%s/%s", args->outdir, LOGO_PS_TMP);
    fp = open_file(epsfile, "w");

    get_ps_header(fp, args);
    get_ps_settings(fp, args, consensus_seq, yscale, y_overall_offset);
    get_logo_height(fp, args, num_pos, logoData);

    close_file(fp);

    verify_logoData(consensus_seq, num_pos, logoData, num_seqs, sequences, args);

    free_cvector(consensus_seq, 0, DUMMY);
}

void get_consensus_seq(long num_pos, struct_logo_info * logoData, char *consensus_seq)
{
    long i;

    for (i = 1; i <= num_pos; i++)
        consensus_seq[i - 1] = logoData[i].height_sorted[TOP_BIDX3].base;

    consensus_seq[num_pos] = '\0';
}

void get_logo_height(FILE * fp, args_logo_generator * args, long num_pos,
                     struct_logo_info * logoData)
{
    char base;
    double height;
    long i, j, k;

    struct_base_height *bh;

    for (i = 1; i <= num_pos; i++) {
        k = args->start_num + i - 1;  /* current base label number */
        if (k < args->sb || k > args->se)
            continue;  /* out of requested range */

        fprintf(fp, "(%ld) startStack\n", k);
        bh = logoData[i].height_sorted;

        if (logoData[i].yoffset > 0)
            fprintf(fp, "  %g (X) makeSymbol\n", logoData[i].yoffset);

        for (j = 0; j < NUM_BASE4; j++) {
            height = bh[j].height;
            base = bh[j].base;

            if (Gvars.RNA)
                base = cvt_base_T2U(base);

            if (height > 0) {
                fprintf(fp, "  %g (%c) makeSymbol\n", height, base);

            } else if (height < 0) {
                fprintf(fp, "  %g (%c) makeSymbol\n", -height, tolower((int) base));

            } else  /* skip the base with ZERO height */
                continue;
        }

        if (args->errorBar)
            fprintf(fp, "  %g setIbeam\n", logoData[i].e);
        fprintf(fp, "endStack\n\n");
    }

    fprintf(fp, "endLine\n");

    if (args->frame)
        fprintf(fp, "x0 0 contentWidth 0 contentWidth contentHeight"
                " x0 contentHeight drawBox\n");
    else
        fprintf(fp, "x0 0 contentWidth 0 drawLine\n");

    fprintf(fp, "endLogo\n");
}

void set_ytick(FILE * fp, long yscale, long y_overall_offset)
{
    if (y_overall_offset > 0)  /* with negative value */
        if (yscale <= 10.0)
            fprintf(fp, "/yaxisTicBits 1 def\n");
        else if (yscale <= 20.0)
            fprintf(fp, "/yaxisTicBits 2 def\n");
        else if (yscale <= 50.0)
            fprintf(fp, "/yaxisTicBits 5 def\n");
        else
            fprintf(fp, "/yaxisTicBits 10 def\n");

    else {
        if (yscale <= 2.0)
            fprintf(fp, "/yaxisTicBits 0.5 def\n");
        else if (yscale <= 5.0)
            fprintf(fp, "/yaxisTicBits 1 def\n");
        else if (yscale <= 10.0)
            fprintf(fp, "/yaxisTicBits 2 def\n");
        else if (yscale <= 20.0)
            fprintf(fp, "/yaxisTicBits 5 def\n");
        else
            fprintf(fp, "/yaxisTicBits 10 def\n");
    }
}

void set_margins(FILE * fp, long label)
{
    long uwidth = 1;

    fprintf(fp, "\n");

    if (label) {
        fprintf(fp, "/topMargin\n"
                "    logoTitle () eq {10} {titleFontsize 4 add} ifelse\n" "def\n");
        fprintf(fp, "/rightMargin  %% add extra room if showing ends\n"
                "    showEnds (-) eq {fontsize} {fontsize 1.5 mul} ifelse\n" "def\n");
        fprintf(fp, "/bottomMargin\n"
                "    fontsize 0.75 mul\n"
                "    xaxis {fontsize 1.75 mul add} if %% add extra room for axis\n"
                "    xaxisLabel () eq {} {fontsize 0.75 mul add} ifelse\n" "def\n");
        fprintf(fp, "/leftMargin\n" "    fontsize 3.6 mul\n" "def\n");

    } else {
        fprintf(fp, "/topMargin %ld def\n", uwidth);
        fprintf(fp, "/rightMargin %ld def\n", uwidth);
        fprintf(fp, "/bottomMargin %ld def\n", uwidth);
        fprintf(fp, "/leftMargin %ld def\n", uwidth);
    }
}

void get_ps_settings(FILE * fp, args_logo_generator * args, char *consensus_seq,
                     long yscale, long y_overall_offset)
{
    char box_style, *p0, *p, parfile[BUF512];
    FILE *fp_ps;

    get_sysfile(parfile, "html", LOGO_PS_PAR);

    fp_ps = open_file(parfile, "r");
    while ((p0 = my_getline(fp_ps)) != NULL) {
        p = rtrim(p0);  /* keep left indentation */

        if (case_strstr(p, "LogoGenerator_CSC") == NULL)
            fprintf(fp, "%s\n", p);

        else {
            fprintf(fp,
                    "%% prolog part #2: Case Specific Settings [LogoGenerator_CSC]\n\n");

            if (args->bw)  /* black-and-white */
                fprintf(fp, "/colorDict <<\n    (X) white\n>> def\n\n");

            fprintf(fp, "/logoWidth %g cm def\n", args->width);
            fprintf(fp, "/logoHeight %g cm def\n", args->height);
            fprintf(fp, "/logoTitle (%s) def\n\n", args->title);

            fprintf(fp, "/yaxis true def\n");
            if (is_equal_string(args->style, "raw"))
                fprintf(fp, "/yaxisLabel (raw data) def\n");

            else if (is_equal_string(args->style, "frequency"))
                fprintf(fp, "/yaxisLabel (frequency) def\n");

            else if (is_equal_string(args->style, "bits_info"))
                fprintf(fp, "/yaxisLabel (bits) def\n");

            else
                fprintf(fp, "/yaxisLabel (ddG/RT) def\n");

            if (args->yaxis == UNSET_LVAL)  /* default */
                fprintf(fp, "/yaxisBits %ld def\n", yscale);
            else
                fprintf(fp, "/yaxisBits %g def\n", args->yaxis);

            if (args->ytick == UNSET_LVAL)  /* default */
                set_ytick(fp, yscale, y_overall_offset);
            else
                fprintf(fp, "/yaxisTicBits %g def\n", args->ytick);

            fprintf(fp, "/yoffset %ld def\n\n", y_overall_offset);

            fprintf(fp, "/xaxis true def\n");
            fprintf(fp, "/xaxisLabel () def\n");
            fprintf(fp, "/showEnds (-) def\n\n");  /* d: DNA; p: Protein; -: none */

            if (args->consensus || !is_empty_string(args->subtitle))
                fprintf(fp, "/showFineprint true def\n");
            else
                fprintf(fp, "/showFineprint false def\n");

            if (!is_empty_string(args->subtitle))
                fprintf(fp, "/fineprint (%s) def\n\n", args->subtitle);
            else
                fprintf(fp, "/fineprint (%s) def\n\n", consensus_seq);

            fprintf(fp, "/charsPerLine %ld def\n", args->se - args->sb + 1);
            fprintf(fp, "/logoLines 1 def\n\n");

            /* "n" to have No boxes around characters;
             * "s" to have boxes around characters, with Shrinking;
             * "f" to have Filled boxes around characters */
            box_style = (args->box == 1) ? 'f' : (args->box == -1) ? 's' : 'n';
            fprintf(fp, "/showingBox (%c) def\n", box_style);
            fprintf(fp, "/shrinking true def\n");
            fprintf(fp, "/shrink 1 def\n");  /* 1 (no shrinking) to 0 (full shrinking) */
            fprintf(fp, "/outline %s def\n\n", args->outline ? "true" : "false");

            fprintf(fp, "/setIbeamFraction  1 def\n");
            fprintf(fp, "/setIbeamGray 0.5 def\n");
            fprintf(fp, "/setIbeamLineWidth 0.5 def\n");

            set_margins(fp, args->label);

            if (args->margin)
                fprintf(fp, "/stackMargin %ld def\n", args->margin);
            else  /* to avoid bug for 0-margin in early versions of gs/gv */
                fprintf(fp, "/stackMargin 0.01 def\n");
        }
        free(p0);
    }

    close_file(fp_ps);
}

void get_ps_header(FILE * fp, args_logo_generator * args)
{
    char ctString[BUF512];
    long llx = 0, lly = 0, urx, ury;  /* bounding box */

    urx = lround(args->width * PPCM);
    ury = lround(args->height * PPCM);

    fprintf(fp, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(fp, "%%%%BoundingBox: %ld %ld %ld %ld\n", llx, lly, urx, ury);
    fprintf(fp, "%%%%Title: (LogoGenerator: %s)\n", VERSION_NUMBER);
    fprintf(fp, "%%%%Creator: (Xiang-Jun Lu [xl2134@columbia.edu; Bussemaker Lab])\n");

    get_currentTimeString(ctString);
    fprintf(fp, "%%%%CreationDate: (%s)\n\n", ctString);
}

void get_logo_data(long num_pos, struct_logo_info * logoData, long num_seqs,
                   struct_seq * sequences, args_logo_generator * args)
{
    char cbase, *p;
    long i, j, num_valid;

    struct_logo_info *pld;  /* pointer to struct_logo_info struct */

    /* count the number of occurrence of each base at each position */
    for (i = 1; i <= num_pos; i++) {
        pld = &logoData[i];  /* to simplify later on reference: 1-index */
        init_struct_logo_info(pld);

        num_valid = 0;
        for (j = 1; j <= num_seqs; j++) {
            cbase = sequences[j].seq[i - 1];  /* 0-index! */
            p = strchr(BASES, toupper((int) cbase));  /* to make it case-insensitive */

            if (p != NULL) {
                pld->raw_value[p - BASES]++;
                num_valid++;
            }
        }

        pld->num_valid = (double) num_valid;

        if (!num_valid)  /* not a single valid base at current position */
            continue;

        get_frequency(pld);
        set_heightByStyle(args, pld);
        sort_height(pld);
    }
}

void get_ssc(long smallSampleCorrection, struct_logo_info * pld)
/* get the small sample correction term */
{
    double ssc_const = TOP_BIDX3 / (2.0 * log(2));

    if (!pld->num_valid)
        return;  /* pld->e = 0.0 by default */

    if (smallSampleCorrection)
        pld->e = ssc_const / pld->num_valid;
}

void get_frequency(struct_logo_info * pld)
{
    long i;

    if (!pld->num_valid)
        return;

    for (i = 0; i < NUM_BASE4; i++)
        pld->frequency[i] = pld->raw_value[i] / pld->num_valid;
}

void get_H(struct_logo_info * pld)
{
    long i;
    double freq, H = 0.0;

    for (i = 0; i < NUM_BASE4; i++) {
        freq = pld->frequency[i];
        if (freq > 0)
            H += freq * my_log2(freq);
    }

    pld->H = -H;
}

void get_R(long stretch, struct_logo_info * pld)
/* amount of information at position 'l' */
{
    double R = 2.0;  /* full height if stretched */

    if (!stretch)
        R -= (pld->H + pld->e);

    pld->R = R;
}

void get_height(struct_logo_info * pld)
/* get the plot height of the each base at position 'l': R * freq */
{
    long i;

    for (i = 0; i < NUM_BASE4; i++)
        pld->height[i] = pld->R * pld->frequency[i];
}

int height_compare(const void *v1, const void *v2)
/* by height; then base */
{
    const struct_base_height *p1, *p2;

    p1 = (const struct_base_height *) v1;
    p2 = (const struct_base_height *) v2;

    if (p1->height > p2->height + DBL_EPSILON)
        return 1;
    else if (p1->height < p2->height - DBL_EPSILON)
        return -1;
    else
        return p2->base - p1->base;  /* reverse order for A..C..G..T */
}

void sort_height(struct_logo_info * pld)
/* sort base according to its height at each position in the logo so
 * that the most significant one is at the top; with equal height,
 * bases are arranged in reverse order (A > C > G > T) */
{
    long i;

    /* firstly make a copy from array height[] */
    for (i = 0; i < NUM_BASE4; i++) {
        pld->height_sorted[i].base = toupper((int) BASES[i]);  /* to upper case */
        pld->height_sorted[i].height = pld->height[i];
    }

    qsort(pld->height_sorted, NUM_BASE4, sizeof(pld->height_sorted[1]), height_compare);
}

void init_struct_logo_info(struct_logo_info * pld)
{
    long i;

    pld->num_valid = 0.0;
    pld->e = 0.0;
    pld->H = 0.0;
    pld->R = 0.0;

    for (i = 0; i < NUM_BASE4; i++) {
        pld->raw_value[i] = 0.0;
        pld->frequency[i] = 0.0;
        pld->height[i] = 0.0;
        pld->height_sorted[i].base = ' ';
        pld->height_sorted[i].height = 0.0;
    }

    pld->yoffset = 0.0;
    pld->ymin = 0;
    pld->ymax = 0;
}

/* ------------------------- command-line processing ------------------------- */
void set_logo_defaults(args_logo_generator * args)
{
    strcpy(args->file, "");
    strcpy(args->format, "png");
    strcpy(args->type, "psam");
    strcpy(args->logo, "");
    strcpy(args->outdir, ".");
    strcpy(args->runlog, "stderr");
    strcpy(args->title, "");
    strcpy(args->subtitle, "");
    strcpy(args->style, "bits_info");

    args->smallSampleCorrection = FALSE;
    args->errorBar = FALSE;
    args->stretch = FALSE;
    args->min_Ka = WEPS;
    args->width = 12;
    args->height = 7.5;
    args->bw = FALSE;
    args->rc = FALSE;
    args->box = FALSE;  /* no rectangular box shown */
    args->frame = FALSE;  /* do not drawing a bounding box */
    args->consensus = FALSE;  /* do not show consensus sequence as fineprint */
    args->outline = FALSE;  /* without outline, but solid */
    args->start_num = 1;
    args->sb = UNSET_LVAL;
    args->se = UNSET_LVAL;  /* end of motif sequence length */
    args->margin = 1;
    args->yaxis = UNSET_LVAL;
    args->ytick = UNSET_LVAL;
    args->ymin = -XBIG;
    args->ymax = XBIG;
    args->label = TRUE;  /* have label */
}

static void write_logo_options(args_logo_generator * args)
{
    char filename[BUF512];
    FILE *fp;

    sprintf(filename, "%s/%s", args->outdir, LOGO_OPTIONS);
    fp = open_file(filename, "w");

    fprintf(fp, "## %s\n\n", getenv("PWD"));  /* current working directory */
    fprintf(fp, "%s \\\n", Gvars.PROGNAME);

    fprintf(fp, "\t-file=%s \\\n", args->file);
    fprintf(fp, "\t-type=%s \\\n", args->type);
    fprintf(fp, "\t-style=%s \\\n", args->style);

    if (is_equal_string(args->style, "ddG"))
        fprintf(fp, "\t-min_Ka=%g \\\n", args->min_Ka);

    if (is_equal_string(args->format, "eps"))
        fprintf(fp, "\t-logo=%s \\\n", args->logo);
    else
        fprintf(fp, "\t-logo=%s \\\n", LOGO_PS_TMP);

    fprintf(fp, "\t-width=%g \\\n", args->width);
    fprintf(fp, "\t-height=%g \\\n", args->height);

    fprintf(fp, "\t-start_num=%ld \\\n", args->start_num);
    if (args->sb != UNSET_LVAL)
        fprintf(fp, "\t-sb=%ld \\\n", args->sb);
    if (args->se != UNSET_LVAL)
        fprintf(fp, "\t-se=%ld \\\n", args->se);

    if (args->smallSampleCorrection)
        fprintf(fp, "\t-smallSampleCorrection \\\n");

    if (args->errorBar)
        fprintf(fp, "\t-errorBar \\\n");

    if (args->stretch)
        fprintf(fp, "\t-stretch \\\n");

    if (args->bw)
        fprintf(fp, "\t-bw \\\n");

    if (args->rc)
        fprintf(fp, "\t-reverse_complementary \\\n");

    if (args->outline)
        fprintf(fp, "\t-outline \\\n");

    if (args->frame)
        fprintf(fp, "\t-frame \\\n");

    if (args->consensus)
        fprintf(fp, "\t-consensus \\\n");
    else if (!is_empty_string(args->subtitle))
        fprintf(fp, "\t-subtitle=%s \\\n", args->subtitle);

    if (!args->label)
        fprintf(fp, "\t-label=off \\\n");

    if (args->box)
        fprintf(fp, "\t-box=%ld \\\n", args->box);

    if (!is_empty_string(args->title))
        fprintf(fp, "\t-title=\"%s\" \\\n", args->title);

    write_option_outdir_runlog(fp, args->outdir, args->runlog);

    fprintf(fp, "\t-format=eps\n");  /* always produce an EPS image */

    close_file(fp);
}

void check_logo_cmdline(args_logo_generator * args)
{
    char str[BUF512];

    check_required_file(args->file, "", "-file=file_name");

    lowerstr(args->type);
    if (str_pmatch(args->type, "fl"))
        strcpy(args->type, "flat");
    else if (str_pmatch(args->type, "fa"))
        strcpy(args->type, "fasta");
    else if (str_pmatch(args->type, "pw"))
        strcpy(args->type, "pwm");
    else
        strcpy(args->type, "psam");

    lowerstr(args->format);
    if (str_pmatch(args->format, "pd"))
        strcpy(args->format, "pdf");
    else if (str_pmatch(args->format, "ps") || str_pmatch(args->format, "e"))
        strcpy(args->format, "eps");
    else if (str_pmatch(args->format, "j"))
        strcpy(args->format, "jpeg");
    else if (str_pmatch(args->format, "g"))
        strcpy(args->format, "gif");
    else
        strcpy(args->format, "png");

    lowerstr(args->style);
    if (*args->style == 'r')
        strcpy(args->style, "raw");
    else if (*args->style == 'f')
        strcpy(args->style, "frequency");
    else if (*args->style == 'b')
        strcpy(args->style, "bits_info");
    else if (*args->style == 'd')
        strcpy(args->style, "ddG");
    else {
        if (is_equal_string(args->type, "psam"))
            strcpy(args->style, "ddG");
        else
            strcpy(args->style, "bits_info");
    }

    if (is_empty_string(args->logo))  /* not specified on command line */
        bname_ext(args->file, args->format, args->logo);

    else if (strrchr(args->logo, '.') == NULL) {  /* no extension */
        sprintf(str, "%s.%s", args->logo, args->format);
        strcpy(args->logo, str);
    }

    if (is_equal_string(args->type, "psam") || is_equal_string(args->type, "pwm")) {
        args->errorBar = FALSE;  /* no small sample error correction */
        args->smallSampleCorrection = FALSE;

    } else  /* -reverse_comp only applies to 'psam' */
        args->rc = FALSE;
}

void logo_cmdline(int argc, char *argv[], args_logo_generator * args)
{
    char helpfile[BUF512];
    long i;

    set_logo_defaults(args);

    bname_ext(Gvars.PROGNAME, "hlp", helpfile);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != DASH)
            break;

        if (check_common_options(argv[i], helpfile, args->runlog, args->outdir))
            continue;

        if (str_pmatch(argv[i], "-fi"))  /* input file name */
            get_strvalue(argv[i], args->file, TRUE);

        else if (str_pmatch(argv[i], "-sty"))  /* input data PSAM | PWM | fasta | flat */
            get_strvalue(argv[i], args->style, FALSE);

        else if (str_pmatch(argv[i], "-fo"))  /* image: EPS | PDF | JPEG | PNG | GIF */
            get_strvalue(argv[i], args->format, FALSE);

        else if (str_pmatch(argv[i], "-ty"))  /* input type: matrix | fasta | flat */
            get_strvalue(argv[i], args->type, FALSE);

        else if (str_pmatch(argv[i], "-lo"))  /* output (logo) image file name */
            get_strvalue(argv[i], args->logo, TRUE);

        else if (str_pmatch(argv[i], "-la"))  /* if to have label on */
            args->label = set_switch_with_dft_true(argv[i]);

        else if (str_pmatch(argv[i], "-ti"))  /* title line */
            get_strvalue(argv[i], args->title, FALSE);

        else if (str_pmatch(argv[i], "-su"))  /* subtitle line */
            get_strvalue(argv[i], args->subtitle, FALSE);

        else if (str_pmatch(argv[i], "-mi"))  /* minimum allowable Ka */
            args->min_Ka = get_dvalue(argv[i], WEPS, 1.0);

        else if (str_pmatch(argv[i], "-w"))  /* logo width in cm */
            args->width = get_dvalue(argv[i], 2, BUF512);

        else if (str_pmatch(argv[i], "-he"))  /* logo height in cm */
            args->height = get_dvalue(argv[i], 2, BUF512);

        else if (str_pmatch(argv[i], "-sm"))
            args->smallSampleCorrection = TRUE;

        else if (str_pmatch(argv[i], "-e"))
            args->errorBar = TRUE;

        else if (str_pmatch(argv[i], "-str"))
            args->stretch = TRUE;

        else if (str_pmatch(argv[i], "-sta"))
            args->start_num = get_lvalue(argv[i], -BUFBIG, BUFBIG);

        else if (str_pmatch(argv[i], "-sb"))
            args->sb = get_lvalue(argv[i], -BUFBIG, BUFBIG);

        else if (str_pmatch(argv[i], "-se"))
            args->se = get_lvalue(argv[i], -BUFBIG, BUFBIG);

        else if (str_pmatch(argv[i], "-bo"))
            args->box = get_lvalue(argv[i], -BUFBIG, BUFBIG);

        else if (str_pmatch(argv[i], "-bw"))
            args->bw = TRUE;

        else if (str_pmatch(argv[i], "-re") || str_pmatch(argv[i], "-rc"))
            args->rc = TRUE;

        else if (str_pmatch(argv[i], "-outl"))
            args->outline = TRUE;

        else if (str_pmatch(argv[i], "-fr"))
            args->frame = TRUE;

        else if (str_pmatch(argv[i], "-c"))
            args->consensus = TRUE;

        else if (str_pmatch(argv[i], "-ma"))
            args->margin = get_lvalue(argv[i], 0, 6);

        else if (str_pmatch(argv[i], "-ya"))
            args->yaxis = get_dvalue(argv[i], 0, BUF512);

        else if (str_pmatch(argv[i], "-yt"))
            args->ytick = get_dvalue(argv[i], 0, BUF512);

        else if (str_pmatch(argv[i], "-ymi"))
            args->ymin = get_dvalue(argv[i], -XBIG, XBIG);

        else if (str_pmatch(argv[i], "-yma"))
            args->ymax = get_dvalue(argv[i], -XBIG, XBIG);

        else
            display_unrecognized_option(argv[i]);
    }

    set_logs_and_check_option(args->runlog, args->outdir, argc, i, argv[i]);

    check_logo_cmdline(args);
    write_logo_options(args);
}
