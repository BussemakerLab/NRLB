#ifndef _REDUCE_SUITE_H
#define _REDUCE_SUITE_H

/* common header files */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>

#define NR_END 1  /* for NRC related functions */
#define FREE_ARG char*

#define UNUSED_PARAMETER(x) (void)(x)  /* to prevent warning about unused parameters */
#define VERSION_NUMBER "REDUCE Suite v2.2.3-2017jun26"

#define TRUE 1L
#define FALSE 0L
#define DUMMY -1L
enum { SILENT = -1, NORMAL = 0, VERBOSE = 1, VDEBUG = 2 };

#define UNSET_LVAL -1024
#define LUXDBG FALSE  /* set to TRUE for more info */

#define WMAX_CUTOFF 5  /* maximum allowed Wi at each position before rescaling */
#define WEPS 0.001  /* should not be too small, e.g., 1.0e-6 */
#define P_VALUE 0.001  /* default p-value */
#define NTOP 10  /* number of top seeds to print */

#define CMPI 2.54  /* cm per inch */
#define PPI 72  /* point per inch */
#define PPCM 28.346  /* point per cm: 72.0 / 2.54 */

#define QUALEN16 16L  /* fast decimal <----> quaternary conversion */
#define MAXNTS_TOP 15L  /* max. no. of unknown nt positions in topo pattern */
#define MAXNTS_DIC 10  /* max. length of ACGT motif in dictfile by topo-form */
#define DEBUG_NCOUNT 9  /* above which more debug info is available */

#define NUM_BASE4 4
#define TOP_BIDX3 3  /* top index in BASES: 0..3 */

#define IUPAC_PATTERN "ACGTUWRKYSMBHDVN-()"
#define IUPAC_D1 "ACGTU-)"
#define IUPAC_D2 "RYKMWS"
#define IUPAC_D3 "DHVB"
#define IUPAC_D4 "N("

#define IUPAC_BASES "ACGTUWRKYSMBHDVN"
#define IUPAC_NON "WRKYSMBHDV"
#define IUPAC_N "WRKYSMBHDVN"
#define BASES "ACGT"  /* converted to 0123 */
#define RC_BASES "TGCA"
#define CASE_BASES "ACGTacgt"
#define CASE_RC_BASES "TGCAtgca"

#define RNA_U 'U'
#define RNA_u 'u'
#define DNA_T 'T'
#define DNA_t 't'

#define BUF8       8L
#define BUF16      16L
#define BUF32      32L
#define BUF510     510L
#define BUF512     512L
#define BUF1024    1024L
#define BUFBIG     8192L
#define SFACTOR    2L  /* scale factor for increasing memory */

#define ITMAX     6000L
#define NR_EPS    DBL_EPSILON
#define NR_FPMIN  DBL_MIN
#define TOLX      1.0e-10
#define XBIG 1.0e+18
#define XBIG_CUTOFF 1.0e+16

#define TABCHAR '\t'  /* tab-delimited fields */
#define WSPACES " \t\n"  /* space, tab & newline */
#define LIST_SEP ",;"  /* list items separator */
#define SKIPS "#\0"  /* lines to be skipped */
#define FASTA_SKIPS ";#\0"  /* skipped lines in FASTA */

#define PKG_CFG "pkg_settings.cfg"

enum {
    SEED_BY_SLOPE_ITSELF = 0,
    SEED_BY_SLOPE_TVAL
};

/* global variables defined here */
struct {
    FILE *RUNLOG;  /* channel for printing diagnostic message */
    FILE *PRGLOG;  /* log of program run */
    long VERBOSITY;  /* switch on/off diagnostic messages to screen */
    long RNA;  /* switch on/off for ACGU vs ACGT alphabets */
    long WOBBLE;  /* switch to count G-U pair in stem region */
    long RAWID;  /* switch to keep id AS IS */
    long RAWSEQ;  /* switch to keep raw sequence */
    long PALINDROME;  /* switch to check for palindrome */
    char *PROGNAME;
    struct {
        long MIN_COUNTS;  /* minimum motif counts in MatrixREDUCE/MotifREDUCE */
        long ALLPOS;  /* switch to optimize all positions, regardless of topo */
        long FIT4;  /* switch to fit 4 instead of 3 Ws at a position */
        long POS[BUF512];  /* positions to optimize in OptimizePSAM, also set FIT4 */
        double IUPAC_CUTOFF;  /* cut off to decide a IUPAC code for a position */
        double SSY;  /* minimum total sum of squares in ls-fitting */
    } misc;
    long seed_criterion;
} Gvars;  /* group global variables in struct for easy handling */

/* ======================================================================== */
typedef struct {
    long idx;  /* serial number matching id */
    char *idstr;  /* string id */
} struct_id;

/* ======================================================================== */
typedef struct {
    char *str;
    char *tag;
} struct_tag;

/* ======================================================================== */
typedef struct {
    long num;  /* total number of motifs */
    long *index;  /* decimal representation of normal ACGT motifs; -1 otherwise */

    long num_nonbit;  /* cases such as non-ACGT IUPAC motifs, or length > MAXNTS_DIC */
    long *idx_nonbit;

    long num_seidx;  /* number of different length in motifs */
    long **seidx;  /* starting ... ending index of motifs per length */

    struct_tag *motifs;
} struct_dic;

/* ======================================================================== */
typedef struct {  /* tab-delimited data file (measurement file) */
    long nrow;  /* number of rows */
    long ncol;  /* number of columns */
    long *seqidx;  /* pointer to corresponding sequence */
    long *oknum;  /* valid number in each column */
    char **ids;  /* IDs -- from the first column */
    char **col_names;  /* header column */
    double **data;  /* actual data matrix */
    double **resid;  /* the residual data matrix; could be NULL */
} struct_data;  /* tab-delimited text data */

/* ======================================================================== */
#define REPEATED_ID  "__REPEAT__"
#define BASE_PER_LINE 60

/* ======================================================================== */
typedef struct {
    char *title;  /* original title line */
    char *id;  /* converted sequence ID */
    long nb;  /* length of the sequence */
    char *seq;  /* actual sequence */
    long num_fragments;  /* total no. of segments with only 'ACGT' */
    long *num_nb;  /* no. of bases in each segment */
    char **seq_fragments;  /* pointer to each segments starting position in seq */
} struct_seq;

typedef struct {
    long num_seqs;
    struct_seq *sequences;
} struct_seqArr;

/* ======================================================================== */
#define XCHR 'X'
#define NCHR 'N'
#define DASH '-'

typedef struct {
    struct_tag str_tag;  /* e.g., XXX-XXXX or X3-X4 */
    long Tcount;  /* total # of characters: 8 for above */
    long Xcount;  /* 7 for above */
    long *offset;  /* 0, 1, 2, 4, 5, 6, 7: size = 7, no position 3 */
    long msize;  /* 4^Xcount */

    long num_parens;  /* number of parenthesis pairs () */
    long parens_bidx[BUF16];  /* 16 should be more than enough */
    long parens_eidx[BUF16];  /* matching offset */
} struct_topo;

typedef struct {
    long num_topo;
    struct_topo *topos;
} struct_topoArr;

/* ======================================================================== */
typedef struct {
    long nrow;  /* total no. of entries of measurement value */
    long num;  /* valid no., i.e., with matched seq and non-NaN */
    long strand;
    double *yval;  /* measurements */
    struct_seqArr *seqs;
    struct_topo *topo;
} struct_lm;  /* data structure for Lev-Mar optimization */

/* ======================================================================== */
typedef struct {
    long num;
    long ma;
    long bonf;
    double chisq;

    double *fval;
    double *tval;
    double *pval;
} struct_fitpars;

/* ======================================================================== */
typedef struct {
    long num;  /* number of valid pairs */
    double sumX;  /* sum of the x variable: no. of motif counts */
    double SSY;  /* initial variance */
    double SSE;  /* residual variance */
    double SSR;  /* variance accounted for */
    double r2;  /* r squared */

    double F;  /* F value, sign follows b_tval */
    double a;  /* intercept */
    double a_tval;  /* its t-value */
    double a_pval;  /* its p-value */
    double b;  /* slope */
    double b_tval;  /* its t-value */
    double b_pval;  /* its p-value */
} struct_ls5fit;  /* linear least-squares fit based on five sums */

/* ======================================================================== */
typedef struct {
    struct_ls5fit fitpars;  /* get all the linear fit parameters */

    double *x;  /* 1-indexed: authentic LM optimized parameters, including C & F */
    double *w;  /* 0-indexed: column-normalized x^2 (4 * Xcount), without C and F */
    double *psam;  /* 0-indexed: full size PSAM (4 * Tcount) */

    long num_mids;  /* number of IDs with matched motif in its sequence */
    long bonf;  /* for bonferroni correction */
    long col_idx;  /* experimental column in the measurement file */
    char *col_name;  /* name of the experiment */

    long topo_idx;  /* topology pattern index */
    struct_topo *topo;  /* pointer to a pattern */

    long index;  /* motif index in decimal */
    char *motif;  /* compact form, w/o gaps */
    char *full;  /* with gaps */
    char *optimal;  /* optimal sequence */

    long strand;
    char *stnd_msg;
    double p_value;
} struct_seed;

typedef struct {
    long col_idx;  /* column index of the measurement file */
    long topo_idx;  /* topo index */
    long index;  /* decimal for topo, or link to struct_dic arrays */
    long strand;
    double b_tval;
} struct_ntop;

typedef struct {
    long ntop;
    long nsel;
    long bonf;

    struct_ntop *ntops;
} struct_ntopArr;

/* ======================================================================== */
#define MAX_MOTIF 20L

#define MATRIXREDUCE_RAW "matrixreduce_results.raw"
#define MATRIXREDUCE_XML "matrixreduce_results.xml"
#define MATRIXREDUCE_HTML "matrixreduce_results.html"

#define MOTIFREDUCE_RAW "motifreduce_results.raw"
#define MOTIFREDUCE_XML "motifreduce_results.xml"
#define MOTIFREDUCE_HTML "motifreduce_results.html"

#define NUM_PROBE_NO_SEQ "number_of_probes_missing_sequence="
#define NUM_PROBE_REPEAT "number_of_repeated_probes="

#define RESIDUAL_FILE "model_residuals.tsv"
#define PREDICT_FILE "model_predicted.tsv"

typedef struct {
    char seqfile[BUF512];  /* sequence file */
    char measfile[BUF512];  /* name of measurements */
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];  /* running log */

    char topo[BUF512];  /* single topology pattern */
    char topo_list[BUF512];  /* name of topology file */
    char dictfile[BUF512];  /* motifs list dictionary  */

    double p_value;  /* p-value cut-off for including a new motif */

    long max_motif;  /* maximum number of motifs */
    long strand;  /* +1 for leading; -1 for rc; 2 for both; 0 dynamically */

    long iupac_pos;
    char iupac_sym[BUF512];

    long ntop;  /* number of top seeds to print out for information */
    long nsel;  /* selected seed number */
} args_reduce;

/* ======================================================================== */
typedef struct {
    char seqfile[BUF512];  /* sequence file */
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];  /* running log */
    char ids_file[BUF512];  /* list of IDs to extract sequences for */
    char filename[BUF512];  /* output file name */

    long rc;  /* write out reverse complementary sequence */
    long reverse;  /* just get the reverse of the sequence */
    long whole_line;  /* write out sequence in one line */
    long sort;  /* switch to sort sequence entries, can't be combined with -ids */
    long base_case;
    long hjb;  /* convert to the simplified HJB format */
} args_process_fasta;

/* ======================================================================== */
#define IDS_CLIST "ids_clist.txt"
typedef struct {
    char tdatfile[BUF512];  /* tab-delimited data file */
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];  /* running log */
    char ids_file[BUF512];  /* list of IDs to extract data for */
    char ids_list[BUFBIG];  /* list of IDs from command line: ',;' delimited */
    char columns[BUF512];  /* list of columns to extract */
    char filename[BUF512];  /* output file name */
    char log_base[BUF512];  /* base of log-transformation */

    long sort;  /* switch to sort entries by ID */
    long id_case;
    long transpose;  /* get transpose of tdat file */
} args_process_tdat;

/* ======================================================================== */
typedef struct {
    char htmlfile[BUF512];
    char psam[BUF512];  /* get an index page of a list of PSAMs */
    char outdir[BUF512];  /* output directory from a MatrixREDUCE run */
    char runlog[BUF512];

    double tval_cutoff;

    long tn_width;  /* width of the thumbnail image in pixel */
    long tn_height;  /* height in pixel */
    long copy;  /* switch to copy .css/.js to output directory */
    long rc;  /* switch to generate logo corresponding to reverse comp strand */
} args_html_summary;

/* ======================================================================== */
typedef struct {
    char source[BUF512];  /* source of the PSAM: i.e. MatrixREDUCE */
    char topology[BUF512];  /* topological pattern this PSAM is based on */
    char seed_motif[BUF512];
    char measurement_file[BUF512];
    char experiment_name[BUF512];

    long experiment_column;
    long bonferroni;
    double p_value;

    long strand;
    long psam_length;
    char optimal_sequence[BUF512];

    double *psam;
} struct_psam;

/* ======================================================================== */
/* only used in HTMLSummary.c: to be consolidated with struct_fitpars? */
typedef struct {
    char experiment_name[BUF512];
    long experiment_column;
    long num;
    double r2;
    double ave;
    double var;
    double fit_var;
    double res_var;
    double F[BUFBIG];
    double t[BUFBIG];
    double P[BUFBIG];
} struct_ftp;

/* ======================================================================== */
/* only used in HTMLSummary.c */
typedef struct {
    char experiment_name[BUF512];
    long experiment_column;
    long num;
    long is_seed;
    double r2;
    double F;
    double t;
    double P;
} struct_spar;

/* ======================================================================== */
typedef struct {
    char file[BUF512];  /* input file */
    char format[BUF512];  /* EPS | PDF | JPEG | PNG | GIF */
    char type[BUF512];  /* PSAM | PWM | fasta | flat */
    char logo[BUF512];  /* name of output image file */
    char outdir[BUF512];  /* name of output directory */
    char runlog[BUF512];
    char title[BUF512];  /* title of the logo image */
    char subtitle[BUF512];  /* subtitle if -consensus not set; following Gabor */
    char style[BUF512];  /* logo style: raw | freq | bits_info | ddG */

    double min_Ka;  /* minimum allowable Ka: apply only when format=matrix */
    double width;  /* logo width in cm */
    double height;  /* logo height in cm */

    long smallSampleCorrection;  /* switch */
    long errorBar;  /* switch */
    long stretch;  /* switch */
    long bw;  /* switch for using black-and-white */
    long box;  /* 0 for none; 1 for fill; -1 for empty rectangular */
    long outline;  /* switch for showing characters in outline */
    long frame;  /* switch for drawing a bounding box around image */
    long consensus;  /* switch to show consensus sequence as fineprint */
    long rc;  /* switch to show reverse complementary strand */

    long start_num;  /* sequence label start number, default to 1  */
    long sb;  /* beginning sequence # to be shown */
    long se;  /* ending sequence # to be shown */

    long margin;  /* stack margin between base letters: default to 1 */
    double yaxis;  /* scale of y-axis: default to auto */
    double ytick;  /* tick of y-axis: defualt to auto */
    double ymin;  /* set the range [ymin, ymax] */
    double ymax;

    long label;  /* switch to turn on/off label */
} args_logo_generator;

typedef struct {
    char base;
    double height;
} struct_base_height;

typedef struct {
    double num_valid;  /* valid number of bases in position/column */

    double raw_value[NUM_BASE4];  /* row count or PSAM value for each of ACGT base */
    double frequency[NUM_BASE4];  /* count[i] / num_seqs */

    double e;  /* amount for small sample correction */
    double H;  /* -sum(f[b,l] * log2(f[b, l])) */
    double R;  /* amount of information at position 'l' */

    double height[NUM_BASE4];  /* base letter height in logo */
    struct_base_height height_sorted[NUM_BASE4];  /* sorted in order */

    double yoffset;  /* offset for negative value; 0 otherwise */
    long ymin;  /* global bottom value for y-axis */
    long ymax;  /* global top value for y-axis */
} struct_logo_info;

/* ======================================================================== */
#define MODEL_FILE_LIST "coeff.tsv", "tvals.tsv", "pvals.tsv"
enum {
    FVAL_IDX = 0,
    TVAL_IDX = 1,
    PVAL_IDX = 2
};

#define LOGO_PS_PAR "LogoGenerator_PS.def"
#define LOGO_PS_TMP "temp_logo.eps"
#define LOGO_OPTIONS "LogoGenerator.opt"
#define LOGO_DATA "logoData.out"

#define GOLDEN_RATIO 1.618
#define TN_WIDTH 145
#define TN_HEIGHT 90

#define MATRIX "MatrixREDUCE"
#define MOTIF "MotifREDUCE"

#define NOTE "<span class=\"note\">Note:</span>"
#define CMENU "\n<p class=\"cmenu\">\
\n<a href=\"#wrap\">Top</a> | \
\n<a href=\"#recap\">Recap</a> | \
\n<a href=\"#thumbnail\">PSAM Thumbnail</a> | \
\n<a href=\"#matrices\">PSAM Logo</a> | \
\n<a href=\"#ftable\">F-value</a> | \
\n<a href=\"#ttable\">t-value</a> | \
\n<a href=\"#ptable\">P-value</a> | \
\n<a href=\"#summary\">Expression Summary</a>\
</p>"

/* ======================================================================== */
typedef struct {
    long strand;
    long column;
    long normalize;
    double threshold;
    long relative;  /* is the threshold relative to the PSAM maximum? */
    long sum_only;  /* only output sum of affinity per gene */
    long rc_match;  /* match to reverse complementary strand: miReduce */
    long max_only;  /* only output maximum affinity per gene, per Ron's request */

    char chromosome[BUF512];
    char seqfile[BUF512];
    char prefix[BUF512];
    char outdir[BUF512];
    char runlog[BUF512];
    char ids[BUFBIG];  /* a comma-delimited list */
    char affsum[BUF512];  /* accumulative affinity per sequence */

    char psam[BUF512];
    char psam_list[BUF512];
    char motif[BUF512];
    char motif_list[BUF512];
} args_affinity_profile;

/* ======================================================================== */
typedef struct {
    long strand;
    long damid;
    long copy;
    long xcheck;  /* check for multiple linear combination of X variables */
    long univariate;  /* univariate fit only */
    long acgt;  /* base composition correction */
    long transpose;  /* swap row and column */
    long intercept;  /* if to show intercept; used with 'transpose' */
    double tval_cutoff;

    char seqfile[BUF512];
    char measfile[BUF512];
    char outdir[BUF512];
    char runlog[BUF512];
    char design_matrix[BUF512];
    char resid_file[BUF512];

    char psam[BUF512];
    char psam_list[BUF512];
    char motif[BUF512];
    char motif_list[BUF512];
} args_transfactivity;

/* ======================================================================== */
typedef struct {
    char seqfile[BUF512];
    char measfile[BUF512];
    char motif[BUF512];
    char psam[BUF512];
    char outdir[BUF512];
    char runlog[BUF512];
    char filename[BUF512];
    char optpos[BUF512];

    long strand;
    double p_value;
} args_optimize_psam;

/* ======================================================================== */
typedef struct {
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];
    char dictfile[BUF512];  /* output file name */

    char topo[BUF512];  /* single topology pattern */
    char topo_list[BUF512];  /* name of topology file */

    char expand[BUF512];  /* expand a single generic pattern */
    char expand_list[BUF512];  /* expand list of generic patterns */

    long base_case;
    long keep_iupac;
} args_topo2dictfile;  /* topology to motif dictionary file */

/* ======================================================================== */
typedef struct {
    char *seqid;  /* sequence ID */
    char *winid;  /* window entry ID */
    long strand;  /* +1 or -1 */
    long ibeg;  /* starting base position */
    long iend;  /* ending base position */
    long width;  /* base fragment length */
} struct_win;

typedef struct {
    long num_wins;
    struct_win *windows;
} struct_winArr;

/* ======================================================================== */
typedef struct {
    char seqfile[BUF512];  /* sequence file */
    char winfile[BUF512];  /* windows file */
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];
    char filename[BUF512];  /* output file name */

    long width;
    long base_case;
} args_extract_windows;

/* ======================================================================== */
typedef struct {
    char source[BUF512];  /* type of source to be converted */
    char inpfile[BUF512];  /* input file name */
    char outdir[BUF512];  /* directory for output */
    char runlog[BUF512];
    char psamfile[BUF512];  /* output PSAM file name */
    char pwmfile[BUF512];  /* for 'standard' PWM format' */
    char style[BUF512];  /* raw | freq | normalized */
    char idstr[BUF512];  /* id string for PWM or <source> tag in PSAM */

    long strand;
} args_convert2psam;

/* ======================================================================== */
/* Levenberg-Marquardt ls-fitting parameters */
typedef struct {
    double delta_g;
    double delta_x;
    long max_iter;
    long min_oknum;
    long num_print;
} struct_levmar;

#include "REDUCE_Suite_fncs.h"

#endif  /* _REDUCE_SUITE_H */
