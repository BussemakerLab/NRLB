
/* AffinityProfile.c */
void affinity_profile_cmdline(int argc, char *argv[], args_affinity_profile * args);
void psam2profile(char *src, struct_seqArr * seqs, struct_topo * t, long strand,
                  double *w, long num_ids, long *ids, double pmax,
                  args_affinity_profile * args, double *csum, double threshold);
void get_psam_count_profile(long num_ids, long *ids, struct_seqArr * seqs,
                            args_affinity_profile * args);
void get_motif_count_profile(long num_ids, long *ids, struct_seqArr * seqs,
                             args_affinity_profile * args);
long extract_id_list(long *ids, char *ids_list, struct_seqArr * seqs);
void check_chromosome(args_affinity_profile * args, time_t time0);
void profiling_chromosome(args_affinity_profile * args);

/* Convert2PSAM.c */
void cvt2psam(args_convert2psam * args);
void cvt_multi_align_seq_psam(args_convert2psam * args, struct_psam * psam_data);
void cvt_pwm_psam(args_convert2psam * args, struct_psam * psam_data);
void cvt_tr_psam(args_convert2psam * args, struct_psam * psam_data);
long populate_transfac_counts(FILE * fp, double *psam);
long populate_jaspar_counts(char *str, int base, double *psam);
long is_jaspar_base_row(int base, char *str);
void cvt_ja_psam(args_convert2psam * args, struct_psam * psam_data);
void check_for_negative(long num, double *psam);
void get_psam_freq(long num, double *psam);
long read_bf_psam(char *inpfile, char *title, double *psam, long *strand);
long check_bf_file(char *inpfile);
long read_v1_psam(char *inpfile, char *title, double *psam, double *p_value, long *bonf,
                  long *strand);
void check_for_space_pound_substring(char *line, char *inpfile, char *line_numID);
long check_v1_file(char *inpfile);
void get_acgt_counts(long num_seqs, struct_seq * sequences, long num_nt, double *psam);
void extract_tamo_source(char *str, char *source);
long populate_tamo_counts(char *str, int base, double *psam);
void cvt_tamo_psam(args_convert2psam * args, struct_psam * psam_data);
void c2p_cmdline(int argc, char *argv[], args_convert2psam * args);

/* ExtractWindows.c */
void ew_cmdline(int argc, char *argv[], args_extract_windows * args);

/* HTMLSummary.c */
void check_psam_files(FILE * fp, long num, struct_tag * str_tags);
void get_psam_logo(long idx, char *psam_file, char *imgfile, char *html_opts,
                   char *bindir, char *outdir);
void table_psam_logo(FILE * fp, long num, char *psam_file, char *imgfile, long width,
                     long height);
void get_psams_index(FILE * fp, char *html_opts, char *img_type, char *bindir,
                     char *outdir, char *psam);
struct_ftp *allocate_memory_for_struct_ftp(long num);
void free_struct_ftp(struct_ftp * statistic_info);
struct_spar *allocate_memory_for_struct_spar(long num);
void free_struct_spar(struct_spar * sinfo);
void get_reduce_html(FILE * fp, char *html_opts, char *img_type, char *bindir,
                     args_html_summary * args);
char output_ftp_header(FILE * fp, long idx);
void set_FtPTable(FILE * fp, long num_expts, long psam_num, struct_ftp * statistic_info,
                  double tc);
void seed_psam_table(FILE * fp, long pnum, char *psam_file, char *imgfile, long width,
                     long height, struct_psam * psam_data, struct_spar * sinfo_idx);
void seed_psam_expt_list(FILE * fp, long idx, long num_expts, struct_spar * sinfo);
void set_psamImages(FILE * fp, long num_expts, long psam_num, char *html_opts,
                    char *img_type, char *bindir, char *outdir);
void get_width_height(char *outdir, long *width, long *height);
void get_html_opts(char *html_opts, char *img_type);
void set_usedOptions(FILE * fp, char *outdir);
void set_exprSummary(FILE * fp, char *outdir, char *measfile, long num_expts,
                     struct_ftp * statistic_info);
void set_thumbnailImages(FILE * fp, long psam_num, char *html_opts, char *img_type,
                         char *bindir, args_html_summary * args);
void check_valid_measurement(char *chk_file, long *num_missing, long *num_duplicate);
void check_valid_sequence(char *seqfile, long *num_seqs, long *num_repeat);
void check_htmlfile(char *outdir, char *htmlfile);
void set_htmlsummary_defaults(args_html_summary * args);
void check_tn_width_height(double w, double h, args_html_summary * args);
void html_summary_cmdline(int argc, char *argv[], args_html_summary * args);
void extract_basic_info(char *outdir, long *num_expts, long *psam_num, char *seqfile,
                        char *measfile);
void get_statisticInfo(char *outdir, long num_expts, long psam_num,
                       struct_ftp * statistic_info);

/* LogoGenerator.c */
void logo_psam(args_logo_generator * args);
void logo_pwm(args_logo_generator * args);
void logo_seq(args_logo_generator * args);
void get_image(args_logo_generator * args);
void get_cfg(char *gs, char *gs_opts, char *convert, char *fmt);
void set_heightByStyle(args_logo_generator * args, struct_logo_info * pld);
void get_y_range(args_logo_generator * args, long num_pos, struct_logo_info * logoData);
double get_maxSum(long num_pos, struct_logo_info * logoData);
double get_minSum(long num_pos, struct_logo_info * logoData);
void get_ddG_height(struct_logo_info * pld, args_logo_generator * args);
void write_matrix_info(FILE * fp, char *consensus_seq, long num_pos, char *style,
                       struct_logo_info * logoData, char *fmt, char *fmt0);
void write_seq_info(FILE * fp, char *consensus_seq, long num_pos, long num_seqs,
                    struct_seq * sequences, struct_logo_info * logoData, char *fmt,
                    char *fmt0);
void verify_logoData(char *consensus_seq, long num_pos, struct_logo_info * logoData,
                     long num_seqs, struct_seq * sequences, args_logo_generator * args);
void print_dval(FILE * fp, double d, char *fmt, char *fmt0);
void generate_eps(long num_pos, struct_logo_info * logoData, long num_seqs,
                  struct_seq * sequences, args_logo_generator * args);
void get_consensus_seq(long num_pos, struct_logo_info * logoData, char *consensus_seq);
void get_logo_height(FILE * fp, args_logo_generator * args, long num_pos,
                     struct_logo_info * logoData);
void set_ytick(FILE * fp, long yscale, long y_overall_offset);
void set_margins(FILE * fp, long label);
void get_ps_settings(FILE * fp, args_logo_generator * args, char *consensus_seq,
                     long yscale, long y_overall_offset);
void get_ps_header(FILE * fp, args_logo_generator * args);
void get_logo_data(long num_pos, struct_logo_info * logoData, long num_seqs,
                   struct_seq * sequences, args_logo_generator * args);
void get_ssc(long smallSampleCorrection, struct_logo_info * pld);
void get_frequency(struct_logo_info * pld);
void get_H(struct_logo_info * pld);
void get_R(long stretch, struct_logo_info * pld);
void get_height(struct_logo_info * pld);
int height_compare(const void *v1, const void *v2);
void sort_height(struct_logo_info * pld);
void init_struct_logo_info(struct_logo_info * pld);
void set_logo_defaults(args_logo_generator * args);
void check_logo_cmdline(args_logo_generator * args);
void logo_cmdline(int argc, char *argv[], args_logo_generator * args);

/* MatrixREDUCE.c */

/* MotifREDUCE.c */

/* OptimizePSAM.c */

/* ProcessFASTA.c */

/* ProcessTdat.c */
void tdat_cmdline(int argc, char *argv[], args_process_tdat * args);

/* Topo2Dictfile.c */

/* Transfactivity.c */
void transfactivity_cmdline(int argc, char *argv[], args_transfactivity * args);
void handle_design_matrix(args_transfactivity * args, struct_data * tdat);
void link_design2tdat_by_ids(char *outdir, char *filename, struct_data * tdat,
                             struct_data * fit_tdat);
void get_matched_design_tdat(struct_data * fit_tdat, struct_data * tdat);
void handle_motif(args_transfactivity * args, struct_seqArr * seqs, struct_data * tdat);
void handle_psam(args_transfactivity * args, struct_seqArr * seqs, struct_data * tdat);
void get_transfactivity_type(args_transfactivity * args, char *type);
void perform_fitting(args_transfactivity * args, long num, struct_tag * str_tags,
                     long num_Xrow, double **Xmtx, struct_data * tdat);

/* fncs_apps.c */
long get_psam_strand(long cmd_strand, long psam_strand);
long get_motif_strand(long cmd_strand);
void log_strand_msg(long strand);
void display_strand_msg(FILE * fp, char *prefix, char *suffix, long strand);
char *strand2msg(long strand);
long msg2strand(char *msg);
void set_residual_file(char *resid_file, char *measfile, char *outdir);
void get_residuals(struct_data * tdat, long pnum, double **X, struct_fitpars * lfpars);
void get_predicted(struct_data * tdat);
double get_slope_abs_tvalue(long n, double sX, double sX2, double sXY, double sY,
                            double sY2);
void lsfit_parameters(long n, double sX, double sX2, double sXY, double sY, double sY2,
                      struct_ls5fit * fitpars);
double get_y_vs_x_slope_abs_tvalue(long nrow, double *x, double *y);
void get_y_vs_x_lsfitpars(long nrow, double *x, double *y, struct_ls5fit * fitpars);
long get_num_mids(long num, double *counts, double *y);
void print_fitpars(char *str, struct_ls5fit * fitpars);
void get_psam_fitpars(struct_data * tdat, struct_seqArr * seqs, struct_seed * seed,
                      struct_ls5fit * fitpars, long debug);
double get_count_per_window(char *p0, struct_topo * t, long strand, double *w);
double get_max_pcount_seq(struct_seq * seq, struct_topo * t, long strand, double *w);
double csum_per_seq(struct_seq * seq, struct_topo * t, long strand, double *w,
                    long normalize, double pmax, double threshold);
void count_per_sequence(struct_seq * seq, struct_topo * t, long strand, double *w,
                        long normalize, double pmax, double threshold, double Nsum,
                        long column_wise, FILE * fp);
double get_pcount_seq(struct_seq * seq, struct_topo * t, long strand, double *w);
void get_psam_counts(struct_seqArr * seqs, struct_topo * t, long strand, double *w,
                     double *counts);
void get_whole_counts(char *motif_seq, long strand, struct_seqArr * seqs,
                      double *motif_count);
long svd_check_psam_counts(long ndat, long npar, char *type, double **X);
void write_fit_yX(long icol, long nrow, long ncol, double *y, double **X, long chk);
long get_number_of_matched_ids(long num, double *counts);
void output_model_summary(args_reduce * args, long num_motifs, struct_seed * smotifs,
                          struct_data * tdat, struct_seqArr * seqs,
                          struct_fitpars * lfpars);
void perform_multifit(struct_data * tdat, long na, double **psam_counts,
                      struct_fitpars * lfpars, long bonf, long chk);
void process_sig_seed(long pnum, struct_data * tdat, struct_seqArr * seqs,
                      struct_seed * seed, double **psam_counts, struct_fitpars * lfpars);
long psam_motif_sig_output(char *type, double E_value, double cutoff);
long is_sig_motif(struct_seed * seed, double p_value);
void free_unsig_seed(double *counts, struct_seed * seed);
long isok_motif(double p_value, struct_data * tdat, struct_seqArr * seqs,
                struct_seed * seed, long *psam_num, double **psam_counts,
                struct_fitpars * lfpars);
long is_sig_ron_correction(double N, long Xcount, double abs_r, double p_value);
long isok_psam(double p_value, struct_data * tdat, struct_seqArr * seqs,
               struct_seed * seed, long *psam_num, double **psam_counts,
               struct_fitpars * lfpars);
void reset_Xmtx(long num_ok, long num0, long *num, long *is_okay, long num_seqs,
                double **counts, struct_tag * str_tags);
void check_for_zeros(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx);
void check_for_degeneracy(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx);
void check_for_multifit(long *num, struct_tag * str_tags, long num_Xrow, double **Xmtx);
void write_design_matrix(struct_data * tdat, long num, double **psam_counts,
                         struct_tag * str_tags, char *type, char *outdir);
void calculate_fitpars(struct_data * tdat, struct_ls5fit * fitpars);
void validate_motif(long *num, struct_tag * str_tags);
void cat_dir_fname(char *outdir, char *fname, char *fullname);
void log_prg(char *msg);
void log_run(char *msg);
void log_msg(char *msg);
void log_msg_vchk(char *msg);
void log_msg_exit(char *msg);
void get_topo_or_dict(char *outdir, char *dictfile, char *topo_one, char *topo_list,
                      struct_topoArr * topo, struct_dic * mlist);
void find_seed(char *dictfile, struct_dic * mlist, struct_topoArr * topo,
               struct_data * tdat, struct_seqArr * seqs, long ntop, long nsel,
               struct_seed * seed);
void set_reduce_defaults(args_reduce * args);
void write_reduce_options(args_reduce * args);
void check_reduce_cmdline(args_reduce * args);
void reduce_cmdline(int argc, char *argv[], args_reduce * args);
void output_reduce_results(args_reduce * args, long psam_num, struct_seed * smotifs,
                           struct_data * tdat, struct_fitpars * lfpars, char *type);

/* fncs_bits.c */
void convert_qua2seq(long *c, long num, char *seq);
void convert_dec2seq(long idx, long num, char *seq);
void fast_d2q(unsigned long n, long *c);
unsigned long fast_q2d(long *c);
long fast_str_q2d(char *p);
void reverse_cmpl(char *seq, long debug);
unsigned long fast_cidx(unsigned long n, long Xcount);
unsigned long turnon_bits(long startbit, long numbits);
unsigned long turnoff_bits(long startbit, long numbits);
unsigned long toggle_bits(unsigned long word, long startbit, long numbits);
long testbit(unsigned long word, long bit_to_test);
void printbits(unsigned long word);

/* fncs_cmns.c */
long is_palindrome(char *str);
long is_rc_palindrome(char *str);
void reverse_string(char *str);
long string_contains_only_those_characters(char *str, char *chars_set);
long is_numeric(char *str);
double cvt2double(char *str);
long cvt2long(char *str);
void quit_with_help_reminder(char *str);
long has_no_equal_sign(char *parstr);
long get_lvalue(char *str, long vmin, long vmax);
double get_dvalue(char *str, double vmin, double vmax);
void get_strvalue(char *str, char *dst, long expand_tilde);
long str_pmatch(char *str, char *sstr);
long case_str_pmatch(char *str, char *sstr);
void print_used_time(time_t time0);
void get_currentTimeString(char *s);
void get_currentYear(char *year);
double my_log2(double d);
char *enlarge_cline(long *maxline, char *line);
char *my_getline(FILE * fp);
long readline_cvt2long(FILE * fp);
void cvtstr_set1toc2(char *str, char *set1, char c2);
void cvtstr_c1toc2(char *str, char c1, char c2);
long csplit(char *str, char *items[], long itemsize, char sepc);
long item_list(char *str, char *items[], long itemsize, char *sep_chars);
FILE *open_tmpfile(void);
long is_std_out_err(char *filename);
FILE *open_file(char *filename, char *filemode);
long close_file(FILE * fp);
long exist_file(char *filename);
void remove_file(char *filename);
long is_same_file(char *src, char *dst);
void rename_file(char *src, char *dst);
void copy_file(char *src, char *dst);
void cat_file(char *src, char *dst);
void get_alnum(char *a, char c);
char *ltrim(char *a);
char *rtrim(char *a);
char *trim(char *a);
void upperstr(char *a);
void lowerstr(char *a);
int case_strcmp(const char *s1, const char *s2);
int case_strncmp(const char *s1, const char *s2, long n);
char *case_strstr(const char *haystack, const char *needle);
char *case_strchr(const char *s, int c);
char *my_strdup(const char *src);
void print_sep(FILE * fp, char c, long num);
void repeat_char_ntimes_string(char c, long n, char *str);
long num_strmatch(char *str, char **strmat, long ib, long ie);
void add_end_slash(char *str);
void delete_end_slash(char *str);
char *basename(char *str);
void del_extension(char *fullname, char *okname);
void bname_noext(char *src, char *dst);
void bname_ext(char *src, char *ext, char *dst);
void get_stropt_wo_end_slash(char *option, char *okstr);
long lround(double d);
long dbl2long(double dval);
void fatal(char *fmt, ...);
void change_case_str_tags(long upper_case, long num, struct_tag * str_tags);
struct_tag *fillup_str_tags(long is_list, char *filename, long *num);
struct_tag *read_unique_strtags(char *idfile, long *num);
void init_strtag_with_strings(char *str, char *tag, struct_tag * s);
void init_strtag(struct_tag * src, struct_tag * dst);
void copy_strtag(struct_tag * src, struct_tag * dst);
void swap_strtag(struct_tag * one, struct_tag * two);
void get_unique_strtag(long *num, struct_tag * str_tags);
void read_strtags(char *filename, struct_tag * str_tags);
void set_strtag(char *line, struct_tag * str_tag);
void skip_lines(FILE * fp, long num);
long get_line_number(char *filename, long skips);
long is_empty_string(const char *str);
long is_equal_string(const char *str1, const char *str2);
long is_equal_case_string(const char *str1, const char *str2);
void assert_equal_longs(char *msg, long a, long b);
void set_parallel_idx(long ib, long ie, long *idx);
void nullify_line_comment(char *str);
long is_comment_line(char *line);
long is_skip_line(char *line);
long is_fasta_skip_line(char *line);
void transpose_matrix(double **src, long nr, long nc, double **dst);
void print_dvector(double *dvec, long nl, long nh, char *msg);
void do_nothing(void);

/* fncs_html.c */
void output_html_spar_fitpars(long psam_num, struct_data * tdat,
                              struct_seqArr * seqs, struct_seed * seed);
void read_spar_fitpars(long psam_num, long num_expts, char *logfile, struct_spar * sinfo);
double get_valFtP(long idx, double F, double t, double P);
void print_valFtP(FILE * fp, long idx, double F, double t, double P);
void write_FtP(char *type, char *outdir, char *expr_bname, struct_data * tdat,
               long psam_num, struct_fitpars * lfpars);
void write2raw(char *filename, args_reduce * args, struct_data * tdat,
               long psam_num, struct_seed * smotifs, struct_fitpars * lfpars);
void write_psams_list(long motif_num, struct_seed * smotifs, char *outdir);
void write_motifs_list(long motif_num, struct_seed * smotifs, char *outdir);
void write_list_resid_predict(long psam_num, struct_seed * smotifs,
                              char *outdir, struct_data * tdat);
void output_mvar_fitpars(struct_data * tdat, long num, char *type, char *outdir,
                         struct_ls5fit * fitpars, struct_fitpars * lfpars);
void write_separate_FtP_HTML(FILE * fp, long intercept, long idx, double F, double t,
                             double P, double tc);
void open_table(FILE * fp, char *class);
void close_table(FILE * fp);
void html_recap_common(FILE * fp, char *options, char *seqfile, char *measfile,
                       struct_data * tdat);
void set_transfactivity_univar_htmlHeader(FILE * fp, char *seqfile, struct_data * tdat,
                                          char *measfile);
void set_transfactivity_htmlHeader(FILE * fp, long num, struct_tag * str_tags,
                                   char *SUBMENU, char *type, char *seqfile,
                                   struct_data * tdat, char *measfile, char *resid_file);
void set_overall_fit(FILE * fp, long num, char *SUBMENU, struct_data * tdat,
                     struct_ls5fit * fitpars, struct_fitpars * lfpars);
void mvar_FtP_transposed(FILE * fp, long FtP_idx, long ncol, long num, long intercept,
                         struct_fitpars * lfpars, struct_tag * str_tags, double tc);
void mvar_FtP_normal(FILE * fp, long FtP_idx, struct_data * tdat,
                     long num, struct_fitpars * lfpars, char *type, double tc);
void output_mvar_fitparsHTML(struct_data * tdat, long num, char *type,
                             struct_ls5fit * fitpars, struct_fitpars * lfpars,
                             struct_tag * str_tags, args_transfactivity * args);
void output_univar(struct_data * tdat, long num, char *type, double **Xmtx,
                   struct_tag * str_tags, char *outdir);
void output_univarHTML(struct_data * tdat, long num, char *type, double **Xmtx,
                       struct_tag * str_tags, args_transfactivity * args);
void output_univarHTML_multiple_expt(FILE * fp, struct_data * tdat, long num, char *type,
                                     double **Xmtx, struct_tag * str_tags);
void output_univarHTML_one_expt(FILE * fp, struct_data * tdat, long num, char *type,
                                double **Xmtx, struct_tag * str_tags);
void set_htmlFileHeader(FILE * fp, char *outdir, long copy, char *pgname);
void set_htmlFileFooter(FILE * fp);
void copy_gifs_css_js(long copy, char *fileLoc, char *outdir, char *css, char *js);
long is_matrix_reduce_run(char *outdir);
void get_raw_data_file(char *outdir, char *filename);
void get_default_html_file(char *outdir, char *filename);

/* fncs_memo.c */
struct_win *allocate_memory_for_windows(long num);
void free_windows(struct_winArr * wins);
struct_id *allocate_memory_for_idstr(long num);
void free_idstr(long ib, long ie, struct_id * idstrs);
struct_tag *allocate_memory_for_strtags(long num);
void free_strtags(long ib, long ie, struct_tag * str_tags);
void free_p_strtag(struct_tag * str_tag);
void free_extra_strtags(long ib, long ie, struct_tag * str_tags);
char **allocate_char_marray(long num);
void free_char_marray(long ib, long ie, char **pp);
long **allocate_long_marray(long num);
void free_long_marray(long ib, long ie, long **pp);
double **allocate_double_marray(long num);
void free_double_marray(long ib, long ie, double **pp);
struct_ntop *allocate_memory_for_ntop_entries(long num);
void free_ntop_entries(struct_ntop * ntops);
struct_seed *allocate_memory_for_smotifs(long num);
void free_smotifs(long ib, long ie, struct_seed * smotifs);
void free_seed(struct_seed * seed);
struct_fitpars *allocate_memory_for_lfpars(long num, long vsize);
void free_lfpars(long ib, long ie, struct_fitpars * lfpars);
struct_ls5fit *allocate_memory_for_fitpars(long num);
void free_fitpars(struct_ls5fit * fitpars);
void init_mlist(struct_dic * mlist);
void free_mlist(struct_dic * mlist);
void free_resources(long psam_num, double **psam_counts, struct_seed * smotifs,
                    struct_data * tdat, struct_seqArr * seqs, struct_topoArr * topo,
                    struct_fitpars * lfpars, struct_dic * mlist);
void initialize_psam_struct(struct_psam * psam_data);
void free_psam_data(struct_psam * psam_data);
struct_logo_info *allocate_memory_for_logoData(long num);

/* fncs_misc.c */
void show_version();
void get_sysdir(char *dstdir, char *dirname, char *filename);
void get_sysfile(char *fullfile, char *dirname, char *basefile);
void display_helpfile(char *helpfile);
long read_pwm(char *pwmfile, double *w);
void populate_pd_struct_from_psam(long num_nt, double *psam, long strand0, double p_value,
                                  struct_psam * psam_data);
void debug_parenthesis(char *pat, long num_parens, long *parens_bidx, long *parens_eidx);
void write2xml(char *filename, args_reduce * args, struct_data * tdat,
               long psam_num, struct_seed * smotifs, struct_fitpars * lfpars);
void xml_parameters(FILE * fp, args_reduce * args, long expt_num, long psam_num);
void xml_experiment(FILE * fp, struct_data * tdat, long psam_num, struct_fitpars * lfpars);
void xml_fitpars(FILE * fp, struct_data * tdat, long psam_num, long icol,
                 struct_fitpars * lfpars);
void xml_FtP_value(FILE * fp, double f, double t, double p);
void xml_psam(FILE * fp, long psam_num, struct_seed * smotifs);

/* fncs_mylm.c */
void populate_A_g_chisq(long ma, double *a, long *ia, struct_lm * lmdata, double *g,
                        double **A, double *chisq);
double my_fdf(struct_lm * lmdata, long idx, long nWs, double *a2, double *da2,
              double F, double *dyda);
void update_A_g_chisq_per_data_point(long ma, long *ia, double *dyda, double dy,
                                     double *g, double **A, double *chisq);
void my_mrqcof(long ma, long mfit, double *a, long *ia, struct_lm * lmdata, double *g,
               double **A, double *chisq, double *ng);
double get_position_Wmax(double *pos);
void do_mylm(long ma, double *a, struct_lm * lmdata);
void init_lmseed(long num_pars, double W0, struct_seed * seed);
void check_opt_allpos(struct_topo * t);
void lm_optimize_seed(struct_seed * seed, struct_data * tdat, struct_seqArr * seqs);
void lm_optimize_expt_psam(struct_seed * seed, struct_data * tdat, struct_seqArr * seqs);

/* fncs_nrcs.c */
void NRC_avevar(double data[], long n, double *ave, double *var);
double NRC_betai(double a, double b, double x);
double NRC_gammln(double xx);
double NRC_betacf(double a, double b, double x);
double NRC_erff(double x);
double NRC_gammp(double a, double x);
double NRC_gammq(double a, double x);
void NRC_gcf(double *gammcf, double a, double x, double *gln);
void NRC_gser(double *gamser, double a, double x, double *gln);
void NRC_fit(double x[], double y[], long ndata, double *a, double *b, double *siga,
             double *sigb, double *chisq);
void NRC_lfit(double **X, double *y, long ndat, long ma, double *a, double **covar,
              double *chisq);
void NRC_my_svdfit(double **X, double y[], long ndat, long ma, double a[],
                   double **covar, double *chisq);
void NRC_svdfit(double **X, double y[], long ndat, long ma, double a[],
                double **u, double **v, double w[], double *chisq);
void NRC_svdvar(double **v, long ma, double w[], double **covar);
void NRC_gaussj(double **a, long n, double **b, long m);
void NRC_covsrt(double **covar, long ma, long ia[], long mfit);
void NRC_pearsn(double x[], double y[], long n, double *r, double *prob, double *z,
                double *tval);
void NRC_ludcmp(double **a, long n, long *indx, double *d);
void NRC_lubksb(double **a, long n, long *indx, double b[]);
void NRC_lu_inv(double **a, long n, double **y);
void NRC_mprove(double **a, double **alud, long n, long indx[], double b[], double x[]);
void NRC_svdcmp(double **a, long m, long n, double w[], double **v);
void NRC_svbksb(double **u, double w[], double **v, long m, long n, double b[],
                double x[]);
void NRC_svdzwi(double *w, long n);
void NRC_svd_inv(double **u, double w[], double **v, long m, long n, double **ginv);
double NRC_pythag(double a, double b);
long NRC_chol_posdef(double **a, long n, double p[]);
void NRC_choldc(double **a, long n, double p[]);
void NRC_cholsl(double **a, long n, double p[], double b[], double x[]);
void NRC_chol_inv(double **a, long n, double p[], double **Ainv);
void NRC_chol_L(double **a, long n, double p[], double **L);
void NRC_chol_Linv(double **a, long n, double p[], double **Linv);
void NRC_chol_Ainv(double **a, long n, double p[], double **Ainv);
void NRC_qrdcmp(double **a, long n, double *c, double *d, long *sing);
void NRC_qrsolv(double **a, long n, double c[], double d[], double b[]);
void NRC_qr_inv(double **a, long n, double c[], double d[], double **ainv);
double NRC_ran1(long *idum);
double NRC_ran2(long *idum);
double NRC_ran3(long *idum);
double NRC_gasdev(long *idum);
void lux_Axsolution(long m, double **A0, double *g, double mu, double *h);
void get_covar(long m, double **A0, double **covar);
void NRC_eigsrt(double *d, double **v, long n);
void NRC_jacobi(double **a, long n, double *d, double **v);

/* fncs_opts.c */
void extract_log_base(char *option, char *log_base);
long extract_option_psam_motif(char *option, char *psam, char *psam_list, char *motif,
                               char *motif_list);
long extract_option_topo(char *option, char *topo, char *topo_list);
long extract_verbose_option(char *option);
long extract_case_option(char *option);
void write_option_case(FILE * fp, long case_setting);
void write_option_outdir_runlog(FILE * fp, char *outdir, char *runlog);
void write_option_psam_motif(FILE * fp, char *psam, char *psam_list, char *motif,
                             char *motif_list);
void write_option_topo(FILE * fp, char *topo, char *topo_list);
void write_option_expand(FILE * fp, char *expand, char *expand_list);
long extract_strand_option(char *option);
long extract_strand_option_no0(char *option);
void check_topo_option(char *topo_list, char *topo);
void set_my_globals(char *pgname);
void cleanup_my_globals();
long check_common_options(char *option, char *helpfile, char *runlog, char *outdir);
void display_unrecognized_option(char *option);
void set_logs_and_check_option(char *runlog, char *outdir, long argc, long end_opt_idx,
                               char *option);
void set_logs(char *runlog, char *outdir);
void check_option_format(long argc, long end_opt_idx, char *option);
void check_either_A_or_B(char *optA, char *optB, char *msg);
void check_required_option(char *option, char *invalid_str, char *msg);
void check_required_file(char *filename, char *invalid_str, char *msg);
long set_switch_with_dft_true(char *option);
long set_switch_with_dft_false(char *option);

/* fncs_psam.c */
void get_full_psam(struct_seed * seed);
char posvec_iupac(double *v0);
void get_psam_optimal_seq(double *psam, long num, char *optimal);
void write_psam_data(FILE * fp, long Tcount, char *opt_seq, char *full, double *psam);
void write_meta(FILE * fp, struct_seed * seed, char *measfile, double p_value);
void print_psam(struct_seed * seed, char *measfile, char *outdir, double p_value,
                long motif_num);
void get_tag_string_pair(char *tag, char *btag, char *etag);
void extract_xml_line_string(FILE * fp, char *tag, char *strval, long required);
long extract_xml_line_long(FILE * fp, char *tag, long dft_value, long required);
double extract_xml_line_double(FILE * fp, char *tag, double dft_value, long required);
void extract_xml_psam(FILE * fp, long npos, double *psam);
void take_expanded_motif_as_topology(char *motif, char *topology);
void set_default_topology(char *topology, long num);
void check_psam_value(long num, double *psam);
void check_psam_topology(long num, char *topology);
void check_psam_optimal_seq(long num, double *psam, char *optseq);
void write_pwm_entry(char *pwm_file, struct_psam * pd);
void write_psam_struct(char *psam_file, struct_psam * pd);
void read_psam_struct(char *psam_file, struct_psam * pd);
void write_psam_struct_lite(char *psam_file, struct_psam * pd);
void read_psam_struct_lite(char *psam_file, struct_psam * pd);
void populate_w_from_psam(double *psam, struct_seed * seed);
void expanded_motif_to_seedtopo(char *motif, struct_seed * seed);
void expanded_motif_to_seed(char *motif, long strand, struct_seed * seed);
void get_motif_counts(char *motif, long strand, struct_seqArr * seqs, double *counts);
void motif2psam_struct(char *motif, long strand, double p_value, struct_psam * psam_data);
void recover_seed_from_psam(struct_psam * pd, struct_seed * seed);
void populate_motif_psam(struct_seed * seed);
void rc_psam(struct_psam * psam_data);
void iupac2psam(char *motif, double *psam);
long is_valid_iupac_motif(char *motif, long debug);
void normalize_psam(long ma, double *a, double *psam);
void normalize_each_psam_position(long num, double *psam, long debug);

/* fncs_seed.c */
long is_motif_to_skip(long idx, double *sumX, char *cklist, long dyn_stnd,
                      long strand, long *rc_lookup);
long is_nonbit_motif_to_skip(long sumX, char *motif, long dyn_stnd, long strand);
void populate_mcount(struct_topo * t, struct_seq * seq, long *mcount);
void get_mcount_per_seq(struct_topo * t, double dval, struct_seq * seq, long msize,
                        long *mcount, long *rc_lookup, long strand, double *sumX,
                        double *sumX2, double *sumXY);
void init_fitpars(struct_ls5fit * fitpars);
void init_seed_motif(struct_seed * seed, long strand, double p_value);
long *populate_rc_lookup(long msize, long Xcount);
void populate_seed_fitpars(struct_data * tdat, struct_seqArr * seqs,
                           struct_topo * topos, struct_seed * seed);
void check_motif_for_iupac_degeneracy(args_reduce * args, struct_seed * seed,
                                      struct_data * tdat, struct_seqArr * seqs);
void deduce_motif_full(struct_topo * t, long idx, char *motif, char *full);
void revert_full_to_motif(struct_seed * seed);
void fulfill_smotif(struct_topo * topos, char **col_names, struct_seed * seed);
void output_linefit_pars(FILE * fp, struct_ls5fit * fitpars);
void print_smotif(FILE * fp, struct_seed * seed, char *stnd);
void print_seed_expt(FILE * fp, struct_seed * seed, char *stnd);
void setup_seed_from_topo_motif(struct_data * tdat, struct_seqArr * seqs,
                                struct_seed * seed, struct_ntop * nte,
                                struct_topo * topos, struct_dic * mlist);
void print_ntop_entry(FILE * fp, long idx, struct_seed * seed, long maxlen);
void check_ntops_populate_seed(struct_ntopArr * ntop_seeds, struct_topo * topos,
                               long strand0, struct_data * tdat, struct_seqArr * seqs,
                               struct_seed * seed, struct_dic * mlist);
void initialize_struct_ntopArr(long ntop, long nsel, struct_ntopArr * ntop_seeds);
void find_seed_by_topos(struct_data * tdat, struct_seqArr * seqs, struct_topoArr * topo,
                        struct_seed * seed, long ntop, long nsel);
void find_seed_by_dictfile(struct_data * tdat, struct_seqArr * seqs,
                           struct_topoArr * topo, struct_dic * mlist, struct_seed * seed,
                           long ntop, long nsel);
void find_seed_experiment(struct_data * tdat, struct_seqArr * seqs, struct_seed * seed);
long check_dictfile_repeat(long num, struct_tag * str_tags);
void validate_expanded_pattern(char *str);
void tidy_dictfile(char *dictfile, char *cln_dictfile);
void read_dictfile(char *dictfile, struct_dic * mlist);
void get_topo_from_topology(char *topology, struct_topo * t);
void categorize_mlist(struct_topoArr * topo, struct_dic * mlist);

/* fncs_seqs.c */
void read_sequences(char *seqfile, struct_seqArr * seqs, long raw_seq, long sort_seq);
void extract_id_header(long num, char *line, struct_seq * seq);
void assign_sequence(long num, char *pb, long nb, struct_seq * sequences, long raw_seq);
void write_oneline_sequences(char *midx, struct_seqArr * seqs, char *filename);
void write_hjb_sequences(char *midx, struct_seqArr * seqs, char *filename);
void format_fasta_sequence(FILE * fp, long nb, char *seq);
void write_one_seq(FILE * fp, struct_seq * seq);
void write_sequences(char *midx, struct_seqArr * seqs, char *filename);
void write_sequences_in_tdat_order(char *filename, struct_seqArr * seqs,
                                   struct_data * tdat);
void convert_sequences_acgt_to_0123(struct_seqArr * seqs);
void reverse_cmpl_sequences(struct_seqArr * seqs);
void reverse_sequences(struct_seqArr * seqs);
long get_number_of_fasta_entries(char *seqfile);
void free_seq_entry(struct_seq * seq);
void free_sequences(char *fidx, struct_seqArr * seqs);
char **allocate_seq_fragments(long num);
struct_seq *allocate_memory_for_sequences(long num);
long base_offset(int base);
void seq2number(char *p);
void number2seq(char *p, long len);
char cvt_base_U2T(char base);
char cvt_base_T2U(char base);
void convert_U_to_T(char *seq);
void convert_T_to_U(char *seq);
int seqid_compare(const void *v1, const void *v2);
void check_redundancy(char *seqfile, long num, struct_seq * sequences);
void match_ids_seqs(long num, struct_tag * str_tags, struct_seqArr * seqs, char *midx);
long get_tdat_ID_match_number(struct_data * tdat);
void get_matched_seq_tdat(struct_seqArr * seqs, struct_data * tdat);
void get_seq_tdat(char *outdir, char *seqfile, struct_seqArr * seqs, char *measfile,
                  struct_data * tdat);
void change_seq_case(struct_seqArr * seqs, long base_case);
void read_windows(char *winfile, struct_winArr * wins, long width);
void initialize_struct_win(struct_win * w);
void populate_windows(char *winfile, struct_win * windows, long width);
void write_windows(char *outdir, char *filename, struct_seqArr * seqs,
                   struct_winArr * wins);
void print_struct_win(FILE * fp, struct_win * w, struct_seq * np);
void free_logo_sequences(long num_seqs, struct_seq * sequences);
void assign_logo_sequence(long num, char *pb, long nb, struct_seq * sequences);
struct_seq *read_logo_fasta(char *seqfile, long *num_seqs, char *outdir);
void check_logo_sequences(long num_seqs, struct_seq * sequences, char *seqfile,
                          char *outdir);
struct_seq *read_logo_flat(char *seqfile, long *num_seqs, char *outdir);

/* fncs_tdat.c */
void read_tdt_data(char *tdtfile, struct_data * tdat);
void free_tdt_data(struct_data * tdat);
void get_tdat_col_info(char *tdtfile, long *header_num, long *header_type, long *ncol);
void get_tdat_col_names(char *tdtfile, long header_num, long header_type, long ncol,
                        char **col_names);
void get_tdat_nrow(char *tdtfile, long header_num, long ncol, long *nrow);
void read_tdat(char *tdtfile, long header_num, char **ids, double **data);
void write_dval(FILE * fp, double dval);
void write_tdat_row(FILE * fp, long row_idx, struct_data * tdat, double **dmtx);
void write_tdat_transpose(struct_data * tdat, char *filename, double **dmtx,
                          struct_id * my_ids);
void write_tdat(struct_data * tdat, char *filename, double **dmtx, struct_id * my_ids);
struct_id *get_sorted_ids(long num_ids, char **ids, char *filename);
void match_ids_tdat(long num, struct_tag * str_tags, long nrow, struct_id * my_ids,
                    long *midx);
void mask_off_na_rows(char *fitmtx_file, struct_data * tdat);
void write_problematic_ids(char *outdir, char *tdtfile, long nrow, char **ids,
                           long num_missing_seq, long num_duplicate, long *check_idx);
void link_expr2seq(struct_data * tdat, struct_seqArr * seqs, char *outdir, char *tdtfile);
void change_id_case(struct_data * tdat, long base_case);
void log_transform_tdat(struct_data * tdat, char *log_base);
void get_tdat_oknum(struct_data * tdat);
long get_number_of_valid_rows(long nvec, long nlen, double **Xmtx);
void check_columns(struct_data * tdat, char *columns);

/* fncs_topo.c */
void read_topos(char *topofile, long list, struct_topoArr * topo);
void unify_parens_chars(char *str);
void expand_abbr_topos(long *num, struct_tag * str_tags);
long expand_char_num_to_full(char *src, char *dst);
long char_count_in_string(char *str, char c);
void populate_topos(long num_topo, struct_topo * topos, struct_tag * str_tags);
void write_topos(struct_topoArr * topo, char *outdir, char *topo_file);
struct_topo *allocate_memory_for_topos(long num);
void free_topos(struct_topoArr * topo);
struct_topo *allocate_memory_for_one_topo(void);
void free_one_topo(struct_topo * topo);
void duplicate_topo(struct_topo * src, struct_topo * dst);
long has_matched_parens(char *str);
long fillup_parens(char *str, long limit, long *parens_bidx, long *parens_eidx);
void match_parens_str(long num_parens, long *parens_bidx, long *parens_eidx, char *str);
long with_matched_parens_bin(struct_topo * t, char *b);

/* fncs_util.c */
void nrerror(char *error_text);
void vector_boundary_check(long nl, long nh, char *fun_name);
void matrix_boundary_check(long nrl, long nrh, long ncl, long nch, char *fun_name);
void init_cvector(char *v, long nl, long nh, char c);
void init_cmatrix(char **m, long nrl, long nrh, long ncl, long nch, char c);
char *cvector(long nl, long nh);
void free_cvector(char *v, long nl, long nh);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void diff_dvector(double *df, double *d1, double *d2, long nl, long nh);
double ave_dvector(double *d, long nl, long nh);
double var_dvector(double *d, long nl, long nh);
double max_dvector(double *d, long nl, long nh);
double min_dvector(double *d, long nl, long nh);
double norm_dvector(double *d, long nl, long nh);
double norm_inf_dvector(double *d, long nl, long nh);
void copy_dvector(double *src, long nl, long nh, double *dst);
void init_dvector(double *v, long nl, long nh, double d);
void init_dmatrix(double **m, long nrl, long nrh, long ncl, long nch, double d);
void identity_dmatrix(double **m, long n);
void copy_dmatrix(double **a, long nrl, long nrh, long ncl, long nch, double **o);
void copy_lvector(long *src, long nl, long nh, long *dst);
void init_lvector(long *v, long nrl, long nrh, long l);
void init_lmatrix(long **m, long nrl, long nrh, long ncl, long nch, long l);
long *lvector(long nl, long nh);
void free_lvector(long *v, long nl, long nh);
long **lmatrix(long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
double dval_sqr(double dval);
double dval_max(double a, double b);
double dval_min(double a, double b);
long lval_max(long a, long b);
long lval_min(long a, long b);
long dval_in_range(double dval, double dlow, double dhigh);
long lval_in_range(long lval, long llow, long lhigh);
long lval_in_array(long lval, long ib, long ie, long *lvec);
void dval_swap(double *pa, double *pb);
void lval_swap(long *pa, long *pb);
void cval_swap(char *pa, char *pb);
void str_swap(char **src, char **dst);
void dvec_swap(double **src, double **dst);
int dval_compare(const void *v1, const void *v2);
int lval_compare(const void *v1, const void *v2);
int cstr_compare(const void *v1, const void *v2);
int strtags_compare(const void *v1, const void *v2);
int idstr_compare(const void *v1, const void *v2);
int ntop_compare(const void *v1, const void *v2);
void call_system(char *cmd);
long decompose_list_ids_to_file(char *ids_str, char *sep_chars, char *filename);
long extract_numlist(char *num_str, char *sep_chars, long lmin, long lmax, long *vnum);
