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
#include <R.h>

long is_palindrome(char *str)
{
    char *p;
    long k;

    p = my_strdup(str);
    reverse_string(p);

    k = is_equal_string(p, str);

    free_cvector(p, 0, DUMMY);

    return k;
}

long is_rc_palindrome(char *str)
{
    char *p;
    long k;

    p = my_strdup(str);
    reverse_cmpl(p, FALSE);

    k = is_equal_string(p, str);

    free_cvector(p, 0, DUMMY);

    return k;
}

void reverse_string(char *str)
/* reverse string 'str' in place */
{
    long i = 0, j = strlen(str) - 1;

    while (i < j) {
        cval_swap(&str[i], &str[j]);
        i++;
        j--;
    }
}

long string_contains_only_those_characters(char *str, char *chars_set)
{
    return strspn(str, chars_set) == strlen(str);
}

static long isEMPTY_or_isNAN_or_isINF(char *str)
{
    lowerstr(str);  /* to lower case for consistency: NAN/nan etc */

    /* empty field/missing value etc: \t\t [trim()]; NaN; Inf/-Inf */
    return *str == '\0' || strstr(str, "nan") || strstr(str, "inf");
}

long is_numeric(char *str)
/* check if string "str" contains only valid numerical values */
{
    char *endp, *p = trim(str);
    double d;

    if (isEMPTY_or_isNAN_or_isINF(p))
        return FALSE;  /* special cases: not taken as a number */

    errno = 0;
    d = strtod(p, &endp);
    UNUSED_PARAMETER(d);

    return !(*endp != '\0' || errno == ERANGE);
}

double cvt2double(char *str)
/* convert string "str" to a double value, with error checking */
{
    char *endp, *p = trim(str);
    double d;

    if (isEMPTY_or_isNAN_or_isINF(p))
        return XBIG;  /* as a practically impossible value */

    errno = 0;
    d = strtod(p, &endp);

    if (*endp != '\0' || errno == ERANGE)
        return XBIG;  /* for NA etc */
    else
        return d;
}

long cvt2long(char *str)
/* convert string "str" to a long value, with error checking */
{
    char *endp, *p = trim(str);
    long d;

    if (isEMPTY_or_isNAN_or_isINF(p))
        return LONG_MAX;  /* as a practically impossible value */

    errno = 0;
    d = strtol(p, &endp, 10);

    if (*endp != '\0' || errno == ERANGE)
        return LONG_MAX;  /* for NA etc */
    else
        return d;
}

void quit_with_help_reminder(char *str)
{
    fprintf(stderr, "%s\n", str);
    fatal("Please use \"%s -h\" for usage information\n", Gvars.PROGNAME);
}

long has_no_equal_sign(char *parstr)
{
    return strchr(parstr, '=') == NULL;
}

static long equalsign_pos(char *parstr)
/* specifically for checking par=val pair: return pos 1-char after '=' */
{
    char *ptr, str[BUF512];

    ptr = strchr(parstr, '=');  /* check from the beginning */

    if (ptr == NULL) {
        sprintf(str, "Please specify option as %s=val\n", parstr);
        quit_with_help_reminder(str);
    }

    return ptr - parstr + 1;
}

long get_lvalue(char *str, long vmin, long vmax)
/* extract long numerical value from a command line option */
{
    char *p;
    long npos, val;

    npos = equalsign_pos(str);
    p = strtok(str, WSPACES);  /* to handle '    -dDEVICEWIDTHPOINTS=340 \' */
    UNUSED_PARAMETER(p);

    val = cvt2long(str + npos);
    if (val == LONG_MAX)
        fatal("wrong option [%s]: not a valid integer value\n", str);

    if (!lval_in_range(val, vmin, vmax))
        fatal("invalid option [%s]: value %ld out of range [%ld %ld]\n",
              str, val, vmin, vmax);

    return val;
}

double get_dvalue(char *str, double vmin, double vmax)
/* extract double numerical value from a command line option */
{
    char *p;
    long npos;
    double val;

    npos = equalsign_pos(str);
    p = strtok(str, WSPACES);
    UNUSED_PARAMETER(p);

    val = cvt2double(str + npos);
    if (val > XBIG_CUTOFF)
        fatal("wrong option [%s]: not a valid numerical value\n", str);

    if (!dval_in_range(val, vmin, vmax))
        fatal("invalid option [%s]: value %f out of range [%g %g]\n",
              str, val, vmin, vmax);

    return val;
}

void get_strvalue(char *str, char *dst, long expand_tilde)
/* extract string from a command line option, substitute ~ to $HOME */
{
    long npos;

    npos = equalsign_pos(str);

    if (expand_tilde && str[npos] == '~') {
        char *p = getenv("HOME");

        if (p == NULL) {
            log_msg("environment variable HOME NOT defined!");
            strcpy(dst, str + npos);

        } else  /* HOME does not have ending slash */
            sprintf(dst, "%s%s", p, str + npos + 1);  /* +1 to skip '~' */

    } else
        strcpy(dst, str + npos);
}

long str_pmatch(char *str, char *sstr)
/* partial command-line option match, from the beginning */
{
    return !strncmp(str, sstr, strlen(sstr));
}

long case_str_pmatch(char *str, char *sstr)
/* case-insensitive command-line option match, from the beginning */
{
    return !case_strncmp(str, sstr, strlen(sstr));
}

void print_used_time(time_t time0)
{
    char msg[BUF512];
    double dtime;
    long minute_secs = 60, hour_secs = 60 * 60, day_secs = 24 * 60 * 60;
    long days, hours, minutes, seconds;

    dtime = difftime(time(NULL), time0);
    seconds = lround(dtime);  /* rounded to seconds in long */

    days = seconds / day_secs;
    seconds %= day_secs;

    hours = seconds / hour_secs;
    seconds %= hour_secs;

    minutes = seconds / minute_secs;
    seconds %= minute_secs;

    sprintf(msg, "\nTime used: %2.2ld:%2.2ld:%2.2ld:%2.2ld", days, hours, minutes,
            seconds);
    log_msg(msg);
}

void get_currentTimeString(char *s)
{
    time_t t0;

    t0 = time(NULL);
    strcpy(s, ctime(&t0));
    s[strlen(s) - 1] = '\0';  /* change \n to \0 */
}

void get_currentYear(char *year)
{
    time_t t0;

    t0 = time(NULL);
    strftime(year, 5, "%Y", localtime(&t0));  /* 5 = 4 + 1 [for '\0'] */
}

double my_log2(double d)
/* name 'log2' is not in ANSI C, but commonly provided, e.g. in GCC */
{
    if (d <= 0.0)
        fatal("d=%g <= 0.0 in my_log2()\n", d);

    return log(d) / log(2);
}

char *enlarge_cline(long *maxline, char *line)
/* for a 0-index char-array */
{
    char *newline;

    *maxline *= SFACTOR;  /* new enlarged size */

    if ((newline = (char *) realloc(line, (*maxline) * sizeof(char))) == NULL)
        fatal("realloc failure in enlarge_cline()\n");

    return newline;
}

static int endofline(FILE * fp, int c)
/* check for and consume \r, \n, \r\n, or EOF */
{
    int eol;

    eol = (c == '\r' || c == '\n');

    if (c == '\r') {
        c = getc(fp);  /* read in one more character */
        if (c != '\n' && c != EOF)
            ungetc(c, fp);  /* read too far; put c back */
    }

    return eol;
}

char *my_getline(FILE * fp)
/* read line of arbitrary length from file pointer *fp */
{
    int c;
    char *line;
    long i, maxline = BUF512;

    line = (char *) malloc(maxline * sizeof(char));
    if (line == NULL)
        fatal("out of memory in my_getline()\n");

    for (i = 0; (c = getc(fp)) != EOF && !endofline(fp, c); i++) {
        if (i >= maxline - 1)  /* note '>=' here */
            line = enlarge_cline(&maxline, line);
        line[i] = c;
    }

    line[i] = '\0';

    if (c == EOF && i == 0) {
        free(line);
        line = NULL;
    }

    return line;
}

long readline_cvt2long(FILE * fp)
{
    char *p;
    long k;

    p = my_getline(fp);
    k = cvt2long(p);

    free(p);

    return k;
}

void cvtstr_set1toc2(char *str, char *set1, char c2)
/* convert all occurrences of any character in set1 to 'c2' in string 'str' */
{
    char *p = str;

    while (*p) {
        if (strchr(set1, *p))
            *p = c2;
        p++;
    }
}

void cvtstr_c1toc2(char *str, char c1, char c2)
/* convert all occurrences of character 'c1' to 'c2' in string 'str' */
{
    char *p = str;

    while (*p) {
        if (*p == c1)
            *p = c2;
        p++;
    }
}

long csplit(char *str, char *items[], long itemsize, char sepc)
/* strtok() is GREEDY where consecutive empty fields will be consumed:
 * e.g. \t\t. This function splits the string at each char 'spec', and
 * return an 1-index pointer array items[], with maximum 'itemsize'. */
{
    char *p0, *p;
    long nitem = 0;

    if (str[0] == '\0')  /* empty string */
        return nitem;

    cvtstr_c1toc2(str, '"', ' ');  /* change '"' to ' ' */

    p0 = str;
    while ((p = strchr(p0, sepc)) != NULL) {
        *p = '\0';  /* change separator character to '\0' */
        items[++nitem] = trim(p0);
        if (nitem >= itemsize)
            return nitem;
        p0 = p + 1;  /* for next field */
    }
    items[++nitem] = trim(p0);  /* last field */

    return nitem;
}

long item_list(char *str, char *items[], long itemsize, char *sep_chars)
/* itemize a string into tokens; return no. of token fields: 1-index */
{
    char *p;
    long nitem = 0;

    for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars)) {
        items[++nitem] = trim(p);
        if (nitem >= itemsize)
            return nitem;
    }

    return nitem;
}

FILE *open_tmpfile(void)
{
    FILE *fp;

    errno = 0;
    fp = tmpfile();
    if (fp == NULL)
        fatal("open_tmpfile() failed: %s\n", strerror(errno));

    return fp;
}

long is_std_out_err(char *filename)
{
    return is_equal_string(filename, "stdout") || is_equal_string(filename, "stderr");
}

FILE *open_file(char *filename, char *filemode)
{
    FILE *fp;

    if (filename == NULL)
        *filename = '\0';  /* filename = "\0": p347 C-ARM 4th-edition */

    if (is_equal_string(filename, "stdin"))
        fp = stdin;

    else if (is_equal_string(filename, "stdout"))
        fp = stdout;

    else if (is_equal_string(filename, "stderr"))
        fp = stderr;

    else {
        errno = 0;
        fp = fopen(filename, filemode);  /* filename can't contains ~ (e.g., ~/fname) */
        if (fp == NULL)
            fatal("open_file() <%s> failed: %s\n", filename, strerror(errno));
    }

    return fp;
}

long close_file(FILE * fp)
{
    long i;

    if (fp == NULL || fp == stdin || fp == stdout || fp == stderr)
        return 0;

    errno = 0;
    i = fclose(fp);
    if (i == EOF)
        fatal("close_file() failed: %s\n", strerror(errno));

    return i;
}

long exist_file(char *filename)
{
    long k;
    FILE *fp;

    fp = fopen(filename, "r");  /* NOT open_file() */
    k = (fp != NULL);
    close_file(fp);

    return k;
}

void remove_file(char *filename)
{
    if (!exist_file(filename))
        return;

    errno = 0;
    if (remove(filename))
        fatal("can not remove file: %s <%s>\n", filename, strerror(errno));
}

long is_same_file(char *src, char *dst)
{
    char msg[BUF512];

    if (is_equal_string(src, dst)) {
        sprintf(msg, "same source/destination file: <%s>", src);
        log_msg(msg);

        return TRUE;
    }

    return FALSE;
}

void rename_file(char *src, char *dst)
{
    if (!exist_file(src))
        fatal("file to be renamed <%s> does NOT exist\n", src);

    if (is_equal_string(src, dst))
        return;

    errno = 0;
    if (rename(src, dst))
        fatal("can not rename file <%s> to <%s>: <%s>\n", src, dst, strerror(errno));
}

void copy_file(char *src, char *dst)
{
    char str[BUF512];
    size_t num_bytes;
    FILE *fpi, *fpo;

    if (is_equal_string(src, dst))
        return;

    fpi = open_file(src, "rb");
    fpo = open_file(dst, "wb");

    while (!feof(fpi)) {
        num_bytes = fread(str, 1, BUF510, fpi);
        if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
            fatal("file copy error from <%s> to <%s>\n", src, dst);
    }

    close_file(fpi);
    close_file(fpo);
}

void cat_file(char *src, char *dst)
/* cat file contents from src to dst: similar to copy_file */
{
    char str[BUF512];
    size_t num_bytes;
    FILE *fpi, *fpo;

    fpi = open_file(src, "rb");
    fpo = open_file(dst, "ab");

    while (!feof(fpi)) {
        num_bytes = fread(str, 1, BUF510, fpi);
        if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
            fatal("file concatenation error from <%s> to <%s>\n", src, dst);
    }

    close_file(fpi);
    close_file(fpo);
}

void get_alnum(char *a, char c)
/* convert chars other than alnum & 'c' to char 'c' */
{
    char *str;
    long i, j = 0, nchar = strlen(a);

    if (isalnum((int) c))
        fatal("char '%c' is already an alpha-numerical character!\n", c);

    str = cvector(0, nchar);

    for (i = 0; i < nchar; i++) {
        if (isalnum((int) a[i]) || a[i] == c)
            str[j++] = a[i];
        else if (j >= 1 && str[j - 1] != c)  /* to avoid two 'c' in a row */
            str[j++] = c;
    }

    str[j] = '\0';
    strcpy(a, str);  /* copy back to a */

    free_cvector(str, 0, DUMMY);
}

char *ltrim(char *a)
/* trim leading (left) white spaces */
{
    int c;

    while (isspace(c = *a))
        a++;

    return a;
}

char *rtrim(char *a)
/* trim trailing (right) white spaces */
{
    long i;

    for (i = strlen(a) - 1; i >= 0; i--)
        if (!isspace((int) a[i]))
            break;

    a[i + 1] = '\0';

    return a;
}

char *trim(char *a)
/* trim leading and trailing white spaces */
{
    return rtrim(ltrim(a));
}

void upperstr(char *a)
{
    while (*a) {
        if (islower((int) *a))
            *a = toupper((int) *a);
        a++;
    }
}

void lowerstr(char *a)
{
    while (*a) {
        if (isupper((int) *a))
            *a = tolower((int) *a);
        a++;
    }
}

int case_strcmp(const char *s1, const char *s2)
{
    int i, c1, c2;

    for (i = 0; c1 = toupper((int) s1[i]), c2 = toupper((int) s2[i]), c1 == c2; i++)
        if (c1 == '\0')
            return 0;

    return c1 - c2;
}

int case_strncmp(const char *s1, const char *s2, long n)
{
    int i, c1, c2;

    for (i = 0; (c1 = toupper((int) s1[i]), c2 = toupper((int) s2[i]), c1 == c2) && i < n;
         i++)
        if (c1 == '\0')
            return 0;

    return (i >= n) ? 0 : c1 - c2;
}

char *case_strstr(const char *haystack, const char *needle)
{
    char *haystack_cp, *needle_cp, *p;

    haystack_cp = my_strdup(haystack);
    needle_cp = my_strdup(needle);

    lowerstr(haystack_cp);
    lowerstr(needle_cp);

    p = strstr(haystack_cp, needle_cp);

    if (p != NULL)  /* make p a pointer w.r.t. the original haystack */
        p = (p - haystack_cp) + (char *) haystack;

    free_cvector(haystack_cp, 0, DUMMY);
    free_cvector(needle_cp, 0, DUMMY);

    return p;
}

char *case_strchr(const char *s, int c)
{
    char *str, *p;

    str = my_strdup(s);
    lowerstr(str);

    p = strchr(str, tolower(c));

    if (p != NULL)
        p = (p - str) + (char *) s;

    free_cvector(str, 0, DUMMY);

    return p;
}

char *my_strdup(const char *src)
/* memory allocated here needs to be freed elsewhere */
{
    char *dst;

    dst = cvector(0, strlen(src));
    strcpy(dst, src);

    return dst;
}

void print_sep(FILE * fp, char c, long num)
/* print char 'c' num-times to stream fp, with '\n' added */
{
    long i;

    for (i = 1; i <= num; i++)
        if (fputc(c, fp) == EOF)
            fatal("error writing character <%c> to the stream\n", c);

    if (fputc('\n', fp) == EOF)
        fatal("error writing '\n' to the stream\n");
}

void repeat_char_ntimes_string(char c, long n, char *str)
{
    long i;

    for (i = 0; i < n; i++)
        str[i] = c;

    str[i] = '\0';  /* here: i = n */
}

long num_strmatch(char *str, char **strmat, long ib, long ie)
/*  return number of matchs of 'str' in 'strmat' */
{
    long i, num = 0;

    for (i = ib; i <= ie; i++)
        if (is_equal_string(str, strmat[i]))
            num++;

    return num;
}

void add_end_slash(char *str)
{
    long n = strlen(str);

    if (str[n - 1] != '/')
        strcat(str, "/");
}

void delete_end_slash(char *str)
{
    long n1 = strlen(str) - 1;

    if (str[n1] == '/')
        str[n1] = '\0';
}

char *basename(char *str)
/* return a pointer 1-char following the last slash or to str w/o '/' */
{
    char *p, *os;

    p = strrchr(str, '\\');  /* handle MinGW: \ as separator and with OS env */
    os = getenv("OS");
    if (p != NULL && os != NULL && case_strstr(os, "WINDOWS") != NULL)
        return p + 1;

    p = strrchr(str, '/');

    return (p == NULL) ? str : p + 1;
}

void del_extension(char *fullname, char *okname)
/* get rid of the extension in a file name, taking into consideration
 * of leading . and .. */
{
    char *p;
    size_t i;

    p = strrchr(fullname, '.');

    if (p == NULL)
        strcpy(okname, fullname);

    else {
        i = p - fullname;
        if (i == 0 || (i == 1 && fullname[0] == '.'))  /* leading '.' or '..' */
            strcpy(okname, fullname);

        else {
            strncpy(okname, fullname, i);
            okname[i] = '\0';  /* overwrite '.' at the end */
        }
    }
}

void bname_noext(char *src, char *dst)
/* extract base name, without extension */
{
    del_extension(basename(src), dst);
}

void bname_ext(char *src, char *ext, char *dst)
/* get a new name with base name from 'src' and extension 'ext' */
{
    char str[BUF512];

    bname_noext(src, str);
    sprintf(dst, "%s.%s", str, ext);
}

void get_stropt_wo_end_slash(char *option, char *okstr)
{
    get_strvalue(option, okstr, TRUE);
    delete_end_slash(okstr);
}

long lround(double d)
{
    return (long) ((d > 0.0) ? d + 0.5 : d - 0.5);
}

long dbl2long(double dval)
/* 5.001 ===> 6; -5.001 ===> -6 */
{
    long lval = (long) dval;

    if (lval < dval - DBL_EPSILON)
        lval++;

    else if (lval > dval + DBL_EPSILON)
        lval--;

    return lval;
}

void fatal(char *fmt, ...)
/* exit program with a customized error message */
{
    va_list args;
    char *error_message = NULL;

    if (strlen(fmt) > 0) {
        va_start(args, fmt);
        vasprintf(&error_message, fmt, args);
        va_end(args);
        error(error_message);
    }

    error("unspecified error\n");
}

void change_case_str_tags(long upper_case, long num, struct_tag * str_tags)
{
    long i;

    if (upper_case) {
        for (i = 1; i <= num; i++)
            upperstr(str_tags[i].str);

    } else {
        for (i = 1; i <= num; i++)
            lowerstr(str_tags[i].str);
    }
}

struct_tag *fillup_str_tags(long is_list, char *filename, long *num)
{
    long num0;
    struct_tag *str_tags;

    if (is_list) {
        num0 = get_line_number(filename, TRUE);
        str_tags = allocate_memory_for_strtags(num0);
        read_strtags(filename, str_tags);

    } else {
        num0 = 1;
        str_tags = allocate_memory_for_strtags(num0);
        set_strtag(filename, &str_tags[num0]);
    }

    *num = num0;

    return str_tags;
}

struct_tag *read_unique_strtags(char *idfile, long *num)
{
    struct_tag *str_tags;

    str_tags = fillup_str_tags(TRUE, idfile, num);
    get_unique_strtag(num, str_tags);

    if (*num == 0)
        fatal("no valid IDs from file [%s]\n", idfile);

    return str_tags;
}

void init_strtag_with_strings(char *str, char *tag, struct_tag * s)
{
    s->str = my_strdup(str);

    if (tag == NULL)
        s->tag = my_strdup(str);
    else
        s->tag = my_strdup(tag);
}

void init_strtag(struct_tag * src, struct_tag * dst)
{
    dst->str = my_strdup(src->str);
    dst->tag = my_strdup(src->tag);
}

void copy_strtag(struct_tag * src, struct_tag * dst)
{
    free_p_strtag(dst);  /* to avoid memory leak */
    init_strtag(src, dst);
}

void swap_strtag(struct_tag * one, struct_tag * two)
{
    struct_tag tmp;

    if (one == two)
        return;

    init_strtag(one, &tmp);
    copy_strtag(two, one);
    copy_strtag(&tmp, two);

    free_p_strtag(&tmp);
}

void get_unique_strtag(long *num, struct_tag * str_tags)
{
    char msg[BUF512];
    long i, j, *masked, num_ok = 0;

    if (*num == 0)
        return;

    log_prg("    getting unique entries...");

    masked = lvector(1, *num);

    for (i = 1; i <= *num; i++)  /* ignore empty strings */
        if (is_empty_string(str_tags[i].str))
            masked[i] = TRUE;

    for (i = 1; i <= *num - 1; i++) {
        if (masked[i])
            continue;

        for (j = i + 1; j <= *num; j++) {
            if (masked[j])
                continue;

            if (is_equal_case_string(str_tags[i].str, str_tags[j].str)) {
                sprintf(msg, "\t***IGNORE string [%s (%s %ld)] -- repeat of [%s %ld]\n",
                        str_tags[j].str, str_tags[j].tag, j, str_tags[i].tag, i);
                log_prg(msg);
                masked[j] = TRUE;  /* masked off later entries */
            }
        }
    }

    for (i = 1; i <= *num; i++) {
        if (masked[i])
            continue;

        num_ok++;
        if (i != num_ok)
            copy_strtag(&str_tags[i], &str_tags[num_ok]);
    }

    free_lvector(masked, 1, DUMMY);

    free_extra_strtags(num_ok + 1, *num, str_tags);

    *num = num_ok;
}

void read_strtags(char *filename, struct_tag * str_tags)
{
    char *line, *p0;
    long num = 0;
    FILE *fp;

    fp = open_file(filename, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        line = trim(p0);  /* keep the original value of p0 */
        if (!is_skip_line(line)) {
            num++;
            set_strtag(line, &str_tags[num]);
        }

        free(p0);
    }

    close_file(fp);
}

void set_strtag(char *line, struct_tag * str_tag)
/* populate one str-tag entry */
{
    char *p_trim, *items[3];
    long nitem, TWO = 2;  /* maximum two items */

    nullify_line_comment(line);

    p_trim = trim(line);
    nitem = item_list(p_trim, items, TWO, WSPACES);

    if (nitem == TWO) {
        str_tag->str = my_strdup(items[1]);  /* 'str' (first item) is primary */
        str_tag->tag = my_strdup(items[2]);

    } else if (nitem == TWO - 1) {
        str_tag->str = my_strdup(items[1]);
        str_tag->tag = my_strdup(items[1]);

    } else  /* unlikely to happen */
        fatal("not a valid input <%s>\n", line);
}

void skip_lines(FILE * fp, long num)
{
    char *p0;
    long i;

    for (i = 1; i <= num; i++) {
        if ((p0 = my_getline(fp)) == NULL)
            fatal("file contains less than required <%ld> lines to skip\n", num);

        free(p0);
    }
}

long get_line_number(char *filename, long skips)
{
    char *p0, *line;
    long num = 0;
    FILE *fp;

    fp = open_file(filename, "r");

    while ((p0 = my_getline(fp)) != NULL) {
        if (skips) {
            line = trim(p0);  /* keep the original address of p0 */
            if (is_skip_line(line)) {
                free(p0);
                continue;
            }
        }

        num++;
        free(p0);
    }

    close_file(fp);

    if (!num)
        fatal("File <%s> contains 0 valid records\n", filename);

    return num;
}

long is_empty_string(const char *str)
{
    return (str == NULL || strcmp(str, "") == 0);
}

long is_equal_string(const char *str1, const char *str2)
{
    return (strcmp(str1, str2) == 0);
}

long is_equal_case_string(const char *str1, const char *str2)
{
    return (case_strcmp(str1, str2) == 0);
}

void assert_equal_longs(char *msg, long a, long b)
{
    if (a != b)
        fatal("%s: %ld =/= %ld\n", msg, a, b);
}

void set_parallel_idx(long ib, long ie, long *idx)
{
    long i;

    for (i = ib; i <= ie; i++)
        idx[i] = i;
}

void nullify_line_comment(char *str)
{
    char *p;

    p = strrchr(str, '#');
    if (p) {
        if (!isspace((int) *(p - 1)))
            fatal("wrong format [%s]: requiring a white space before #\n", str);
        *p = '\0';
    }
}

long is_comment_line(char *line)
{
    char *p = ltrim(line);

    return (*p == '#');
}

/* line should be trim()-ed by the caller */
long is_skip_line(char *line)
{
    return (strchr(SKIPS, *line) != NULL);  /* line starting with # or empty */
}

long is_fasta_skip_line(char *line)
{
    return (strchr(FASTA_SKIPS, *line) != NULL);  /* line starting with #/; or empty */
}

void transpose_matrix(double **src, long nr, long nc, double **dst)
{
    long i, j;

    for (i = 1; i <= nc; i++)
        for (j = 1; j <= nr; j++)
            dst[i][j] = src[j][i];
}

void print_dvector(double *dvec, long nl, long nh, char *msg)
{
    long i;

    fprintf(stderr, "%s\n", msg);

    for (i = nl; i <= nh; i++)
        fprintf(stderr, "%ld\t%20.12f\n", i, dvec[i]);
}

void do_nothing(void)
{
    return;
}
