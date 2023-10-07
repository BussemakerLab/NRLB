/*
Copyright (c) 2014, Andrew Delong.
All rights reserved.

BSD 3-clause license:
Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:

  1. Redistributions of source code must retain the
     above copyright notice, this list of conditions
     and the following disclaimer.

  2. Redistributions in binary form must reproduce
     the above copyright notice, this list of conditions
     and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the
     names of its contributors may be used to endorse or
     promote products derived from this software without
     specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#ifdef _WIN32
#include <windows.h>
#define PATH_SEP '\\'
#else
#include <unistd.h>
#include <sys/stat.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#define _stricmp strcasecmp
#define _strnicmp strncasecmp
#define PATH_SEP '/'
#endif


/* Print a formatted error message and quit the program */
void panic(char* msg, ...)
{
	va_list va;
	va_start(va, msg);
	vprintf(msg, va);
	va_end(va);
	printf("\n");
	exit(-1);
}


void trim_trailing_whitespace(char* str)
{
	int i = (int)strlen(str);
	while (i > 0 && (str[i-1] == '\n' || str[i-1] == ' ' || str[i-1] == '\t' || str[i-1] == '\r'))
		i--;
	str[i] = '\0';
}


#define INVALID_BASE -1
#define UNKNOWN_BASE 4

/* The base2index table is used to:
   Convert chars 'A','C','G','T'/'U' to integers 0,1,2,3 respectively.
   Convert char 'N' to UNKNOWN_BASE.
   Convert anything else to INVALID_BASE. */
int base2index[256];

void init_base2index_table()
{
	int i;
	for (i = 0; i < 256; ++i)
		base2index[i] = INVALID_BASE;
	base2index['a'] = base2index['A'] = 0;
	base2index['c'] = base2index['C'] = 1;
	base2index['g'] = base2index['G'] = 2;
	base2index['t'] = base2index['T'] = 3;
	base2index['u'] = base2index['U'] = 3;
	base2index['n'] = base2index['N'] = UNKNOWN_BASE;
}

/* The base2comp table is used to:
   Convert chars ACGT to their complements TGCA.
   U is converted to A, and A is always converted to T.
   N is converted to N.
*/
char base2comp[256];

void init_base2comp_table()
{
	memset(base2comp, '\0', 256);
	base2comp['a'] = base2comp['A'] = 'T';
	base2comp['c'] = base2comp['C'] = 'G';
	base2comp['g'] = base2comp['G'] = 'C';
	base2comp['t'] = base2comp['T'] = 'A';
	base2comp['u'] = base2comp['U'] = 'A';
	base2comp['n'] = base2comp['N'] = 'N';
}


char db_dir[512];
char db_params_dir[512];


void check_is_dir(char* dir_name)
{
	int found_dir = 0;
#ifdef _WIN32
	{
	DWORD attrib = GetFileAttributesA(dir_name);
	found_dir = (attrib != INVALID_FILE_ATTRIBUTES && (attrib & FILE_ATTRIBUTE_DIRECTORY));
	}
#else
	{
	struct stat st;
	found_dir = (!stat(dir_name, &st) && S_ISDIR(st.st_mode));
	}
#endif
	if (!found_dir)
		panic("Could not locate directory \"%s\"", dir_name);
}


void find_db_dirs()
{
	char db_subdir[] = {
		PATH_SEP, 'd', 'b',
		PATH_SEP,
		'\0'
	};

	char db_params_subdir[] = {
		'p', 'a', 'r', 'a', 'm', 's',
		PATH_SEP,
		'\0'
	};

	char* pos = 0;
	int exe_path_len = 0;
#ifdef _WIN32
	exe_path_len = (int)GetModuleFileNameA(NULL, db_dir, 512);
#elif defined(__APPLE__)
	exe_path_len = 512;
	if (_NSGetExecutablePath(db_dir, (uint32_t*)&exe_path_len) != 0)
		panic("Error in NSGetExecutablePath. Quitting.");
	exe_path_len = strlen(db_dir);
#else
	exe_path_len = readlink("/proc/self/exe", db_dir, 512);
#endif
	if (exe_path_len >= 512)
		panic("Path is too long for path variable! Quitting.");
	if (exe_path_len <= 0)
		panic("Error reading path of executable; cannot find deepbind database. Quitting");

	/* Get rid of \deepbind.exe suffix */
	pos = strrchr(db_dir, PATH_SEP);
	if (!pos) panic("Could not parse path name of this executable. Quitting.");
	*pos = '\0';

	/* Get rid of \Debug or \Release suffix, if any */
	pos = strrchr(db_dir, PATH_SEP);
	if (!pos) panic("Could not parse path name of this executable. Quitting.");
	if (!_stricmp(pos+1, "debug") || !_stricmp(pos+1, "release"))
		*pos = '\0';

	/* Append /db/ to the current db_dir */
	strcat(db_dir, db_subdir);
	check_is_dir(db_dir);

	/* Append params/ to the current db_params_dir */
	strcpy(db_params_dir, db_dir);
	strcat(db_params_dir, db_params_subdir);
	check_is_dir(db_params_dir);
}


int has_extension(char* file_name, char* ext)
{
	char* pos = strrchr(file_name, '.');
	if (pos && pos > file_name)
		return _stricmp(pos+1, ext) ? 0 : 1;
	return 0;
}


int is_seq_filename(char* filename)
{
	return (has_extension(filename, "fa") || has_extension(filename, "seq")) ? 1 : 0;
}


int is_fasta_comment(char* line)
{
	return (line[0] == '>' || line[0] == ';') ? 1 : 0;
}


#define MODEL_ID_LEN 10

typedef struct {
	int major;
	int minor;
} model_id_t;


model_id_t str2id(char* str)
{
	model_id_t id = { 0, 0 };
	char tmp[6];
	if (!str)
		panic("Invalid model id NULL");
	if (strlen(str) < MODEL_ID_LEN || str[0] != 'D' || str[6] != '.')
		panic("Invalid model id \"%s\"; should be of form D#####.###", str);
	memcpy(tmp, str+1, 5); tmp[5] = '\0';
	id.major = atoi(tmp);
	id.minor = atoi(str+7);
	if (id.major <= 0 || id.minor <= 0)
		panic("Invalid model id \"%s\"; should be of form DXXXXX.YYY where XXXXX >= 1 and YYY >= 1", str);
	return id;
}


void id2str(model_id_t id, char* dst)
{
	if (id.major <= 0 || id.minor <= 0)
		panic("Invalid model id");
	sprintf(dst, "D%05d.%03d", id.major, id.minor);
}


void usage()
{
	panic("%s",
	  "deepbind [--window-size n] [--average] [--echo] [--no-head] [--dump-info]\n"
	  "         <id> [<id>...] [<sequence_file>]\n"
	  "  Each <id> can be an actual ID (e.g. D00054.002), or it can be a file\n"
	  "  name with extension .ids or .txt containing a single ID on each line.\n"
	  "  \n"
	  "  deepbind reads DNA/RNA sequences either from the <sequence_file>\n"
	  "  or from stdin. If the first input character is '>', then FASTA input\n"
	  "  format is assumed. Otherwise, the input is assumed to contain just raw\n"
	  "  sequences, with one sequence on each line.\n"
	  "  Each output line contains tab-separated scores, with\n"
	  "  one score for each ID specified, in the order specified.\n"
	  "  \n"
	  "  --window-size n tells the model to use a length-n subsequence for each\n"
	  "    per-position along the larger original input sequence; the default\n"
	  "    is to use the length of the model's motif detectors (filters) x 1.5.\n"
	  "  \n"
	  "  --average causes the per-position scores to be averaged rather than maxed.\n"
	  "  \n"
	  "  --echo causes each output line echo the corresponding sequence;\n"
	  "    ideal for pasting into a spreadsheet.\n"
	  "  \n"
	  "  --no-head omits the column names (IDs) from the first line of output.\n"
	  "  \n"
	  "  --dump-info causes information about each specified ID to be dumped to\n"
	  "  stdout, e.g. the protein name, type, species, experiment, etc.\n"
	  );
}


void dump_info(model_id_t* ids, int num_ids)
{
	char dbfile_name[512];
	char line[1024];
	FILE* dbfile = 0;
	model_id_t id = { 0, 0 };
	int i;

	/* Assumes db_dir has been initialized */
	strcpy(dbfile_name, db_dir);
	strcat(dbfile_name, "db.tsv");

	for (i = 0; i < num_ids; ++i) {
		int found_id = 0;

		/* This reads the file several times, as the easiest way to dump
		   the ids in the same order they were specified.
		   This isn't something that needs to be fast, so it should be ok. */
		dbfile = fopen(dbfile_name, "r");
		if (!dbfile)
			panic("Could not open \"%s\".", dbfile_name);

		while (fgets(line, 1024, dbfile)) {
			/* find the next line starting with an id, skipping comments */
			line[1023] = '\0';
			if (line[0] == '#')
				continue;
			if (strlen(line) > 2 && line[0] == 'I' && line[1] == 'D') {
				if (i == 0)
					fputs(line, stdout);
				continue;
			}
			id = str2id(line);
			if (id.major == ids[i].major && id.minor == ids[i].minor) {
				fputs(line, stdout);
				found_id = 1;
				break;
			}
		}

		fclose(dbfile);

		if (!found_id)
			panic("Could not find id D%05d.%03d in database.", id.major, id.minor);
	}
	exit(0);
}


/* Check that a sequence is valid.
   If not, quit with error message. */
void check_valid_seq(char* seq, int seq_len)
{
	int i;
	if (seq_len == 0)
		panic("Sequence must have length > 0");

	for (i = 0; i < seq_len; ++i) {
		char c = seq[i];
		if (base2index[(unsigned char)c] == INVALID_BASE)
			panic("Invalid DNA/RNA base character '%c' in sequence \"%s\"", c, seq);
	}
}


/* All the information and coefficients needed to define a trained deepbind model. */
typedef struct {
	model_id_t id;
	int reverse_complement;
	int num_detectors;
	int detector_len;
	int has_avg_pooling;
	int num_hidden;
	float* detectors;
	float* thresholds;
	float* weights1;
	float* biases1;
	float* weights2;
	float* biases2;
} deepbind_model_t;


/* The number of inputs/outputs for the first layer
   after pooling depends on the type of pooling and
   whether there are hidden units in the next layer */
int get_num_hidden1(deepbind_model_t* model) { return model->has_avg_pooling ? model->num_detectors*2 : model->num_detectors; }
int get_num_hidden2(deepbind_model_t* model) { return model->num_hidden ? model->num_hidden : 1; }


/* Returns index a specific detector coefficient */
int indexof_detector_coeff(int num_detector, int detector, int pos, int base)
{
	assert(detector >= 0);
	assert(pos >= 0);
	assert(base >= 0 && base < 4);
	return detector + num_detector*(base + 4*pos);
}


/* Returns index a specific featuremap coefficient */
int indexof_featuremap_coeff(int num_detector, int detector, int pos)
{
	assert(detector >= 0);
	assert(pos >= 0);
	return detector + num_detector*pos;
}


/* Load a list of comma-separated floating point
   values from a param file, allocating 'float' memory
   to store it, and assigning _dst to point to that
   new memory */
void load_model_paramlist(FILE* file, char* param_file, char* param_name, float** _dst, int num_params)
{
	int i;
	char buffer[64];
	float* dst = 0;
	*_dst = 0;
	strcpy(buffer, param_name);
	if (num_params > 0) {
		float* dst = (float*)malloc(sizeof(float)*num_params);
		strcat(buffer, " = %f");
		if (fscanf(file, buffer, &dst[0]) != 1)
			panic("Failed parsing %s in file %s", param_file, param_name);
		for (i = 1; i < num_params; ++i)
			if (fscanf(file, ",%f", &dst[i]) != 1)
				panic("Failed parsing %s in file %s", param_file, param_name);
		fscanf(file, "\n"); // eat up the carriage return
		*_dst = dst;
	} else {
		strcat(buffer, " = ");
		fscanf(file, buffer); // eat up the line
	}
}


/* Allocate a new, uninitialized model instance,
   with enough memory for all the required parameters, and
   then read the parameters from the given file */
deepbind_model_t* load_model(model_id_t id)
{
	char param_file[512];
	int major_version = 0;
	int minor_version = 0;
	deepbind_model_t* model = 0;
	FILE* file = 0;

	/* Find path to file*/
	strcpy(param_file, db_params_dir);
	id2str(id, param_file+strlen(param_file));
	strcat(param_file, ".txt");

	/* Open the file and parse each line */
	file = fopen(param_file, "r");
	if (!file)
		panic("Could not open param file %s.", param_file);

	if (fscanf(file, "# deepbind %d.%d\n", &major_version, &minor_version) != 2)
		panic("Expected parameter file to start with \"# deepbind 0.1\"");
	if (major_version != 0 || minor_version != 1)
		panic("Param file is for deepbind %d.%d, not 0.1", major_version, minor_version);

	model = (deepbind_model_t*)malloc(sizeof(deepbind_model_t));
	model->id = id;

	if (fscanf(file, "reverse_complement = %d\n", &model->reverse_complement) != 1) panic("Failed parsing reverse_complement in %s", param_file);
	if (fscanf(file, "num_detectors = %d\n",      &model->num_detectors) != 1)      panic("Failed parsing num_detectors in %s", param_file);
	if (fscanf(file, "detector_len = %d\n",       &model->detector_len) != 1)       panic("Failed parsing detector_len in %s", param_file);
	if (fscanf(file, "has_avg_pooling = %d\n",    &model->has_avg_pooling) != 1)    panic("Failed parsing has_avg_pooling in %s", param_file);
	if (fscanf(file, "num_hidden = %d\n",         &model->num_hidden) != 1)         panic("Failed parsing num_hidden in %s", param_file);

	load_model_paramlist(file, param_file, "detectors",  &model->detectors,  model->num_detectors*model->detector_len*4);
	load_model_paramlist(file, param_file, "thresholds", &model->thresholds, model->num_detectors);
	load_model_paramlist(file, param_file, "weights1",   &model->weights1,   get_num_hidden1(model)*get_num_hidden2(model));
	load_model_paramlist(file, param_file, "biases1",    &model->biases1,    get_num_hidden2(model));
	load_model_paramlist(file, param_file, "weights2",   &model->weights2,   model->num_hidden ? model->num_hidden : 0);
	load_model_paramlist(file, param_file, "biases2",    &model->biases2,    model->num_hidden ? 1 : 0);

	fclose(file);
	return model;
}


/* Free all memory associated with this model,
   including its parameters */
void free_model(deepbind_model_t* model)
{
	if (model->biases2)  free(model->biases2);
	if (model->weights2) free(model->weights2);
	free(model->biases1);
	free(model->weights1);
	free(model->thresholds);
	free(model->detectors);
	free(model);
}


/* Reverse complement a string in-place */
void reverse_complement(char* seq, int seq_len)
{
	int i, j;
	for (i = 0, j = seq_len-1; i <= j; ++i, --j) {
		unsigned char ci = (unsigned char)seq[i];
		unsigned char cj = (unsigned char)seq[j];
		seq[i] = base2comp[cj];
		seq[j] = base2comp[ci];
	}
}


/* Apply the model directly to the given DNA/RNA sequence string,
   e.g. "ACCGTCG". Return the score of the entire sequence. */
float apply_model(deepbind_model_t* model, char* seq, int seq_len)
{
	int n = seq_len;
	int m = model->detector_len;
	int d = model->num_detectors;
	int num_hidden1 = get_num_hidden1(model);
	int num_hidden2 = get_num_hidden2(model);
	int chunk0 = d*(n+m-1), chunk1 = num_hidden1, chunk2 = num_hidden2;
	float* mem = (float*)malloc(sizeof(float)*(chunk0 + chunk1 + chunk2));
	float* featuremaps = mem;
	float* hidden1     = mem+chunk0;
	float* hidden2     = mem+chunk0+chunk1;
	float* detectors   = model->detectors;
	float* thresholds  = model->thresholds;
	float* weights1    = model->weights1;
	float* biases1     = model->biases1;
	float* weights2    = model->weights2;
	float* biases2     = model->biases2;
	float  p;
	int i, j, k;

	/* Convolution, rectification */
	for (k = 0; k < d; ++k) {
		for (i = 0; i < n+m-1; ++i) {

			/* Convolve */
			float featuremap_ik = 0;
			for (j = 0; j < m; ++j) {
				char c = ((i-m+1)+j >= 0 && (i-m+1)+j < n) ? seq[(i-m+1)+j] : 'N';
				int index = base2index[(unsigned char)c];
				if (index == UNKNOWN_BASE) {
					for (index = 0; index < 4; ++index)
						featuremap_ik += .25f*detectors[indexof_detector_coeff(d, k, j, index)];
				} else {
					featuremap_ik += detectors[indexof_detector_coeff(d, k, j, index)];
				}
			}

			/* Shift and rectify */
			featuremap_ik += thresholds[k];
			if (featuremap_ik < 0)
				featuremap_ik = 0;

			featuremaps[indexof_featuremap_coeff(d, k, i)] = featuremap_ik;
		}
	}

	/* Pooling */
	if (model->has_avg_pooling) {
		for (k = 0; k < d; ++k) {
			float z_max = 0;
			float z_sum = 0;
			for (i = 0; i < n+m-1; ++i) {
				float featuremap_ik = featuremaps[indexof_featuremap_coeff(d, k, i)];
				z_sum += featuremap_ik;
				if (z_max < featuremap_ik)
					z_max = featuremap_ik;
			}
			hidden1[2*k+0] = z_max;
			hidden1[2*k+1] = z_sum /(float)(n+m-1);
		}
	} else {
		for (k = 0; k < d; ++k) {
			float z_max = 0;
			for (i = 0; i < n+m-1; ++i) {
				float featuremap_ik = featuremaps[indexof_featuremap_coeff(d, k, i)];
				if (z_max < featuremap_ik)
					z_max = featuremap_ik;
			}
			hidden1[k] = z_max;
		}
	}

	/* First hidden layer after convolution and pooling */
	for (j = 0; j < num_hidden2; ++j) {
		float h_j = biases1[j];
		for (i = 0; i < num_hidden1; ++i) {
			h_j += weights1[i*num_hidden2+j] * hidden1[i];
		}
		hidden2[j] = h_j;
	}

	if (num_hidden2 == 1) {
		/* No second hidden layer, so the lone hidden value is the final score */
		p = hidden2[0];
	} else {
		/* Second hidden layer, has its own biases, rectification, and weights */
		p = biases2[0];
		for (j = 0; j < num_hidden2; ++j) {
			float h_j = hidden2[j];
			if (h_j < 0)
				h_j = 0;
			p += weights2[j] * h_j;
		}
	}

	free(mem);
	return p;
}


/* Apply the model to the given DNA/RNA sequence string,
   e.g. "ACCGTCG". Return the score of the entire sequence. */
float scan_model(deepbind_model_t* model, char* seq, int seq_len, int window_size) {
	float scan_score = -10000.0f;
	int i;
	if (seq_len <= window_size)
		return apply_model(model, seq, seq_len);
	for (i = 0; i < seq_len-window_size+1; i++) {
		float score_i = apply_model(model, seq+i, window_size);
		if (score_i > scan_score) 
			scan_score = score_i;  
	}
	return scan_score;
}


int get_commandline_ids(model_id_t* ids, char** argv, char** argv_end)
{
	model_id_t id = { 0, 0 };
	int num_ids = 0;
	char buffer[1024];
	for (; argv < argv_end; ++argv) {
		char* arg = *argv;

		if (has_extension(arg, "ids")) {

			/* arg is a file full of ids */
			FILE* file = fopen(arg, "r");
			if (!file)
				panic("Could not open file %s.", arg);

			while (fgets(buffer, 1024, file)) {
				/* extract the next id, skipping comments */
				if (buffer[0] == '#')
					continue;
				buffer[1023] = '\0';
				trim_trailing_whitespace(buffer);
				id = str2id(buffer);
				if (ids)
					ids[num_ids] = id;
				num_ids++;
			}
			fclose(file);

		} else {

			if (argv+1 == argv_end) {
				/* the last arg is allowed to be an input file */
				if (is_seq_filename(arg))
					continue;
			}

			/* arg is directly an id */
			id = str2id(arg);
			if (ids)
				ids[num_ids] = id;
			num_ids++;
		}
	}
	return num_ids;
}


char* bases = "ACGT";

void randomSequence(int k, char* inputSequence, char* outputSequence) {
	int changePos = arc4random_uniform(k);
	char currNuc  = inputSequence[changePos];
	char newNuc   = bases[arc4random_uniform(4)];
	while (TRUE) {
		if (currNuc != newNuc) {
			break;
		}
		newNuc = bases[arc4random_uniform(4)];
	}
	strcpy(outputSequence, inputSequence);
	outputSequence[changePos] = newNuc;
}

void fullRandomSequence(int k, char* inputSequence, char* outputSequence) {
	strcpy(outputSequence, inputSequence);	
	for (int i=0; i<k; i++) {
		outputSequence[i] = bases[arc4random_uniform(4)];
	}
}

float scoreSeq(deepbind_model_t* model, char* seq, int seq_len, int window_size) {
	float score = scan_model(model, seq, seq_len, window_size);
	if (model->reverse_complement) {
		/* Reverse complement also needs to be scored. Take the max. */
		reverse_complement(seq, seq_len);
		float rscore = scan_model(model, seq, seq_len, window_size);
		if (rscore > score)
			score = rscore;
	}
	return score;
}

double mean(int* input, int* mod, int nElems) {
	double out = 0, count = 0;
	
	for (int i=0; i<nElems; ++i) {
		if (mod[i]) {
			out += input[i];
			count += 1;
		}
	}
	return (out/count);
}

int min(int* input, int* mod, int nElems) {
	int minval = INT_MAX;
	
	for (int i=0; i<nElems; ++i) {
		if (mod[i]) {
			if ( input[i] < minval) {
				minval = input[i];
			}
		}
	}
	return minval;
}

int main(int argc, char* argv[])
{
	int  i;
	int  head_flag = 1;
	int  echo_flag = 0;
	int  info_flag = 0;
	int  window_size = 0;
	int  average_flag = 0;
	int  num_models = 1;
	model_id_t* model_ids = 0;
	deepbind_model_t* currModel = 0;
	int   seq_len = 0;
	char* seq = 0;
	char* line = 0;
	FILE* in_file = 0;
	int   is_fasta = -1;
	int   done = 0;
	
	/* DOS vars */
	int varLen, nInit, minCount, nDivs, Io, In, totLen, iters;
	double f = 1, a, Eo, En, Emin, Emax, divSize, tol, fCrit;
	char* lFlank;
	char* rFlank;
	char* xo;
	char* xn;
	char* newTestSeq;
	int* H;
	int* RH;
	double* S;

	/* Need to read in args: modelID, varLen, lFlank, rFlank, Emin, Emax, divSize, nInit, tol, minCount, fCrit */
	if (argc != 12)
		panic("Incorrect number of arguments");  /* Expected at least one id */

	/* Initialize ACGT -> 0123 and ACGT -> TGCA conversion tables. */
	init_base2index_table();
	init_base2comp_table();

	/* Find out where the db/db.tsv and db/params/DXXXXX.YYY.txt files are stored */
	find_db_dirs();

	/* print out command line args */
/*	for (i = 1; i<12; ++i) {
		printf(argv[i]);
		printf("\n");
	}  */
	
	/* Parse command line inputs */
	model_ids = (model_id_t*)malloc(sizeof(model_id_t)*num_models);
	get_commandline_ids(model_ids, argv+1, argv+2);

	sscanf(argv[2], "%d", &varLen);
	lFlank = argv[3];
	rFlank = argv[4];
	sscanf(argv[5], "%lf", &Emin);
	sscanf(argv[6], "%lf", &Emax);
	sscanf(argv[7], "%lf", &divSize);
	sscanf(argv[8], "%d", &nInit);
	sscanf(argv[9], "%lf", &tol);
	sscanf(argv[10], "%d", &minCount);
	sscanf(argv[11], "%lf", &fCrit);
	
	/* Set up current model */
	currModel = load_model(model_ids[0]);

	/* Initialize variables */	
	nDivs = (int) (Emax-Emin)/divSize + 1;
	H  = (int*)malloc(nDivs*sizeof(int));
	RH = (int*)malloc(nDivs*sizeof(int));
	S  = (double*)malloc(nDivs*sizeof(double)); 
	for (i=0; i<nDivs; ++i) {
		H[i] = 0;
		RH[i]= 0;
		S[i] = 0;
	}
	window_size = (int)(currModel->detector_len*1.5);
	
	xo = (char*)malloc(varLen+1);
	xn = (char*)malloc(varLen+1);
	memset(xo, 'A', varLen);			/* initialize xo to poly-A */
	xo[varLen] = '\0';
	totLen = varLen + strlen(lFlank) + strlen(rFlank);
	newTestSeq = (char*)malloc(totLen+1);
	
	/* Start computing the reference histogram */
	strcpy(newTestSeq, lFlank);
	strcat(strcat(newTestSeq, xo), rFlank);
	Eo = (double) scoreSeq(currModel, newTestSeq, totLen, window_size);
	Io = (int) ceil((Eo+Emax)/divSize);
	
	/* Loop until appropriate div size is found */
	while (TRUE) {
		for (i=0; i<(nDivs*nInit); ++i) {
			/* propose new state */
			randomSequence(varLen, xo, xn);
			strcpy(newTestSeq, lFlank);
			strcat(strcat(newTestSeq, xn), rFlank);
			En = (double) scoreSeq(currModel, newTestSeq, totLen, window_size);
			In = (int) ceil((En+Emax)/divSize);
			/* calculate acceptance */
			a = exp(S[Io] - S[In]);
			if ( ((float) arc4random()/UINT_MAX) <= a) {
				/* accept current proposal */
				strcpy(xo, xn);
				Eo = En;
				Io = In;
			}
			/* perform updates */
			S[Io] += f;
			H[Io] += 1;
		}
		
		/* Check the number of divisions */
		int minIdx, maxIdx, isMinHit = 0;
		for (i=0; i<nDivs; ++i) {
			if (H[i]>((int).1*nInit)) {
				if (isMinHit==0) {
					minIdx = i;
					isMinHit = 1;
				} else {
					maxIdx = i;
				}
			}
		}
				
		/* Minimum number of divisions met. Exit */
		if (maxIdx-minIdx+1 >= 100) {
			break;
		}
		
		fprintf(stderr, "minIdx: %d\tmaxIdx: %d\n", minIdx, maxIdx);
		fprintf(stderr, "Emin: %lf\tEmax: %lf\n", Emin+minIdx*divSize, Emin+maxIdx*divSize);
				
		/* Else, reallocate */
		Emax = Emin+(maxIdx+10)*divSize;
		Emin = Emin+(minIdx-10)*divSize;
		divSize = divSize/2;
		if (fabs(Emin) > fabs(Emax)) {
			Emin = -fabs(Emin);
			Emax = fabs(Emin);
		} else {
			Emin = -fabs(Emax);
			Emax = fabs(Emax);
		}
		
		nDivs = (int) (Emax-Emin)/divSize + 1;
		free(H);
		free(RH);
		free(S);
		
		fprintf(stderr, "Emin: %lf\tEmax: %lf\tdivSize: %lf\tnDivs: %d\n", Emin, Emax, divSize, nDivs);
				
		H  = (int*)malloc(nDivs*sizeof(int));
		RH = (int*)malloc(nDivs*sizeof(int));
		S  = (double*)malloc(nDivs*sizeof(double)); 
		for (i=0; i<nDivs; ++i) {
			H[i] = 0;
			RH[i]= 0;
			S[i] = 0;
		}
	}
	
	fprintf(stderr, "Reference Histogram Computed for n = %d samples, Emin = %lf, Emax = %lf, divSize = %lf\n", nDivs*nInit, Emin, Emax, divSize);
	/* Set reference histogram and reset H, S */
	for (i=0; i<nDivs; ++i) {
		RH[i] = (H[i]>((int).1*nInit)) ? 1 : 0;
		H[i] = 0;
		S[i] = 0;
	}

	/* Begin DOS computation */
	randomSequence(varLen, xo, xo);
	strcpy(newTestSeq, lFlank);
	strcat(strcat(newTestSeq, xo), rFlank);
	Eo = (double) scoreSeq(currModel, newTestSeq, totLen, window_size);
	Io = (int) ceil((Eo+Emax)/divSize);
	iters = 0;
	
	while (f > tol) {
		randomSequence(varLen, xo, xn);
		strcpy(newTestSeq, lFlank);
		strcat(strcat(newTestSeq, xn), rFlank);
		En = (double) scoreSeq(currModel, newTestSeq, totLen, window_size);
		In = (int) ceil((En+Emax)/divSize);
		/* calculate acceptance */
		a = exp(S[Io] - S[In]);
		if ( ((float) arc4random()/UINT_MAX) <= a) {
			/* accept current proposal */
			strcpy(xo, xn);
			Eo = En;
			Io = In;
		}
		/* perform updates */
		S[Io] += f;
		H[Io] += 1;
		iters += 1;
		if (min(H, RH, nDivs) > minCount) {
			if (min(H, RH, nDivs) >= fCrit*mean(H, RH, nDivs)) {
				fprintf(stderr, "Flatness criteria reached for f = %lf in %d iterations.\n", f, iters);
				f = f/2;
				iters = 0;
				for (i=0; i<nDivs; ++i) {
					H[i] = 0;
				}
			}
		}
		
		if (iters % 5000 ==0) {
			fprintf(stderr, "Iters: %d\t%d\t%lf\t%lf\n", iters, min(H, RH, nDivs), mean(H, RH, nDivs), min(H, RH, nDivs)/mean(H, RH, nDivs));
		}
	}
	fprintf(stderr, "Density of States computed. Entropy S is given below:\n");
	for (i=0; i<nDivs; ++i) {
		fprintf(stdout, "%lf\t%lf\n", Emin+i*divSize, S[i]);
	}

	/* Free memory */
	free(seq);
	free(line);
	free(model_ids);
	free(currModel);
	free(xo);
	free(xn);
	free(newTestSeq);
	free(H);
	free(RH);
	free(S);
	return 0;
}