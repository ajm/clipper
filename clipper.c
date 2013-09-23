#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>


#define CLIPPER_DEFAULT_SUFFIX ".filtered"

typedef enum {
    CLIPPER_FASTA,
    CLIPPER_FASTQ
} filetype_t;

typedef struct {
    int min_length;         // filter out short reads post-soft clip
    int qual_trim;          // parameter to soft clip algorithm
    int qual_filt;        // number of bases with QUAL < 3

    filetype_t type;
    char* output_dir;
    char* suffix;

    int paired_flg;         // maintain pairs, maybe print orphans?
    int remove_ambig_flg;   // filter reads containing 'N'
    int remove_adapter_flg; 
    int paranoid_flg;

    int phred_offset;

    char* adapter;

    int verbose_flg;

} options_t;

typedef struct {
    char* buf;          // buffer 
    size_t buf_cap;     // length of buffer

    char* str;          // start of string, on read - reset to buf[0]
    ssize_t str_len;    // length of string

} strbuf_t;

typedef struct {
    strbuf_t id;
    strbuf_t seq;
    strbuf_t qual;

} fq_entry_t;

typedef struct {
    FILE* in;
    FILE* out;

    fq_entry_t entry;
    int error;

    char* input_name;
    char* output_name;

} fq_file_t;

/*
 sick of typing fprintf(stderr, ... etc...
 */
void err(const char* fmt, ...) {
    va_list v;
    va_start(v, fmt);

    vfprintf(stderr, fmt, v);

    va_end(v);
}

void init_options(options_t* opt) {
    opt->min_length = 0;
    opt->qual_trim = -1;
    opt->qual_filt = -1;
    opt->type = CLIPPER_FASTQ;
    opt->output_dir = NULL;
    opt->suffix = NULL;
    opt->paired_flg = 0;
    opt->remove_ambig_flg = 0;
    opt->remove_adapter_flg = 0;
    opt->paranoid_flg = 1; //0;
    opt->adapter = NULL;
    opt->verbose_flg = 0;
    opt->phred_offset = 33;
}

void destroy_options(options_t* opt) {
    if(opt->output_dir)
        free(opt->output_dir);

    if(opt->suffix)
        free(opt->suffix);

    if(opt->adapter)
        free(opt->adapter);
}

void print_options(options_t* opt) {
    fprintf(stdout,
            "min_length = %d\n"
            "qual_trim = %d\n"
            "qual_filt = %d\n"
            ,
            opt->min_length,
            opt->qual_trim,
            opt->qual_filt
            );
}

void init_strbuf(strbuf_t* sb) {
    sb->buf = NULL;
    sb->str = NULL;
    sb->buf_cap = 0;
    sb->str_len = 0;
}

void destroy_strbuf(strbuf_t* sb) {
    if(sb->buf)
        free(sb->buf);
}

void init_entry(fq_entry_t* fq) {
    init_strbuf(&fq->id);
    init_strbuf(&fq->seq);
    init_strbuf(&fq->qual);
}

void destroy_entry(fq_entry_t* fq) {
    destroy_strbuf(&fq->id);
    destroy_strbuf(&fq->seq);
    destroy_strbuf(&fq->qual);
}

int init_fq(fq_file_t* fq, char* dir, char* filename, char* suffix) {
    fq->error = 0;

    if((fq->input_name = strdup(filename)) == NULL) {
        err("Error: during strdup: %s\n", strerror(errno));
        return -1;
    }

    if(!dir) {
        fq->output_name = (char*) malloc(strlen(filename) + strlen(suffix) + 1);
        if(fq->output_name == NULL) {
            err("Error: during strdup: %s\n", strerror(errno));
            return -1;
        }

        strcpy(fq->output_name, filename);
        strcat(fq->output_name, suffix);
    }
    else {
        char* tmp = basename(filename);
        
        fq->output_name = (char*) malloc(strlen(dir) + strlen(tmp) + strlen(suffix) + 1);
        if(fq->output_name == NULL) {
            err("Error: during strdup: %s\n", strerror(errno));
            return -1;
        }

        strcpy(fq->output_name, dir);
        strcat(fq->output_name, tmp);
        strcat(fq->output_name, suffix);
    }

    fq->in = fopen(fq->input_name, "r");
    if(fq->in == NULL) {
        err("Error: could not open '%s' for reading: %s\n", 
                fq->input_name, strerror(errno));
        return -1;
    }

    fq->out = fopen(fq->output_name, "w");
    if(fq->out == NULL) {
        err("Error: could not open '%s' for writing: %s\n", 
                fq->output_name, strerror(errno));
        return -1;
    }

    init_entry(&fq->entry);

    return 0 ;
}

void destroy_fq(fq_file_t* fq) {
    if(fq->in && fclose(fq->in) == EOF) {
        err("Error: closing file %s : %s\n", fq->in, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(fq->out && fclose(fq->out) == EOF) {
        err("Error: closing file %s : %s\n", fq->out, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(fq->input_name)
        free(fq->input_name);

    if(fq->output_name)
        free(fq->output_name);

    destroy_entry(&fq->entry);
}

void rtrim_whitespace(char* s, ssize_t* len) {
    int i, end;

    for(i = *len-1, end = 0; i > -1; --i) {
        if(! isspace(s[i]))
            break;
        
        end++;
    }

    s[*len-end] = '\0';
    *len -= end;
}

int get_entry(FILE* f, fq_entry_t* fq) {
    int i;
    strbuf_t* sb;

    for(i = 0; i < 4; ++i) {

        switch(i) {
            case 0: sb = &fq->id; break;
            case 1: sb = &fq->seq; break;
            case 2: 
            case 3: sb = &fq->qual; break;
            default:
                return -1;
        }

        if(!((sb->str_len = getline(&sb->buf, &sb->buf_cap, f)) > 0))
            return -1;

        if((i == 0) && (sb->buf[0] != '@'))
            return -1;

        if((i == 2) && (sb->buf[0] != '+'))
            return -1;

        // getline always includes the delimiter (in this case \n)
        // so null terminate the string one character early
        //   - changed to all whitespace
        //sb->buf[sb->str_len-1] = '\0';
        //sb->str_len -= 1;

        rtrim_whitespace(sb->buf, &sb->str_len);

        // reset start of string
        sb->str = sb->buf;
    }

/*
    // apply offset to quality values
    for(i = 0; i < fq->qual.str_len; ++i) {
        fq->qual.buf[i] -= 33;
    }
*/
    
    if(fq->seq.str_len != fq->qual.str_len) {
        err("Error: sequence and qualities were different lengths! (%s)\n", fq->id.str);
        exit(EXIT_FAILURE);
    }

    return 0;
}

int readnext_fq(fq_file_t* fq) {
    return get_entry(fq->in, &fq->entry);
}

int done_fq(fq_file_t* fq) {
    return fgetc(fq->in) == EOF;
}

int paranoid_dna(char* s, ssize_t len) {
    int i;

    for(i = 0; i < len; ++i) {        
        switch(s[i]) {
            case 'A':
            case 'G':
            case 'C':
            case 'T':
            case 'N':
                break;
            default:
                err("Warning: found '%c' in sequence!\n", s[i]);
                return -1;
        }
    }

    return 0;
}

int paranoid_quality(char* s, ssize_t len, int offset) {
    int i, hi, lo;
    
    // https://en.wikibooks.org/wiki/Next_Generation_Sequencing_%28NGS%29/Pre-processing
    if(offset == 33) {
        lo = 33;
        hi = 74;
    }
    else {
        lo = 59;
        hi = 104;
    }

    for(i = 0; i < len; ++i) {
        if(s[i] < lo || s[i] > hi) {
            err("Warning: found '%c' in quality scores! (phred%d = ('%c' to '%c'))\n", 
                s[i], offset, lo, hi);
            return -1;
        }
    }

    return 0;
}

int process_current_fq(fq_file_t* fq, options_t* opt) {
    fq_entry_t* e = &fq->entry;
    int i,j,qual,max_val,argmax;
    char* c;

    i = -1;
    j = e->seq.str_len;
    c = e->seq.str;

    // trim leading + trailing Ns
    while(c[++i] == 'N');
    while(c[--j] == 'N');

    //printf("%d - %d : %s\n", i, j, e->seq.str); 

    e->seq.str[j+1] = '\0';
    e->seq.str += i;
    e->seq.str_len = j - i;

    e->qual.str[j+1] = '\0';
    e->qual.str += i;
    e->qual.str_len = j - i;

    //printf("%d - %d : %s\n", i, j, e->seq.str);
    //printf("\n");

    // paranoid?
    if(opt->paranoid_flg) {
        if(paranoid_dna(e->seq.str, e->seq.str_len) < 0)
            return -1;

        if(paranoid_quality(e->qual.str, e->qual.str_len, opt->phred_offset) < 0) 
            return -1;
    }

    // remove adapter - adapter AA : XXXAAxxx --> XXXxxx
    // XXX not implemented

    // primer check - return -1 if present
    // XXX not implemented

    // contains 'N' - return -1 if bad
    if(opt->remove_ambig_flg) {
        c = e->seq.str;
        
        for(i = 0; i < e->seq.str_len; ++i) {
            if(c[i] == 'N')
                return -1;
        }
    }

    // quality trim - only changes the location of '\0'
    // argmax_x{
    //    \sum_{i=x+1}^l
    //       (INT-q_i)
    // } if q_l<INT where l is the original read length.
    if(opt->qual_trim > -1) {
        c = e->qual.str;
        j = max_val = 0;
        argmax = e->qual.str_len;

        if((c[argmax] - opt->phred_offset) < opt->qual_trim) {
            for(i = e->qual.str_len; i > -1; --i) {
                qual = c[i] - opt->phred_offset;
                j += (opt->qual_trim - qual);

            // done in BWA, but not SGA
//          if(j < 0)
//              break;

                if(j > max_val) {
                    max_val = j;
                    argmax = i;
                }
            }

            e->seq.str[argmax] = '\0';
            e->seq.str_len = argmax;
            e->qual.str[argmax] = '\0';
            e->qual.str_len = argmax;
        }
    }

    // quality filter - return -1 if bad
    if(opt->qual_filt > -1) { 
        j = 0;
        c = e->qual.str;

        for(i = 0; i < e->seq.str_len; ++i) {
            qual = c[i] - opt->phred_offset;
            if(qual <= 3)
                j++;
        }

        if(j > opt->qual_filt)
            return -1;
    }

    // length filter - return -1 if bad
    if(e->seq.str_len < opt->min_length)
        return -1;

    return 0;
}

int output_fq(fq_file_t* fq) {
    fq_entry_t* e = &fq->entry;

    return fprintf(fq->out, 
                "%s\n%s\n+\n%s\n", 
                e->id.str, 
                e->seq.str, 
                e->qual.str);
}

void usage(const char* prog) {
    fprintf(stderr, 
            "Usage: %s [options] file1, [file2 ...]\n"
            "  -q INT\t--softclip=INT\t\tquality trim threshold\n"
            "  -f INT\t--qualityfilter=INT\tquality filter threshold\n"
            "  -m INT\t--minlength=INT\t\tminimum length threshold\n"
            "  -p\t\t--paired\t\tpaired end reads\n"
            "  -n\t\t--filterambiguous\tfilter reads containing 'N's\n"
            "  -a STR\t--removeadapter=STR\tremove fwd adapter\n"
            "  -t STR\t--type=STR\t\tfile type (fasta,fa,fastq,fq)\n"
            "  -z\t\t--phred64\t\tphred offset of 64 (default 33)\n"
            //"  -x\t\t--paranoid\t\tcheck all characters in seq + qual are within expected ranges\n"
            "  -d DIR\t--outputdir\t\tspecify output directory (default '.')\n"
            "  -s STR\t--suffix=STR\tspecify suffix appended to output file (default = " CLIPPER_DEFAULT_SUFFIX ")\n"
            "  -v\t\t--verbose\t\tverbose\n"
            "  -h\t\t--help\t\t\tshow this help message\n"
            "\n"
            , prog);
    exit(EXIT_FAILURE);
}

int parse_int(const char* s, int* dst) {
    char* endptr;
    
    *dst = strtol(s, &endptr, 10);

    return *endptr == '\0';
}

void set_output_dir(options_t* opt, const char*s) {
    struct stat tmp;
    int len;

    if(stat(s, &tmp) < 0) {
        err("Error: could not stat '%s': %s\n", s, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if(! S_ISDIR(tmp.st_mode)) {
        err("Error: '%s' is not a valid directory!\n", s);
        exit(EXIT_FAILURE);
    }

    // ensure there is a trailing '/'
    len = strlen(s);
    opt->output_dir = (char*) malloc(s[len] == '/' ? len : len+1);
    if(!opt->output_dir) {
        err("Error: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    strcpy(opt->output_dir, s);
    if(s[len] != '/') {
        opt->output_dir[len] = '/';
        opt->output_dir[len+1] = '\0';
    }

#ifdef DEBUG
    printf("output dir = %s\n", opt->output_dir);
#endif
}

void handle_cli(int* argc, char*** argv, options_t* opt) {
    int ch;
    int exit_flg = 0;

    static struct option longopts[] = {
        { "softclip",        required_argument, NULL, 'q' },
        { "qualityfilter",   required_argument, NULL, 'f' },
        { "minlength",       required_argument, NULL, 'f' },
        { "output",          required_argument, NULL, 'o' },
        { "filterambiguous", no_argument,       NULL, 'n' },
        { "removeadapter",   required_argument, NULL, 'a' },
        { "type",            required_argument, NULL, 't' },
        { "paired",          no_argument,       NULL, 'p' },
        { "verbose",         no_argument,       NULL, 'v' },
        { "help",            no_argument,       NULL, 'h' },
        { "phred64",         no_argument,       NULL, 'z' },
        { "paranoid",        no_argument,       NULL, 'x' },
        { "outputdir",       required_argument, NULL, 'd' },
        { "suffix",          required_argument, NULL, 's' },
        { NULL, 0, NULL, 0 }
    };

    init_options(opt);

    while((ch = getopt_long(*argc, *argv, ":q:f:m:t:a:o:pnhvzxd:s:", longopts, NULL)) != -1) {
        switch(ch) {
            case 'q':
                if(! parse_int(optarg, &opt->qual_trim)) {
                    err("Error: option '%c' expects INT (read '%s')\n", ch, optarg);
                    exit_flg = 1;
                }
                break;

            case 'f':
                if(! parse_int(optarg, &opt->qual_filt)) {
                    err("Error: option '%c' expects INT (read '%s')\n", ch, optarg);
                    exit_flg = 1;
                }
                break;

            case 'm':
                if(! parse_int(optarg, &opt->min_length)) {
                    err("Error: option '-%c' expects INT (read '%s')\n", ch, optarg);
                    exit_flg = 1;
                }
                break;

            case 't':
                err("Error: not implemented! must be fastq!\n");
                exit(EXIT_FAILURE);

                break;
            case 'a':
                err("Error: not implemented!\n");
                exit(EXIT_FAILURE);

                break;
            case 'o':
                err("Error: not implemented!\n");
                exit(EXIT_FAILURE);

                break;
            case 'd':
                set_output_dir(opt, optarg);
                break;
            case 's':
                if(strchr(optarg, '/') != NULL) {
                    err("Error: suffix cannot contain '/'!\n");
                    exit(EXIT_FAILURE);
                }
                opt->suffix = strdup(optarg);
                if(!opt->suffix) {
                    err("Error: %s\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
                break;
            case 'z':
                opt->phred_offset = 64;
                break;
            case 'x':
                opt->paranoid_flg = 1;
                break;
            case 'p':
                opt->paired_flg = 1;
                break;
            case 'n':
                opt->remove_ambig_flg = 1;
                break;
            case 'v':
                opt->verbose_flg = 1;
                break;
            case 'h':
                usage(*argv[0]);
                break;
            case ':':
                err("Error: missing parameter for option '-%c'\n", optopt);
                exit_flg = 1;
                break;
            case '?':
            default:
                usage(*argv[0]);
        }
    }

    *argc -= optind;
    *argv += optind;

    if(opt->paired_flg && *argc != 2) {
        err("Error: paired mode assumes two file names, you provided %d\n", *argc);
        exit(EXIT_FAILURE);
    }
/*
    if(!((opt->qual_trim > -1) && (opt->qual_filt > -1))) {
        err("Error: both options -q and -f must be set!\n");
        exit(EXIT_FAILURE);
    }
*/
    if(*argc == 0) {
        err("Error: you must specify at least one file!\n");
        exit(EXIT_FAILURE);
    }

    if(exit_flg) 
        exit(EXIT_FAILURE);
}

void open_all(int num_files, char** file_names, char* dir, char* suffix, fq_file_t* files[]) {
    int i;
    fq_file_t* f;
    
    *files = (fq_file_t*) malloc(sizeof(fq_file_t) * num_files);
    if(! *files) {
        err("Error: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < num_files; ++i) {
        f = *files + i;

        if(init_fq(f, dir, file_names[i], suffix) < 0) {
            //err("Error: could not init fastq file object for %s: %s\n", 
            //        file_names[i], strerror(errno));
            exit(EXIT_FAILURE);
        }

#ifdef DEBUG
        fprintf(stderr, 
                "%s -> %s\n", 
                f->input_name, 
                f->output_name);
#endif
    }
}

void close_all(fq_file_t* files[], int num_files) {
    int i;

    for(i = 0; i < num_files; ++i) {
        destroy_fq(*files + i);
    }

    free(*files);
}

void process_paired(fq_file_t* files[], int num_files, options_t* opt) {
    fq_file_t* left = *files;
    fq_file_t* right = *files + 1;

    while(readnext_fq(left) == 0  && readnext_fq(right) == 0) {
        if(process_current_fq(left, opt) == 0 && process_current_fq(right, opt) == 0) {
            if(output_fq(left) < 0 || output_fq(right) < 0) {
                err("Error: write to %s failed: %s\n", 
                        left->error ? left->output_name : right->output_name, 
                        strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    }

    // check if both files are finished, the while condition is false
    // if only one does
    if(! done_fq(left)) {
        err("Error: %s has more lines than %s!\n", left->input_name, right->input_name);
        exit(EXIT_FAILURE);
    }

    if(! done_fq(right)) {
        err("Error: %s has more lines than %s!\n", right->input_name, left->input_name);
        exit(EXIT_FAILURE);
    }
}

void process_nonpaired(fq_file_t* files[], int num_files, options_t* opt) {
    int i;
    fq_file_t* current;

    for(i = 0; i < num_files; ++i) {
        current = *files + i;

        while(readnext_fq(current) == 0) {
            if(process_current_fq(current, opt) == 0) {
                if(output_fq(current) < 0) {
                    err("Error: write to %s failed: %s\n", 
                            current->output_name, 
                            strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    options_t opt;
    fq_file_t* files;

    handle_cli(&argc, &argv, &opt);
    open_all(argc, 
             argv, 
             opt.output_dir, 
             opt.suffix ? opt.suffix : CLIPPER_DEFAULT_SUFFIX, 
             &files);

    if(opt.paired_flg)
        process_paired(&files, argc, &opt);
    else
        process_nonpaired(&files, argc, &opt);

    close_all(&files, argc);
    destroy_options(&opt);

    return EXIT_SUCCESS;
}

