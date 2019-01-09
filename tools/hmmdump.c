#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define STB_DEFINE
#include "stb.h"

void print_usage(const char *program_name)
{
	printf("Usage:\n");
	printf("%s -f {Pfam-A.hmm} \n", program_name);
	printf("Options: -h display usage info\n");
	printf("         -f file containing the HMM library\n");
}

int main(int argc, char const *argv[])
{
	char **opts = NULL;
    opts = stb_getopt_param(&argc, argv, "fh");
    int length = 0;
    char *infilename = NULL;
    char *search = NULL;
    char **inhmmfile;
    FILE *fout;
    char *filename;

    /* get all arguments and options, if specified (in param)*/
    int i = 0;
    for (i = 0; i < stb_arr_len(opts) - 1; i++) {
        switch(opts[i][0]) {
            case 'f':
                infilename = opts[i]+1;
                break;
            case 'h':
            	print_usage(argv[0]);
                return 0;
            default:
                printf("Unsupported option %c\n", opts[i][0]);
                return 0;
        }
    }

    if (!infilename) {
        print_usage(argv[0]);
        return 0;
    }

    inhmmfile = stb_stringfile(infilename, &length);

    for (i = 0; i < length; i++) {
        /* Print record */
        if (strncmp(inhmmfile[i], "HMMER3/", 7) != 0) {
            fprintf(stderr, "File does not seem to be a valid HMM file\n");
            return 1;
        }
        stb_asprintf(&filename, "%s.hmm", &inhmmfile[i+1][6]);
        fout = fopen(filename, "w");
        while (inhmmfile[i][0] != '/') { 
            fprintf(fout, "%s\n", inhmmfile[i]);
            i++;
        }
        fprintf(fout, "//\n", inhmmfile[i]);
        fclose(fout);
        free(filename);
    }
    

    stb_getopt_free(opts);

	return 0;
}