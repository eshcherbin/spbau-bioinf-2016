#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "htslib/sam.h"

void print_help()
{
    printf("usage: coverage_graph alignment\n");
    exit(0);
}

int main(int argc, char **argv)
{
    if (argc == 1 || !strcmp(argv[1], "-h"))
        print_help();

    char *fn = argv[1];
    samFile *fp = sam_open(fn, "r");

    bam_hdr_t *hdr = sam_hdr_read(fp);
    printf("%d\n", hdr->target_len[0]);

    bam_hdr_destroy(hdr);
    sam_close(fp);

    return 0;
}
