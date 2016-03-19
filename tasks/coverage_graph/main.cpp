#include <cstdio>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include "htslib/sam.h"

using namespace std;

const int DEFAULT_GRAPH_STEP_SIZE = 1000;

int **cnt, step_size;
char *fn;

samFile *fp;
bam_hdr_t *hdr;
bam1_t *record;

void print_help()
{
    printf("usage: coverage_graph alignment [step]\n");
    exit(0);
}

void process_args(int argc, char **argv)
{
    if (argc == 1 || !strcmp(argv[1], "-h"))
        print_help();

    if (argc == 3)
        step_size = atoi(argv[2]);
    else
        step_size = DEFAULT_GRAPH_STEP_SIZE;

    fn = argv[1];
}

void init()
{
    fp = sam_open(fn, "r");
    hdr = sam_hdr_read(fp);

    cnt = new int*[hdr->n_targets];
    for (int tg = 0; tg < hdr->n_targets; tg++)
    {
        cnt[tg] = new int[hdr->target_len[tg]];
        fill(cnt[tg], cnt[tg] + hdr->target_len[tg], 0);
    }

    record = bam_init1();
}

void proc_record()
{
    if (record->core.flag & BAM_FUNMAP) 
        return;
    int pos = record->core.pos;
    for (int i = 0; i < record->core.n_cigar; i++)
    {
        uint32_t c = bam_get_cigar(record)[i];
        if (bam_cigar_type(c) & 1)
        {
            if (bam_cigar_type(c) & 2)
            {
                for (int j = 0; j < (int) bam_cigar_oplen(c); j++)
                    cnt[record->core.tid][pos++]++;
            }
            else
                pos += bam_cigar_oplen(c);
        }
    }
}

void output_stats()
{
    for (int tg = 0; tg < hdr->n_targets; tg++)
    {
        string fn = "plot_" + string(hdr->target_name[tg]) + ".dat";
        ofstream data(fn.c_str());
        printf("reference sequence %s:\n", hdr->target_name[tg]);
        data << hdr->target_name[tg] << '\n';

        int max_cvg = 0;
        double cvg = 0, avg = 0;

        for (int i = 0; i < (int) hdr->target_len[tg]; i++)
        {
            if (cnt[tg][i] > max_cvg)
                max_cvg = cnt[tg][i];
            avg += cnt[tg][i];
            if (cnt[tg][i])
                cvg += 1;
        }
        avg /= hdr->target_len[tg];
        cvg /= hdr->target_len[tg];

        printf("maximal coverage: %d reads\n", max_cvg);
        printf("average coverage: %.2f reads/position\n", avg);
        printf("%.2f%% positions covered\n", cvg * 100.);

        data << max_cvg << ' ' << fixed << setprecision(2) << avg << ' ' << cvg << '\n';
        for (int l = 0; l < (int) hdr->target_len[tg]; l += step_size)
        {
            int r = min(l + step_size, (int) hdr->target_len[tg]);
            double avg_step = 0;
            for (int i = l; i < r; i++)
                avg_step += cnt[tg][i];
            avg_step /= (r - l);
            data << l << ' ' << fixed << setprecision(5) << avg_step << '\n';
        }
    }
}

void destroy()
{
    for (int tg = 0; tg < hdr->n_targets; tg++)
        delete [] cnt[tg];
    delete [] cnt;
    bam_destroy1(record);
    bam_hdr_destroy(hdr);
    sam_close(fp);
}

int main(int argc, char **argv)
{
    process_args(argc, argv);

    init();

    while (!sam_read1(fp, hdr, record))
        proc_record();

    output_stats();

    destroy();

    return 0;
}
