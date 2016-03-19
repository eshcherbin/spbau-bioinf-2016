#include <cstdio>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include "htslib/sam.h"

using namespace std;

const int DEFAULT_GRAPH_STEP_SIZE = 1000;

int len, *cnt, step_size;
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

    len = hdr->target_len[0];
    cnt = new int[len];
    fill(cnt, cnt + len, 0);

    record = bam_init1();
}

void proc_record()
{
    for (int i = record->core.pos; i < bam_endpos(record); i++)
        cnt[i]++;
}

void output_stats()
{
    int max_cvg = 0;
    double cvg = 0, avg = 0;

    for (int i = 0; i < len; i++)
    {
        if (cnt[i] > max_cvg)
            max_cvg = cnt[i];
        avg += cnt[i];
        if (cnt[i])
            cvg += 1;
    }
    avg /= len;
    cvg /= len;

    printf("maximal coverage: %d reads\n", max_cvg);
    printf("average coverage: %.2f reads/position\n", avg);
    printf("%.2f%% positions covered\n", cvg * 100.);

    ofstream data("plot.dat");
    data << max_cvg << ' ' << fixed << setprecision(2) << avg << ' ' << cvg << '\n';
    for (int l = 0; l < len; l += step_size)
    {
        int r = min(l + step_size, len);
        double avg_step = 0;
        for (int i = l; i < r; i++)
            avg_step += cnt[i];
        avg_step /= (r - l);
        data << l << ' ' << fixed << setprecision(5) << avg_step << '\n';
    }
}

void destroy()
{
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
