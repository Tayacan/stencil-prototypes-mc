#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "common.c"
#include "smooth_3d.c"

void smooth_3d_experiment(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    const uint64_t TW,
    const uint64_t TH,
    const uint64_t TD
) {
    float *xs_in, *xs_out, *xs_res;

    const size_t mem_size = sizeof(float) * W * H * D;
    xs_in  = (float*)malloc(mem_size);
    xs_out = (float*)malloc(mem_size);
    xs_res = (float*)malloc(mem_size);
    assert(xs_in != NULL);
    assert(xs_out != NULL);
    assert(xs_res != NULL);
    write_random_1d(W * H * D, xs_in);

    smooth_3d_untiled(W, H, D, xs_in, xs_res);
    measure_smooth_3d_tiled(W, H, D, TW, TH, TD, xs_in, xs_out, xs_res);
    measure_smooth_3d_tiled_partitioned(W, H, D, TW, TH, TD, xs_in, xs_out, xs_res);

    free(xs_in);
    free(xs_out);
    free(xs_res);
}

int main(
    int          argc,
    const char** argv
) {
    if(argc < 4) {
        smooth_3d_experiment(1000, 1000, 100, 16, 64, 64);
    } else {
        int tw = atoi(argv[1]);
        int th = atoi(argv[2]);
        int td = atoi(argv[3]);
        smooth_3d_experiment(1000, 1000, 100, tw, th, td);
    }
}
