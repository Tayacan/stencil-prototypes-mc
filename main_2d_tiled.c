#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NUM_RUNS 10

#include "common.c"
#include "smooth_2d.c"

void smooth_2d_experiment(
    const uint64_t W,
    const uint64_t H,
    const uint64_t TW,
    const uint64_t TH
) {
    float *xs_in, *xs_out, *xs_res; // , **d_in, **d_out, **d_tmp;

    const size_t mem_size = sizeof(float) * W * H;
    xs_in  = (float*)malloc(W * H * sizeof(float *));
    xs_out = (float*)malloc(W * H * sizeof(float *));
    xs_res = (float*)malloc(W * H * sizeof(float *));
    assert(xs_in != NULL);
    assert(xs_out != NULL);
    assert(xs_res != NULL);


    write_random_1d(W * H, xs_in);

    smooth_2d_untiled(W, H, xs_in, xs_res);
    measure_smooth_2d_tiled(W, H, TW, TH, xs_in, xs_out, xs_res);
    measure_smooth_2d_tiled_partitioned(W, H, TW, TH, xs_in, xs_out, xs_res);

    free(xs_in);
    free(xs_out);
    free(xs_res);
}

int main(
    int          argc,
    const char** argv
) {
    if(argc < 3) {
      //smooth_2d_experiment(1000000, 100, 256, 256);
      smooth_2d_experiment(10000, 10000, 256, 256);
    } else {
      int tw = atoi(argv[1]);
      int th = atoi(argv[2]);
      //smooth_2d_experiment(1000000, 100, tw, th);
      smooth_2d_experiment(10000, 10000, tw, th);
    }
}

