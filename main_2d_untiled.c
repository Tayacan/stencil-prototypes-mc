#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "common.c"
#include "smooth_2d.c"

void smooth_2d_experiment(
    const uint64_t W,
    const uint64_t H
) {
    float *xs_in, *xs_out, *xs_res; // , **d_in, **d_out, **d_tmp;

    const size_t mem_size = sizeof(float) * W * H;
    xs_in  = (float*)malloc(mem_size);
    xs_out = (float*)malloc(mem_size);
    xs_res = (float*)malloc(mem_size);
    assert(xs_in != NULL);
    assert(xs_out != NULL);
    assert(xs_res != NULL);


    write_random_1d(W * H, xs_in);

    measure_smooth_2d_untiled(W, H, xs_in, xs_res);
    measure_smooth_2d_untiled_partitioned(W, H, xs_in, xs_out, xs_res);

    free(xs_in);
    free(xs_out);
    free(xs_res);
}

int main(
    int          argc,
    const char** argv
) {
    smooth_2d_experiment(10000, 10000);
}

