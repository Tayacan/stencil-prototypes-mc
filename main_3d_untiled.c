#include <assert.h>
#include <stdint.h>

#include "common.c"
#include "smooth_3d.c"

void smooth_3d_experiment(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D
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

    measure_smooth_3d_untiled(W, H, D, xs_in, xs_res);
    measure_smooth_3d_untiled_partitioned(W, H, D, xs_in, xs_out, xs_res);

    free(xs_in);
    free(xs_out);
    free(xs_res);
}

int main(
    int          argc,
    const char** argv
) {
    smooth_3d_experiment(1000, 1000, 100);
}
