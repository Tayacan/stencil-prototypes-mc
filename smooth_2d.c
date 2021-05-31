#ifndef SMOOTH_2D_SEQ
#define SMOOTH_2D_SEQ

#include "common.c"

#define RADX 4
#define RADY 4

void smooth_2d_untiled_partitioned(
    const int W,
    const int H,
    float* xs_in,
    float* xs_out
) {
    const int l = (2 * RADX + 1) * (2 * RADY + 1);

    // Top edge
#pragma omp parallel for collapse(2)
    for(int i = 0; i < RADY; i++) {
        for(int j = 0; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

    // Bottom edge
#pragma omp parallel for collapse(2)
    for(int i = H-RADY; i < H; i++) {
        for(int j = 0; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }


#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i++) {
        // Left edge
        for(int j = 0; j < RADX; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i++) {
        // Right edge
        for(int j = W-RADX; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i++) {
        for(int j = RADX; j < W-RADX; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[(i + i1 - RADY) * W + (j + j1 - RADX)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }
}

void smooth_2d_untiled(
    const int W,
    const int H,
    float* xs_in,
    float* xs_out
) {
    const int l = (2 * RADX + 1) * (2 * RADY + 1);

#pragma omp parallel for collapse(2)
    for(int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }
}

void smooth_2d_tiled(
    const uint64_t W,
    const uint64_t H,
    int TW,
    int TH,
    float* xs_in,
    float* xs_out
) {
    const uint64_t TILES_X = (W + TW - 1) / TW;
    const uint64_t TILES_Y = (H + TH - 1) / TH;
    const int l = (2 * RADX + 1) * (2 * RADY + 1);

    int imax, jmax;

#pragma omp parallel for collapse(2)
    for(int i = 0; i < H; i += TH) {
        for(int j = 0; j < W; j += TW) {
            imax = MIN(H, i+TH);
            jmax = MIN(W, j+TW);
            for(int ii = i; ii < imax; ii++) {
                for(int jj = j; jj < jmax; jj++) {
                    float acc = 0;
                    for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                        for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                            acc += xs_in[MIN(MAX(ii + i1 - RADY, 0), H-1) * W + MIN(MAX(jj + j1 - RADX, 0), W-1)];
                        }
                    }
                    xs_out[ii * W + jj] = acc / l;
                }
            }
        }
    }

}

void smooth_2d_tiled_partitioned(
    const uint64_t W,
    const uint64_t H,
    int TW,
    int TH,
    float* xs_in,
    float* xs_out
) {
    const uint64_t TILES_X = (W + TW - 1) / TW;
    const uint64_t TILES_Y = (H + TH - 1) / TH;
    const int l = (2 * RADX + 1) * (2 * RADY + 1);

    int imax, jmax;

    // Top edge
#pragma omp parallel for collapse(2)
    for(int i = 0; i < RADY; i++) {
        for(int j = 0; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

    // Bottom edge
#pragma omp parallel for collapse(2)
    for(int i = H-RADY; i < H; i++) {
        for(int j = 0; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }


#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i++) {
        // Left edge
        for(int j = 0; j < RADX; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i++) {
        // Right edge
        for(int j = W-RADX; j < W; j++) {
            float acc = 0;
            for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                    acc += xs_in[MIN(MAX(i + i1 - RADY, 0), H-1) * W + MIN(MAX(j + j1 - RADX, 0), W-1)];
                }
            }
            xs_out[i * W + j] = acc / l;
        }
    }

#pragma omp parallel for collapse(2)
    for(int i = RADY; i < H-RADY; i += TH) {
        for(int j = RADX; j < W-RADX; j += TW) {
            imax = MIN(H-RADY, i+TH);
            jmax = MIN(W-RADX, j+TW);
            for(int ii = i; ii < imax; ii++) {
                for(int jj = j; jj < jmax; jj++) {
                    float acc = 0;
                    for(int i1 = 0; i1 < (2 * RADY + 1); i1++) {
                        for(int j1 = 0; j1 < (2 * RADX + 1); j1++) {
                            acc += xs_in[(ii + i1 - RADY) * W + (jj + j1 - RADX)];
                        }
                    }
                    xs_out[ii * W + jj] = acc / l;
                }
            }
        }
    }

}

void measure_smooth_2d_untiled(
    const uint64_t W,
    const uint64_t H,
    float* xs_in,
    float* xs_out
) {
    // dry run
    smooth_2d_untiled(W, H, xs_in, xs_out);

    uint64_t start = start_measurement();

    //for(int i = 0; i < 3; i++) {
        smooth_2d_untiled(W, H, xs_in, xs_out);
    //}

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = W * H * ((RADX * 2 + 1) * (RADY * 2 + 1) + 1) * 4;

    char name[100];
    sprintf(name, "UNTILED WITH %d x %d STENCIL", RADX * 2 + 1, RADY * 2 + 1);
    report(name, accesses, elapsed, RESULT_NONE);
}

void measure_smooth_2d_untiled_partitioned(
    const uint64_t W,
    const uint64_t H,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {
    // dry run
    smooth_2d_untiled_partitioned(W, H, xs_in, xs_out);

    uint64_t start = start_measurement();

    //for(int i = 0; i < 3; i++) {
        smooth_2d_untiled_partitioned(W, H, xs_in, xs_out);
    //}

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = W * H * ((RADX * 2 + 1) * (RADY * 2 + 1) + 1) * 4;

    Result res = arraysEqual(W, H, xs_res, xs_out);

    char name[100];
    sprintf(name, "UNTILED P WITH %d x %d STENCIL", RADX * 2 + 1, RADY * 2 + 1);
    report(name, accesses, elapsed, res);
}

void measure_smooth_2d_tiled(
    const uint64_t W,
    const uint64_t H,
    int TW,
    int TH,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {
    // dry run
    //smooth_2d_tiled(W, H, TW, TH, xs_in, xs_out);

    uint64_t start = start_measurement();

    //for(int i = 0; i < 3; i++) {
        smooth_2d_tiled(W, H, TW, TH, xs_in, xs_out);
    //}

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = W * H * ((RADX * 2 + 1) * (RADY * 2 + 1) + 1) * 4;

    Result res = arraysEqual(W, H, xs_res, xs_out);

    char name[100];
    sprintf(name, "TILED %d x %d WITH %d x %d STENCIL", TW, TH, RADX * 2 + 1, RADY * 2 + 1);
    //report(name, elapsed/3, res);
    report(name, accesses, elapsed, res);
}

void measure_smooth_2d_tiled_partitioned(
    const uint64_t W,
    const uint64_t H,
    int TW,
    int TH,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {
    // dry run
    //smooth_2d_tiled_partitioned(W, H, TW, TH, xs_in, xs_out);

    uint64_t start = start_measurement();

    //for(int i = 0; i < 3; i++) {
        smooth_2d_tiled_partitioned(W, H, TW, TH, xs_in, xs_out);
    //}

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = W * H * ((RADX * 2 + 1) * (RADY * 2 + 1) + 1) * 4;

    Result res = arraysEqual(W, H, xs_res, xs_out);

    char name[100];
    sprintf(name, "TILED P %d x %d WITH %d x %d STENCIL", TW, TH, RADX * 2 + 1, RADY * 2 + 1);
    //report(name, elapsed/3, res);
    report(name, accesses, elapsed, res);
}

#endif // SMOOTH_2D_SEQ
