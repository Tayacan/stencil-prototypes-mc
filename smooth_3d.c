#include <math.h>
#include "common.c"

#define RAD 1

void smooth_3d_untiled(
    const int W,
    const int H,
    const int D,
    float* xs_in,
    float* xs_out
) {
    const int side = 2 * RAD + 1;
    const int l = pow(side, 3);
#pragma omp parallel for collapse(3)
    for (int k = 0; k < D; k++) {
        for(int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
}

void smooth_3d_untiled_partitioned(
    const int W,
    const int H,
    const int D,
    float* xs_in,
    float* xs_out
) {
    const int side = 2 * RAD + 1;
    const int l = pow(side, 3);
#pragma omp parallel for collapse(3)
    for (int k = 0; k < RAD; k++) {
        for(int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = D-RAD; k < D; k++) {
        for(int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = 0; i < RAD; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = H-RAD; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = RAD; i < H-RAD; i++) {
            for(int j = 0; j < RAD; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = RAD; i < H-RAD; i++) {
            for(int j = W-RAD; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = RAD; i < H-RAD; i++) {
            for(int j = RAD; j < W-RAD; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[(k + k1 - RAD) * W * H +
                                         (i + i1 - RAD) * W +
                                         (j + j1 - RAD)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
}

void smooth_3d_tiled(
    const int W,
    const int H,
    const int D,
    int TW,
    int TH,
    int TD,
    float* xs_in,
    float* xs_out
) {
    const int side = 2 * RAD + 1;
    const int l = pow(side, 3);
    int imax, jmax, kmax;
#pragma omp parallel for collapse(3)
    for (int k = 0; k < D; k += TD) {
        for(int i = 0; i < H; i += TH) {
            for(int j = 0; j < W; j += TW) {
                imax = MIN(H, i+TH);
                jmax = MIN(W, j+TW);
                kmax = MIN(D, k+TD);

                for(int kk = k; kk < kmax; kk++) {
                    for(int ii = i; ii < imax; ii++) {
                        for(int jj = j; jj < jmax; jj++) {
                            float acc = 0;
                            for(int k1 = 0; k1 < side; k1++) {
                                for(int i1 = 0; i1 < side; i1++) {
                                    for(int j1 = 0; j1 < side; j1++) {
                                        acc += xs_in[MIN(MAX(kk + k1 - RAD, 0), D-1) * W * H +
                                                     MIN(MAX(ii + i1 - RAD, 0), H-1) * W +
                                                     MIN(MAX(jj + j1 - RAD, 0), W-1)];
                                    }
                                }
                            }
                            xs_out[kk * W * H + ii * W + jj] = acc / l;
                        }
                    }
                }
            }
        }
    }
}

void smooth_3d_tiled_partitioned(
    const int W,
    const int H,
    const int D,
    int TW,
    int TH,
    int TD,
    float* xs_in,
    float* xs_out
) {
    const int side = 2 * RAD + 1;
    const int l = pow(side, 3);

#pragma omp parallel for collapse(3)
    for (int k = 0; k < RAD; k++) {
        for(int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = D-RAD; k < D; k++) {
        for(int i = 0; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = 0; i < RAD; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = H-RAD; i < H; i++) {
            for(int j = 0; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = RAD; i < H-RAD; i++) {
            for(int j = 0; j < RAD; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k++) {
        for(int i = RAD; i < H-RAD; i++) {
            for(int j = W-RAD; j < W; j++) {
                float acc = 0;
                for(int k1 = 0; k1 < side; k1++) {
                    for(int i1 = 0; i1 < side; i1++) {
                        for(int j1 = 0; j1 < side; j1++) {
                            acc += xs_in[MIN(MAX(k + k1 - RAD, 0), D-1) * W * H +
                                         MIN(MAX(i + i1 - RAD, 0), H-1) * W +
                                         MIN(MAX(j + j1 - RAD, 0), W-1)];
                        }
                    }
                }
                xs_out[k * W * H + i * W + j] = acc / l;
            }
        }
    }

    int imax, jmax, kmax;
#pragma omp parallel for collapse(3)
    for (int k = RAD; k < D-RAD; k += TD) {
        for(int i = RAD; i < H-RAD; i += TH) {
            for(int j = RAD; j < W-RAD; j += TW) {
                imax = MIN(H-RAD, i+TH);
                jmax = MIN(W-RAD, j+TW);
                kmax = MIN(D-RAD, k+TD);

                for(int kk = k; kk < kmax; kk++) {
                    for(int ii = i; ii < imax; ii++) {
                        for(int jj = j; jj < jmax; jj++) {
                            float acc = 0;
                            for(int k1 = 0; k1 < side; k1++) {
                                for(int i1 = 0; i1 < side; i1++) {
                                    for(int j1 = 0; j1 < side; j1++) {
                                        acc += xs_in[(kk + k1 - RAD) * W * H +
                                                     (ii + i1 - RAD) * W +
                                                     (jj + j1 - RAD)];
                                    }
                                }
                            }
                            xs_out[kk * W * H + ii * W + jj] = acc / l;
                        }
                    }
                }
            }
        }
    }
}

void measure_smooth_3d_untiled(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    float* xs_in,
    float* xs_out
) {

    uint64_t start = start_measurement();

    smooth_3d_untiled(W, H, D, xs_in, xs_out);

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = pow(2 * RAD + 1, 3) * W * H * D * 4;

    report("UNTILED", accesses, elapsed, RESULT_NONE);
}

void measure_smooth_3d_tiled(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    int TW,
    int TH,
    int TD,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {

    uint64_t start = start_measurement();

    smooth_3d_tiled(W, H, D, TW, TH, TD, xs_in, xs_out);

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = pow(2 * RAD + 1, 3) * W * H * D * 4;
    Result res = arraysEqual3d(W, H, D, xs_res, xs_out);

    char name[100];
    sprintf(name, "TILED %d x %d x %d", TW, TH, TD);
    report(name, accesses, elapsed, res);
}

void measure_smooth_3d_untiled_partitioned(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {

    uint64_t start = start_measurement();

    smooth_3d_untiled_partitioned(W, H, D, xs_in, xs_out);

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = pow(2 * RAD + 1, 3) * W * H * D * 4;
    Result res = arraysEqual3d(W, H, D, xs_res, xs_out);

    report("UNTILED P", accesses, elapsed, res);
}

void measure_smooth_3d_tiled_partitioned(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    int TW,
    int TH,
    int TD,
    float* xs_in,
    float* xs_out,
    float* xs_res
) {

    uint64_t start = start_measurement();

    smooth_3d_tiled_partitioned(W, H, D, TW, TH, TD, xs_in, xs_out);

    uint64_t elapsed = end_measurement(start);

    const uint64_t accesses = pow(2 * RAD + 1, 3) * W * H * D * 4;
    Result res = arraysEqual3d(W, H, D, xs_res, xs_out);

    char name[100];
    sprintf(name, "TILED P %d x %d x %d", TW, TH, TD);
    report(name, accesses, elapsed, res);
}
