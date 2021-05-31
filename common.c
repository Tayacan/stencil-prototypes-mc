#ifndef COMMON
#define COMMON
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define RESULT_EQUAL   0
#define RESULT_UNEQUAL 1
#define RESULT_NONE    2

typedef uint8_t Result;

// Validation
Result arraysEqual(
    const uint64_t W,
    const uint64_t H,
    float *a,
    float *b
) {
    for(uint64_t i = 0; i < H; i++) {
        for(uint64_t j = 0; j < W; j++) {
            if(fabs(a[i * W + j] - b[i * W + j]) > 0.01) {
                printf("(%lu, %lu): %8.2f != %8.2f\n", j, i, a[i * W + j], b[i * W + j]);
                return RESULT_UNEQUAL;
            }
        }
    }
    return RESULT_EQUAL;
}

Result arraysEqual3d(
    const uint64_t W,
    const uint64_t H,
    const uint64_t D,
    float *a,
    float *b
) {
    for(uint64_t k = 0; k < D; k++) {
        for(uint64_t i = 0; i < H; i++) {
            for(uint64_t j = 0; j < W; j++) {
                if(fabs(a[k * W * H + i * W + j] - b[k * W * H + i * W + j]) > 0.01) {
                    printf("(%lu, %lu, %lu): %8.2f != %8.2f\n", j, i, k, a[k * W * H + i * W + j], b[k * W * H + i * W + j]);
                    return RESULT_UNEQUAL;
                }
            }
        }
    }
    return RESULT_EQUAL;
}

// measurements
const uint64_t sec2usec = 1000000;

uint64_t start_measurement() {
    struct timeval time;
    gettimeofday(&time, NULL);

    // current time in microseconds
    return time.tv_usec + sec2usec * time.tv_sec;
}

uint64_t end_measurement(uint64_t start_usec) {
    const uint64_t end_usec = start_measurement();

    // measurement in microseconds
    return end_usec - start_usec;
}

void report(
    const char *name,
    const uint64_t accesses,
    const float elapsed,
    const Result result
) {
    const char* valid;

    if (result == RESULT_EQUAL) {
        valid = "VALID";
    } else if (result == RESULT_UNEQUAL) {
        valid = "INVALID";
    } else if (result == RESULT_NONE) {
        valid = "-------";
    }
    const float gbs = (accesses / (1000 * 1000 * 1000)) / (elapsed / 1000000);
    printf("%40s: %7s - RAN FOR %10.2f usec (%8.2f ms) on average (%8.2f GB/s)\n", name, valid, elapsed, elapsed / 1000.0, gbs);
}

float rand_float(float min, float max) {
    float range = max - min;
    float r = ((float)rand()) / ((float)RAND_MAX/range);
    return r + min;
}

void write_random_1d(
    const uint64_t N,
    float* dst
) {
    for(uint64_t i = 0; i < N; i++) {
        dst[i] = rand_float(-100.0, 100.0);
    }
}

void write_random_2d(
    const uint64_t W,
    const uint64_t H,
    float** dst
) {
    for(uint64_t i = 0; i < H; i++) {
        write_random_1d(W, dst[i]);
    }
}
#endif // COMMON
