all: 2d 3d

bench3d:
	./tilesize3d.sh 16 50 10
	make 3d
	./measure_3d_untiled.sh $(OUT_FILE)
	./measure_3d_tiled.sh $(OUT_FILE)
	./tilesize3d.sh 16 10 50
	make 3d
	./measure_3d_tiled.sh $(OUT_FILE)
	./tilesize3d.sh 10 16 50
	make 3d
	./measure_3d_tiled.sh $(OUT_FILE)
	./tilesize3d.sh 10 50 16
	make 3d
	./measure_3d_tiled.sh $(OUT_FILE)
	./tilesize3d.sh 50 10 16
	make 3d
	./measure_3d_tiled.sh $(OUT_FILE)
	./tilesize3d.sh 50 16 10
	make 3d
	./measure_3d_tiled.sh $(OUT_FILE)

bench2d:
	./stensize2d.sh 4 4
	make 2d
	./measure_2d_untiled.sh bench9x9-gbs.txt
	./measure_2d_tiled.sh bench9x9-gbs.txt 8 8
	./measure_2d_tiled.sh bench9x9-gbs.txt 24 24
	./measure_2d_tiled.sh bench9x9-gbs.txt 56 56
	./measure_2d_tiled.sh bench9x9-gbs.txt 120 120
	./measure_2d_tiled.sh bench9x9-gbs.txt 248 248
	./measure_2d_tiled.sh bench9x9-gbs.txt 504 504
	./measure_2d_tiled.sh bench9x9-gbs.txt 24 504
	./measure_2d_tiled.sh bench9x9-gbs.txt 504 24
	./measure_2d_tiled.sh bench9x9-gbs.txt 24 248
	./measure_2d_tiled.sh bench9x9-gbs.txt 248 24

2d: main_2d_untiled.c main_2d_tiled.c
	gcc -std=c99 -O3 -fopenmp -o main_2d_untiled_par main_2d_untiled.c
	gcc -std=c99 -O3 -fopenmp -o main_2d_tiled_par main_2d_tiled.c
	gcc -std=c99 -O3 -o main_2d_untiled_seq main_2d_untiled.c
	gcc -std=c99 -O3 -o main_2d_tiled_seq main_2d_tiled.c

3d: main_3d_untiled.c main_3d_tiled.c
	gcc -std=c99 -O3 -fopenmp -o main_3d_untiled_par main_3d_untiled.c
	gcc -std=c99 -O3 -o main_3d_untiled_seq main_3d_untiled.c
	gcc -std=c99 -O3 -fopenmp -o main_3d_tiled_par main_3d_tiled.c
	gcc -std=c99 -O3 -o main_3d_tiled_seq main_3d_tiled.c
