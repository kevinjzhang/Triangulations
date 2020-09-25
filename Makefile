all: main.cc
	mpic++ -O3 -fopenmp -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
serial:
	mpic++ -O3 -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
debug:
	g++ -g -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
prof:
	g++ -pg -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
gprof:
	gprof triangulation gmon.out > analysis.txt
slurm:
	rm slurm*
	