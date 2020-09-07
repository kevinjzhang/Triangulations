all: main.cc
	g++ -O3 -fopenmp -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
debug:
	g++ -g -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
prof:
	g++ -pg -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
gprof:
	gprof triangulation gmon.out > analysis.txt
slurm:
	rm slurm*
	