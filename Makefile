all: main.cc
	g++ -O3 -fopenmp -std=c++17 `regina-engine-config --cflags --libs` main.cc -o triangulation
slurm:
	rm slurm*
	