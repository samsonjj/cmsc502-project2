program.exe : main.cpp tsp.cpp header.h
	mpic++ main.cpp -o program.exe -pthread -g -std=c++11

run:
	mpirun -np 26 ./program.exe
