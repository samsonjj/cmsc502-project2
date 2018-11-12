program.exe : main.cpp tsp.cpp header.h
	mpic++ main.cpp -o program.exe -pthread -g -std=c++11

run:
	mpirun -np 25 ./program.exe
