CXX	= g++
CXXFLAGS = -std=c++11 -O2


COMPASS: COMPASS.o Scores.o Node.o Tree.o Inference.o input.o
	${CXX} ${CXXFLAGS} -fopenmp -o COMPASS COMPASS.o Scores.o Node.o Tree.o Inference.o input.o


COMPASS.o: COMPASS.cpp Inference.h Tree.h input.h Scores.h Structures.h
	${CXX} ${CXXFLAGS} -fopenmp -c COMPASS.cpp
Scores.o: Scores.cpp Scores.h Structures.h
	${CXX} ${CXXFLAGS} -c Scores.cpp
Node.o:  Node.cpp Node.h Scores.h Structures.h
	${CXX} ${CXXFLAGS} -c Node.cpp
Tree.o: Tree.cpp Tree.h Node.h Scores.h Structures.h input.h
	${CXX} ${CXXFLAGS} -c Tree.cpp
Inference.o: Inference.cpp Inference.h Structures.h Scores.h Tree.h
	${CXX} ${CXXFLAGS} -c Inference.cpp
input.o: input.cpp input.h Structures.h Scores.h
	${CXX} ${CXXFLAGS} -c input.cpp 
clean:
	rm *.o