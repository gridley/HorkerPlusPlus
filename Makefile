OBJECTS = input_parsing.o position.o
DEBUG = -O0 -g
OPT = -O4 -ffast-math -fopenmp
FLAGS = $(OPT) -std=c++17

%.o: %.cc %.h
	g++ -c $(FLAGS) $< -o $@

"horker++": $(OBJECTS)
	g++ $(OBJECTS) $(FLAGS) horker.cc -o horker++ -lstdc++fs
