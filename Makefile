OBJECTS = input_parsing.o
DEBUG = -O0 -g
OPT = -O4 -ffast-math -fopenmp
FLAGS = $(OPT) -std=c++17

%.o: %.cc %.h
	g++ -c $(FLAGS) $< -o $@

"horker++": $(OBJECTS) input_parsing.cc
	g++ $(OBJECTS) $(FLAGS) horker.cc -o horker++ -lopenmc -lstdc++fs
