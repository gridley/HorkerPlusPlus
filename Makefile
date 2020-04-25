OBJECTS = input_parsing.o
FLAGS = -O0 -g -std=c++17

%.o: %.cc %.h
	g++ -c $(FLAGS) $< -o $@

"horker++": $(OBJECTS)
	g++ $(OBJECTS) $(FLAGS) horker.cc -o horker++ -lopenmc
