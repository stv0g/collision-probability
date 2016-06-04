CXX = clang++

CXXFLAGS = -I/opt/local/include/eigen3/

.PHONY: plot all clean

all: collision_probabilty

collision_probabilty: collision_probabilty.cpp
	$(CXX) $(CXXFLAGS) collision_probabilty.cpp -o collision_probabilty
	
samples.dat: collision_probabilty
	./collision_probabilty > samples.dat
	
plot: samples.dat
	gnuplot -c gnuplot.script

clean:
	rm -f collision_probability
	rm -f samples.dat