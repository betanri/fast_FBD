CXX ?= c++

all: ffbd ffbd_sim

ffbd: ffbd.cpp
	$(CXX) -std=c++17 -O3 -o ffbd ffbd.cpp

ffbd_sim: ffbd_sim.cpp
	$(CXX) -std=c++17 -O2 -o ffbd_sim ffbd_sim.cpp

clean:
	rm -f ffbd ffbd_sim
