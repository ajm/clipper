CC := gcc
CXX := g++

CFLAGS := -g -O2 -Wall -pedantic -std=c99
CXXFLAGS := -g -O2
LDFLAGS :=
INCLUDES := -I. 
LIBS := -lm

OBJECTS	:= \
    clipper.o

clipper: $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) clipper

