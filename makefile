CC = g++
CFLAGS = -std=c++11 -Wall
LDFLAGS = -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LIBS = -lhdf5 -lhdf5_cpp

SRCS = main.cpp Parameters.cpp Utilities.cpp Meshgeneration.cpp Element.cpp Lagrangebasis.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = Vortex-Transport

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
	rm -r grid
	rm -r output