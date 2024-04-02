# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall

# Source files
SRCS = main.cpp Parameters.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = Vortex-Transport 

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Rule to build object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

