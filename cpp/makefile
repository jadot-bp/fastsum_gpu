CC=nvc++
CFLAGS= -Wall -Wextra -pedantic -lm -mp=multicore

SRCDIR=src
BINDIR=bin
SOURCES=$(wildcard $(SRCDIR)/*.cpp)
OBJECTS=$(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))
EXECUTABLE=$(BINDIR)/project

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE).cpp -o $(EXECUTABLE)_cpp.x $(OBJECTS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)_cpp.x
