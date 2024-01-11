CC=g++
CFLAGS= -O3 -Wall -Wextra -pedantic -g -lm
SRCDIR=src
BINDIR=bin
SOURCES=$(wildcard $(SRCDIR)/*.cpp)
OBJECTS=$(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))
EXECUTABLE=$(BINDIR)/project

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE).cpp -o $(EXECUTABLE).x $(OBJECTS)

$(BINDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE).x
