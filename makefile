CC=gcc
CFLAGS= -std=c99 -Wall -Wextra -pedantic -g -lm
SRCDIR=src
BINDIR=bin
SOURCES=$(wildcard $(SRCDIR)/*.c)
OBJECTS=$(patsubst $(SRCDIR)/%.c,$(BINDIR)/%.o,$(SOURCES))
EXECUTABLE=$(BINDIR)/project

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE).c -o $(EXECUTABLE).x $(OBJECTS)

$(BINDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE).x
