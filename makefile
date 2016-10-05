CC=g++
CC+= -Wno-deprecated-declarations
CCFLAGS = -I./ -I./eigen/ -std=c++11

SRCDIR = ./src
OBJDIR = ./obj

main: 
	$(CC) $(CCFLAGS) -o $@.exe $(SRCDIR)/*.cpp

# Utility targets
.PHONY: clean
clean:
	@rm -f *.exe