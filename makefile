CC=g++
CC+= -Wno-deprecated-declarations

CFLAGS = -I./ -I./eigen/ -std=c++14
LFLAGS = -lgsl -lgslcblas -lm -std=c++14

SRCDIR = ./src
OBJDIR = ./obj

object:
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $(SRCDIR)/*.cpp
	@mv -f ./*.o $(OBJDIR)

main: object
	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.exe
	@rm -rf $(OBJDIR)