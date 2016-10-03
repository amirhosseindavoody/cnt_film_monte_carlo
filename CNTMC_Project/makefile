CC=g++
CC+= -Wno-deprecated-declarations
CCFLAGS = -I./ -I../Eigen/ -std=c++11

SRCDIR = ./CNTMC_Project
OBJDIR = ./obj

main: 
	$(CC) $(CCFLAGS) -o $@.exe $(SRCDIR)/*.cpp

# Utility targets
.PHONY: clean
clean:
	@rm -f *.exe