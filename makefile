
OS_NAME=$(shell uname -s)

ifeq ($(OS_NAME),Linux)
	CC=g++-7
	CFLAGS = -std=c++17
	LFLAGS = -lm -larmadillo -Wl,-rpath=/home/amirhossein/anaconda3/lib/ -lpthread -lstdc++fs -std=c++17
endif

ifeq ($(OS_NAME),Darwin)
	CC=g++-8
	CFLAGS = -std=c++17
	LFLAGS = -lm -larmadillo  -lstdc++fs -std=c++17
endif

CC+= -O3 -Wall
# CC+= -O0 -g -Wall # for debuging using gdb

SRCDIR = ./src/*.cpp ./src/exciton_transfer/*.cpp ./src/discrete_forster/*.cpp ./src/helper/*.cpp
OBJDIR = ./obj

BUILD_DIR ?= ./build


object:
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $(SRCDIR)
	@mv -f ./*.o $(OBJDIR)

main: object
	$(CC) -o $@.exe $(OBJDIR)/*.o $(LFLAGS)

# Utility targets
.PHONY: clean
clean:
	@rm -f *.o *.exe
	@rm -rf $(OBJDIR)