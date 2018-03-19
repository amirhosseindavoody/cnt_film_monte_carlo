CC=g++-7
CC+= -O3 -Wall

OS_NAME=$(shell uname -s)



ifeq ($(OS_NAME),Linux)
	CFLAGS = -I./ -I./eigen/ -std=c++17
	LFLAGS = -lm -larmadillo -Wl,-rpath=/home/amirhossein/anaconda3/lib/  -lstdc++fs -std=c++17
endif

ifeq ($(OS_NAME),Darwin)
	CFLAGS = -I./ -I./eigen/ -std=c++17
	LFLAGS = -lm -larmadillo  -lstdc++fs -std=c++17
endif

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