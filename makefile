CC=g++-7
CC+= -O3 -Wall

OS_NAME=$(shell uname -s)
# CC+= -g

ifeq ($(OS_NAME),Linux)
	CFLAGS = -I./ -I./eigen/ -std=c++17
	LFLAGS = -lm -larmadillo -lblas -llapack -Wl,-rpath=/home/amirhossein/anaconda3/lib/  -lstdc++fs -std=c++17
endif

ifeq ($(OS_NAME),Darwin)
	CFLAGS = -I./ -I./eigen/ -std=c++17
	LFLAGS = -lm -larmadillo  -lstdc++fs -std=c++17
endif

SRCDIR = ./src/*.cpp ./src/exciton_transfer/*.cpp ./src/discrete_forster/*.cpp ./src/helper/*.cpp
OBJDIR = ./obj


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

# OSFLAG 				:=
# ifeq ($(OS),Windows_NT)
# 	OSFLAG += -D WIN32
# 	ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
# 		OSFLAG += -D AMD64
# 	endif
# 	ifeq ($(PROCESSOR_ARCHITECTURE),x86)
# 		OSFLAG += -D IA32
# 	endif
# else
# 	UNAME_S := $(shell uname -s)
# 	ifeq ($(UNAME_S),Linux)
# 		OSFLAG += -D LINUX
# 	endif
# 	ifeq ($(UNAME_S),Darwin)
# 		OSFLAG += -D OSX
# 	endif
# 		UNAME_P := $(shell uname -p)
# 	ifeq ($(UNAME_P),x86_64)
# 		OSFLAG += -D AMD64
# 	endif
# 		ifneq ($(filter %86,$(UNAME_P)),)
# 	OSFLAG += -D IA32
# 		endif
# 	ifneq ($(filter arm%,$(UNAME_P)),)
# 		OSFLAG += -D ARM
# 	endif
# endif

# all:
# 	@echo $(OS_NAME)
# 	# @echo $(OSFLAG)