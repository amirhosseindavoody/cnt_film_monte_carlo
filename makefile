CC=g++-7
CC+= -O3 -Wall
# CC+= -g

CFLAGS = -I./ -I./eigen/ -std=c++17
LFLAGS = -lm -larmadillo -lstdc++fs -std=c++17

# SRCDIR = ./src/*.cpp ./src/*/*.cpp
SRCDIR = ./src/*.cpp ./src/exciton_transfer/*.cpp ./src/discrete_forster/*.cpp ./src/helper/*.cpp
# SRCDIR = ./src
OBJDIR = ./obj

# object:
# 	@mkdir -p $(OBJDIR)
# 	$(CC) $(CFLAGS) -c $(SRCDIR)/*.cpp
# 	@mv -f ./*.o $(OBJDIR)

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


# # Look at this page for more
# # https://stackoverflow.com/questions/8654548/using-make-to-move-o-files-to-a-separate-directory
#
# CC=g++-7
# CC+= -O3
#
# CPPFLAGS = -I./ -I./eigen/ -std=c++17
# LDFLAGS = -lm -lstdc++fs -std=c++17
#
# SRCDIR = ./src
# OBJDIR = ./obj
#
# OBJ = $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
# TARGET = main
#
# .PHONY: all clean
#
# all: $(OBJDIR) $(TARGET)
#
# $(OBJDIR):
# 	mkdir $(OBJDIR)
#
# $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
# 	$(CC) $(CPPFLAGS) -c $< -o $@
#
# $(TARGET): $(OBJ)
# 	$(CC) $(LDFLAGS) -o $@ $^
#
# clean:
# 	@rm -f $(TARGET) $(wildcard *.o)
# 	@rm -rf $(OBJDIR)
