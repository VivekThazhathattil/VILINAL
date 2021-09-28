PROG_NAME := test_matrix
CC ?= gcc
SRC := ./src/*.c ./tests/*.c
OBJ_DIR := ./obj
OBJ_FILES := ./obj/*.o
CFLAGS := -Wall -Wextra -Werror -std=c99 -pedantic \
          -DPROG_NAME=\"$(PROG_NAME)\" -Iinclude/
ifdef DEBUG
	CFLAGS += -g
endif
ifeq ($(OS), Windows_NT)
	CFLAGS += -mwindows
endif
#LINK_FLAGS := $(shell pkg-config --libs gtk+-3.0 librsvg-2.0) -lcurl -pthread
build:
	$(CC) -g $(CFLAGS) -c $(SRC)
	mkdir -p $(OBJ_DIR) && mv ./*.o $(OBJ_DIR)
	$(CC) $(OBJ_FILES) -o $(PROG_NAME) -lm
#$(CC) $(OBJ_FILES) -o $(PROG_NAME) $(LINK_FLAGS)

clean:
	rm -f $(PROG_NAME)
	rm -rf $(OBJ_DIR)

format:
	clang-format -i */*.h */*.c

run: 
	chmod 755 $(PROG_NAME)
	./$(PROG_NAME)
