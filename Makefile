CC = g++

C_FLAGS = -g -O3 -lm -fopenmp -lstdc++fs -std=c++17
C_FLAGS_DEBUG = -g -O0 -lm -fopenmp -lstdc++fs -std=c++17

CUR_SOURCE = ${wildcard *.cpp}
CUR_OBJS = ${patsubst %.cpp, %.o, $(CUR_SOURCE)}

ROOT_DIR=$(shell pwd)

BIN=IBM
BIN_DEBUG=IBM_debug

OBJS_DIR=build/release/obj
OBJS_DIR_DEBUG=build/debug/obj
BIN_DIR=build/release/bin
BIN_DIR_DEBUG=build/debug/bin

INCLUDE_DIR = -I$(ROOT_DIR)

all: $(CUR_OBJS) final_bin

final_bin:$(CUR_OBJS)
	$(CC) $(C_FLAGS) -o $(ROOT_DIR)/$(BIN_DIR)/$(BIN) $(OBJS_DIR)/*.o -lstdc++fs
	$(CC) $(C_FLAGS_DEBUG) -o $(ROOT_DIR)/$(BIN_DIR_DEBUG)/$(BIN_DEBUG) $(OBJS_DIR_DEBUG)/*.o -lstdc++fs

$(CUR_OBJS):%.o:%.cpp
	$(CC) $(INCLUDE_DIR) $(C_FLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@
	$(CC) $(INCLUDE_DIR) $(C_FLAGS_DEBUG) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR_DEBUG)/$@

.PHONY: clean

clean:
	rm -rf build
	mkdir build
	mkdir build/debug
	mkdir build/release
	mkdir ${OBJS_DIR}
	mkdir ${OBJS_DIR_DEBUG}
	mkdir ${BIN_DIR}
	mkdir ${BIN_DIR_DEBUG}