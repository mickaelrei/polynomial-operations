BIN := ./bin
FLAGS := -Wall

all: $(BIN) $(BIN)/main

$(BIN):
	@if [ ! -d $(BIN) ]; then mkdir $(BIN); fi

$(BIN)/main: main.cpp polynomial.hpp
	$(CXX) main.cpp -o $(BIN)/main

clean:
	@if [ -d $(BIN) ]; then rm $(BIN)/*; fi
	@if [ -d $(BIN) ]; then rmdir $(BIN); fi