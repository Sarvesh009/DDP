all: model model-debug

model: model.cpp
	g++ -o model model.cpp -litpp

model-debug: model.cpp
	g++ -ggdb -o model-debug model.cpp -litpp

.PHONY: clean
clean:
	$(RM) model model-debug
