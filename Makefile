all: model

model: model.cpp
	g++ -o model model.cpp -litpp

.PHONY: clean
clean:
	$(RM) model
