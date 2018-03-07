CC       = g++
CXXFLAGS = -std=c++0x -pedantic -Wall -Wextra

PROGRAM = inm

#TARGETS = ./INM.cpp ./min.cpp ./main.cpp
LISTEOBJ = grad.o INM.o min.o

all : $(PROGRAM)

objs : $(LISTEOBJ)

nvector.o : ../nvector.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

$(PROGRAM) : $(LISTEOBJ) main.cpp
	$(CC) $(CXXFLAGS) $^ -o $@

%.o : %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean :
	rm *.o ; rm $(PROGRAM)
