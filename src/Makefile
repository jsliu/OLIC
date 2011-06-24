# compiler
CC = g++

# compiler flags go here
INCLUDE = -I../include
CCFLAGS = -g -c
SHARE = -shared -L/usr/local/lib

# source files
SOURCE = BCF.o simBCF.o mat_utsil.o rand.o

all:main
main:$(SOURCE)
	$(CC) -g -o bcf $(SOURCE) 
BCF.o:BCF.cpp
	$(CC) $(INCLUDE) $(CCFLAGS) BCF.cpp -o BCF.o
simBCF.o:simBCF.cpp
	$(CC) $(INCLUDE) $(CCFLAGS) simBCF.cpp -o simBCF.o
mat_utsil.o:mat_utsil.cpp
	$(CC) $(INCLUDE) $(CCFLAGS) mat_utsil.cpp -o mat_utsil.o
rand.o:rand.cpp
	$(CC) $(INCLUDE) $(CCFLAGS) rand.cpp -o rand.o

clean:
	rm -rf *.o bcf
