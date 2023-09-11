all: a.out


a.out: tests.o 
	g++ tests.o -o a.out

tests.o: tests.cpp 
	g++ -c tests.cpp  -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -O3

clean: 
	rm -rf *.o a.out


