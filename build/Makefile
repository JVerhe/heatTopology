CXX=g++
CPPFLAGS=-std=c++20 -Wall -Wextra -Werror

all: adjointMethod

adjointMethod: ../src/adjointMethod.cpp
	$(CXX) $(CPPFLAGS) -o $@ $<

clean:
	find . -type f -executable -delete
	rm -f *.o