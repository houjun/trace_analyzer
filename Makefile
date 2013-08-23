all:
	make main.c pattern.c -o analyzer -g

clean:
	rm -f *.o analyzer
