all:
	gcc main.c pattern.c -o analyzer

clean:
	rm -f analyzer
