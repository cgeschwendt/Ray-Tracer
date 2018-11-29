default: main 

valgrind:
	gcc -g rayTracer.c -o ray -lm

main:
	gcc rayTracer.c -o ray -lm
	
reference:
	./ray reference

custom:
	./ray custom

clean:
	rm -f ray reference.png custom.png
		 
