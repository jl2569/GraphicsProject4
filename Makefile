all: JSONParser.c
	gcc  JSONParser.c -o raycast

clean:
	rm -rf JSONParser *~
