

all:
	g++ -o acretee acretee.cc

debug:
	g++ -g -o acretee acretee.cc

clean:
	rm *.o 
	rm acretee
