LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main : main.cc circuito.cc randomgen.h circuito.h
	g++ -o $@ main.cc circuito.cc ${INCS} ${LIBS} 

clean:
	rm main