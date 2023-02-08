LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main : main.cc integrale.h funzione_base.h 
	g++ -o main main.cc ${INCS} ${LIBS} 

clean:
	rm main