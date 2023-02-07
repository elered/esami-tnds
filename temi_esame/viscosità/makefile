LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main : main.cc visco.cc randomgen.h visco.h
	g++ -o main main.cc visco.cc ${INCS} ${LIBS} 

clean:
	rm main.cc