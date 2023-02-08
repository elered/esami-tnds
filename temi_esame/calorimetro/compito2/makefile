LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main : main.cc calorimetro.cc randomgen.h calorimetro.h
	g++ -o main main.cc calorimetro.cc ${INCS} ${LIBS} 

clean:
	rm main.cc