LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main : main.cc resistenza.cc randomgen.h resistenza.h
	g++ -o $@ main.cc resistenza.cc ${INCS} ${LIBS} 

clean:
	rm main