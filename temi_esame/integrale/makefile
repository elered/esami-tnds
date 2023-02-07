LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

reticolo : reticolo.cpp funzioni.h FunzioneBase.h RandomGen.h EsperimentoReticolo.h
	g++ -o reticolo reticolo.cpp ${INCS} ${LIBS} 

integrale : integrale.cpp funzioni.h FunzioneBase.h Integral.h
	g++ -o integrale integrale.cpp ${INCS} ${LIBS} 

viscosita : viscosita.cpp funzioni.h FunzioneBase.h RandomGen.h Viscosita.h
	g++ -o viscosita viscosita.cpp ${INCS} ${LIBS}

calorimetro : calorimetro.cpp funzioni.h FunzioneBase.h RandomGen.h EquazioniDifferenziali.h
	g++ -o calorimetro calorimetro.cpp ${INCS} ${LIBS}

cspecifico : cspecifico.cpp funzioni.h FunzioneBase.h RandomGen.h EsperimentoCalorimetro.h
	g++ -o cspecifico cspecifico.cpp ${INCS} ${LIBS}

resistenza : resistenza.cpp funzioni.h FunzioneBase.h RandomGen.h EsperimentoResistenza.h
	g++ -o resistenza resistenza.cpp ${INCS} ${LIBS}

clean_reticolo:
	rm reticolo

clean_integrale:
	rm integrale

clean_viscosita:
	rm viscosita

clean_calorimetro:
	rm calorimetro

clean_cspecifico:
	rm cspecifico

clean_resistenza:
	rm resistenza