all: main

SOURCES = makefile\
	mapa/modul.f90\
	mapa/main.f90\

lnx: $(SOURCES)
	gfortran -c mapa/modul.f90
	gfortran -c mapa/main.f90
	gfortran -o main.out main.o modul.o
	rm main.o modul.o rutine.mod


win: $(SOURCES)
	gfortran -c mapa/modul.f90
	gfortran -c mapa/main.f90
	gfortran -o main.exe main.o modul.o
	del main.o modul.o rutine.mod
