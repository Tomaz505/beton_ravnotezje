all: main

SOURCES = makefile\
	modul.f90\
	main.f90\

lnx: $(SOURCES)
	gfortran -c modul.f90
	gfortran -c main.f90
	gfortran -o main.out main.o modul.o
	rm main.o modul.o rutine.mod


win: $(SOURCES)
	gfortran -c modul.f90
	gfortran -c main.f90
	gfortran -o main.exe main.o modul.o
	del main.o modul.o rutine.o
