FC = gfortran
OPT = -g

#---------executables------------------

# Programme principal Gauss-Seidel
GS : GS.o MatHilbert.o CalculB.o matriceCarre.o diagonaleDom.o GaussSeidelFunction.o \
     inverse_matrice.o diagonale.o triangulaireInf.o triangulaireSupp.o \
     matrice_multiplication.o matrice_vect_multiplication.o
	$(FC) $(OPT) $^ -o GS

#---------objets------------------

GS.o : GS.f
	$(FC) $(OPT) -c GS.f

MatHilbert.o : MatHilbert.f
	$(FC) $(OPT) -c MatHilbert.f

CalculB.o : CalculB.f
	$(FC) $(OPT) -c CalculB.f

matriceCarre.o : matriceCarre.f
	$(FC) $(OPT) -c matriceCarre.f

diagonaleDom.o : diagonaleDom.f
	$(FC) $(OPT) -c diagonaleDom.f

GaussSeidelFunction.o : GaussSeidelFunction.f
	$(FC) $(OPT) -c GaussSeidelFunction.f

inverse_matrice.o : inverse_matrice.f
	$(FC) $(OPT) -c inverse_matrice.f

diagonale.o : diagonale.f
	$(FC) $(OPT) -c diagonale.f

triangulaireInf.o : triangulaireInf.f
	$(FC) $(OPT) -c triangulaireInf.f

triangulaireSupp.o : triangulaireSupp.f
	$(FC) $(OPT) -c triangulaireSupp.f

matrice_multiplication.o : matrice_multiplication.f
	$(FC) $(OPT) -c matrice_multiplication.f

matrice_vect_multiplication.o : matrice_vect_multiplication.f
	$(FC) $(OPT) -c matrice_vect_multiplication.f

#---------outils------------------

print_GS :
	cat *.f > print_GS

clean :
	rm -f GS *.o *.mod

