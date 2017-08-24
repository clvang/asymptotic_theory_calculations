#makefile for asymptotic theory calculations
OBJS=	Fx_Eval.o 	\
	asymTheoryDiffusivity.o 	\
	bisection.o	\
	linspace.o	\
	partialF_partial_x.o	\
	partialF_partial_dY.o	\
	partialF_partial_dc.o	\
	partialF_partial_do.o	\
	uncertainty_diffusivity.o

FC= gfortran 

asymTheoryDiffusivity: $(OBJS)
	$(FC) -o asymTheoryDiffusivity $(OBJS)

Fx_Eval.o:Fx_Eval.f90
	$(FC) -c Fx_Eval.f90

asymTheoryDiffusivity.o:asymTheoryDiffusivity.f90
	$(FC) -c asymTheoryDiffusivity.f90

bisection.o:bisection.f90
	$(FC) -c bisection.f90

linspace.o:linspace.f90
	$(FC) -c linspace.f90

partialF_partial_x.o:partialF_partial_x.f90
	$(FC) -c partialF_partial_x.f90

partialF_partial_dY.o:partialF_partial_dY.f90
	$(FC) -c partialF_partial_dY.f90

partialF_partial_dc.o:partialF_partial_dc.f90
	$(FC) -c partialF_partial_dc.f90

partialF_partial_do.o:partialF_partial_do.f90
	$(FC) -c partialF_partial_do.f90

uncertainty_diffusivity.o:uncertainty_diffusivity.f90
	$(FC) -c uncertainty_diffusivity.f90

clean:
	rm $(OBJS)

#End of makefile