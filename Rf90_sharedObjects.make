#makefile for asymptotic theory calculations. this creates
#shared objects between R and Fortran
OBJS=	Fx_Eval.o 	\
	bisection.o 	\
	mc_uncertainty.o

FC= gfortran 

mc_uncertainty.so: $(OBJS)
	$(FC) -shared -o mc_uncertainty.so $(OBJS)

Fx_Eval.o:Fx_Eval.f90
	$(FC) -c Fx_Eval.f90

bisection.o:bisection.f90
	$(FC) -c bisection.f90

mc_uncertainty.o:mc_uncertainty.f90 
	$(FC) -c mc_uncertainty.f90

clean:
	rm $(OBJS) 

#End of makefile