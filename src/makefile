# Start by clearing the list of suffixes
.SUFFIXES:

# Default suffix list
.SUFFIXES: .f90 .f .o


# List of object files
OBJS =  modules.o              \
        main.o			\
        read_input.o            \
        gen_mesh.o              \
        ext_potential.o         \
	calc_rs.o		\
        starting_density.o \
        fillup_electrons.o \
        rho_construction.o  \
        check_rho.o       \
	solve_3dhamilton.o \
	nonzero_sparse_counter.o \
	gen_hamilton.o		\
	gen_kinetic_hamilton.o \
	dvd_setup.o		\
	dvdson.o                \
	mat_vec_multi1.o        \
	dvd_err_chk.o		\
	write_eigenvecup_for_vmd.o \
        write_eigenvecdn_for_vmd.o \
        poisson_allocate_deallocate.o \
        poisson.o			\
        exchange.o   \
        correlation.o \
	charge_density.o \
	double_grid_trans.o \
	fft_grid_trans.o      \
	fft.o             \
	gen_hartree_pot.o \
        analytic_pot.o \
        charge_den_mixing.o \
        KE_spin.o \
        del2Psi.o \
        kinetic_en.o\
        disaster_management.o \
	analytic_dft_energy.o \
	del2.o                \
	charge_den_for_vmd.o \
        write_for_plot.o


EXECNAME =../bin/3d_dft.x
FC = ifort 

CFLAGS = -r8  -O3 -i8 
LFLAGS = -L /opt/intel/mkl/9.1.023/lib/32/ \
	 -L .  \
	-lmkl_lapack -lmkl -lpthread -lguide \
	-lfftw3 -lm  -lsvml \
	-lxc

	   
# Default rule for making a .o file from a .f90 file
.f90.o:
	$(FC) -c $(CFLAGS) $<

# Perform a complete compilation (**default action**)
all: $(EXECNAME)

$(EXECNAME): $(OBJS)
	$(FC) -o $(EXECNAME) $(OBJS) $(LFLAGS)


clean:
	rm -f **.o global_var.mod  poisson_var.mod kind_param.mod
distclean:
	rm -f *.o global_var.mod  poisson_var.mod kind_param.mod *.dat *.VASP_CHGCAR *.info *.x	

