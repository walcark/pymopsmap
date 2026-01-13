FC=gfortran
FFLAGS=-Ofast -ffree-line-length-none -lnetcdff -I/usr/include -Wall -std=f2008

SOURCEFILES=module_input.f90 module_aux_var.f90 module_results.f90 module_contributions.f90 module_dbindex.f90 mopsmap.f90 add_contribution.f90 init_shape.f90 init_size.f90 log_distr.f90 mod_gamma.f90 size_distr_interpolate.f90 calc_hygroscopic_growth.f90 write_output_ascii.f90 calc_optics.f90 integrate.f90 make_contributions.f90 get_nc_filename.f90 interpolate_linear.f90 matr.f90 read_dbindex.f90 write_output_netcdf.f90 init_wavelength_refr.f90 log_distr_ar.f90 mix_contributions.f90 read_input.f90
BINARY=../mopsmap

CREATE_INDEX_SOURCEFILES=module_input.f90 create_index.f90 create_index_check_file.f90 get_nc_filename.f90
CREATE_INDEX_BINARY=../create_index

mopsmap: $(SOURCEFILES)
	$(FC) $(SOURCEFILES) $(FFLAGS) -o $(BINARY); rm *.mod

create_index: $(CREATE_INDEX_SOURCEFILES)
	$(FC) $(CREATE_INDEX_SOURCEFILES) $(FFLAGS) -o $(CREATE_INDEX_BINARY); rm *.mod
