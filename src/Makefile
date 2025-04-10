ABDIR   = ${MYDIR}/srcSphereAB
IODIR   = ${DATAREADINDICESDIR}
GLDIR   = ${MYDIR}/srcGlow

SEARCH_EXTRA = -I${ABDIR} -I${IODIR} -I${EMPIRICALUADIR} -I${EMPIRICALIEDIR}

-include ../srcMake/Makefile.conf

-include Makefile.DEPEND
-include Makefile.RULES

MODULES = \
	ModSize.o\
	ModTime.o\
	ModConstants.o\
	ModPlanet.o\
	ModGITM.o\
	ModSources.o\
	ModVertical.o\
	ModInputs.o\
	ModSphereInterface.o\
	ModTimingGITM.o\
	ModRates.o\
	ModEUV.o\
	ModElectrodynamics.o\
	ModSatellites.o\
	ModLimiter.o\
	ModTides.o\
	ModNewell.o\
	ModMagTrace.o\
	ModChemistry.o\
	ModGITMImplicit.o\
	ModUser.o

OBJECTS = \
	planet.o\
	satellites.o\
	calc_electron_temperature.o\
	aurora.o\
	get_potential.o\
	calc_efield.o\
	set_inputs.o\
	init_energy_deposition.o\
	get_glow.o\
	init_grid.o\
	time_routines.o\
	vertical_solver.o\
	library.o\
	stretch_grid.o\
	init_mpi.o\
	restart.o\
	initialize.o\
	init_msis.o\
	init_iri.o\
	init_b0.o\
	init_altitude.o\
	finalize.o\
	calc_neutral_friction.o\
	lubksb.o\
	ludcmp.o\
	calc_electrodynamics.o\
	divergence.o\
	calc_pressure.o\
	calc_rates.o\
	calc_timestep.o\
	exchange_messages_sphere.o\
	timing.o\
	advance.o\
	advance_vertical.o\
	advance_horizontal.o\
	set_vertical_bcs.o\
	set_horizontal_bcs.o\
	output.o\
	output_common.o\
	logfile.o\
	calc_conduction.o\
	calc_physics.o\
	calc_avesza.o\
	calc_euv.o\
	calc_chemistry.o\
	calc_ion_v.o\
	chapman_new.o\
	stop_file.o\
	read_inputs.o\
	add_sources.o\
	calc_sources.o\
	apex.o\
	write_output.o\
	gradient.o\
	tides.o\
	get_location.o\
	user.o

OBJECTS_EXE = main.o ${OBJECTS} ${MODULES}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS_EXE}

EXE = ${BINDIR}/GITM.exe

GITM:	DEPEND
	@make ${EXE}
	@make POST
	@echo "${EXE} has been created"

${EXE}:	${OBJECTS_EXE}
	@(if [ ${OS} = "Darwin" ]; then make GITM_for_Darwin; \
	else                            make GITM_general;     fi)

GITM_general:
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} \
		-L${ABDIR} -lSphere \
		-L${LIBDIR} -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE \
		-L${GLDIR} -lGLOW ${LBLAS} ${LflagMpi}

GITM_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${ABDIR}/libSphere.a; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libEMPIRICALUA.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a; \
		ar -x ${GLDIR}/libGLOW.a
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${LflagMpi}
	rm -rf Tmp_

POST:	PostProcess.o
	${LINK.f90} ${Lflag} -o ${BINDIR}/PostProcess.exe PostProcess.o


MY_LIB = libGITM.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBPREV} ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

clean: cleanfiles

distclean: clean
	rm -f 	ModSize.f90 \
		ModPlanet.f90 \
		ModChemistry.f90 \
		ModRates.f90 \
		aurora.f90 \
		calc_chemistry.f90 \
		calc_electron_temperature.f90 \
		calc_ion_density.f90 \
		calc_neutral_density.f90 \
		calc_rates.f90 \
		init_iri.f90 \
		init_msis.f90 \
		output.f90 \
		planet.f90 \
		set_vertical_bcs.f90
