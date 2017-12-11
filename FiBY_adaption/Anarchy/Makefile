#============================================================================
# Basic operation mode of code
#============================================================================

# Set this if you want to have periodic boundary conditions.
OPT += -DPERIODIC

# Set this if particles of different type have different softening.
#OPT += -DUNEQUALSOFTENINGS

# Set this to swap endianess of ICs.
#OPT += -DAUTO_SWAP_ENDIAN_READIC

#OPT += -DNO_ISEND_IRECV_IN_DOMAIN

#OPT += -DNOCALLSOFSYSTEM


#============================================================================
# TreePM Options
#============================================================================

# This enables the TreePM method, i.e. the long-range force is
# computed with a PM-algoritthm, and the short range force with the
# tree. The parameter has to be set to the size of the mesh that
# should be used, (e.g. 64, 96, 128, etc). The mesh dimensions need
# not necessarily be a power of two.  NOTE: if the simulation is not
# in a periodic box, then a FFT method for vacuum boundaries is
# employed, using a mesh with dimension twice that specified by PMGRID.
#OPT += -DPMGRID=512

# This values are set by default in allvars.h. Use the Makefile flags
# only if you want to change their default value.
#OPT += -DASMTH=1.25
#OPT += -DRCUT=4.5

# If this option is set (will only work together with PMGRID), then
# the long range force is computed in two stages: One Fourier-grid is
# used to cover the whole simulation volume, allowing the computation
# of the large-scale force. A second Fourier mesh is placed on the
# region occupied by "high-resolution" particles, allowing the
# computation of an intermediate scale force. Finally, the force on
# very small scales is supplemented by the tree. This procedure can be
# useful for "zoom-simulations", where the majority of particles (the
# high-res particles) are occupying only a small fraction of the
# volume. To activate this option, the parameter needs to be set to an
# integer that encodes the particle types that represent the high-res
# particles in the form of a bit mask. For example, if types 0, 1, and
# 4 form the high-res particles, set the parameter to
# PLACEHIGHRESREGION=1+2+16. The spatial region covered by the
# high-res grid is determined automatically from the initial
# conditions. NOTE: If a periodic box is used, the high-res zone may
# not intersect the box boundaries.
#OPT += -DPLACEHIGHRESREGION=3

# The spatial region covered by the high-res zone has a fixed size
# during the simulation, which initially is set to the smallest region
# that encompasses all high-res particles. Normally, the simulation
# will be interrupted, if high-res particles leave this region in the
# course of the run. However, by setting this parameter to a value
# larger than one, the high-res region can be expanded. For example,
# setting it to 1.4 will enlarge its side-length by 40% (it remains
# centered on the high-res particles). Hence, with such a setting, the
# high-res region may expand or move by a limited amount. The code
# will be able to continue even if high-res particles leave the
# initial high-res grid. In this case, the code will update the size
# and position of the grid that is placed onto the high-resolution
# region automatically. To prevent that this potentially happens every
# single PM step, one should nevertheless assign a value slightly
# larger than 1 to ENLARGEREGION.
#OPT += -DENLARGEREGION=1.2

#OPT += -DONLY_PM


#============================================================================
# Multi-Domain and Top-Level Tree options
#============================================================================

#OPT += -DMULTIPLEDOMAINS=8
#OPT += -DTOPNODEFACTOR=3.0


#============================================================================
# Single/Double Precision
#============================================================================

# If not set the code runs in single precision. If set to 1 all arrays
# are in double precision. If set to 2 only position, velocity and
# entropy are in double precision. Note that output files are
# nevertheless written by converting to single precision.
OPT += -DDOUBLEPRECISION=2

#OPT += -DDOUBLEPRECISION_FFTW

# If this is set, the fftw-header/libraries are accessed without type
# prefix (adopting whatever was chosen as default at compile of
# fftw). Otherwise, the type prefix 'd' for double is used.
#OPT += -DNOTYPEPREFIX_FFTW

# Enables round off reduction in particle sums. If DOUBLEPRECISION is
# set, these sums are done in 'long double'. If single precision is
# used, they are done in 'double'. This should in principle allow to
# make computations *exactly* invariant to different numbers of CPUs.
#OPT += -DFLTROUNDOFFREDUCTION

# When this is set, a software implementation of 128bit double-double
# addition is used, implemented as a c++ class. Hence, this option
# requires compilation with a c++ compiler.
#OPT += -DSOFTDOUBLEDOUBLE


#============================================================================
# On the fly FOF group finder 
#============================================================================

# Enable FOF code
#OPT += -DFOF

# Output FOF file every time the catalogues are computed
#OPT += -DFOF_OUTPUT

# Add extra OWLS specific arrays to FOF output files
#OPT += -DFOF_OUTPUT_OWLS_ARRAYS

# Add extra ParticleID  arrays to FOF output files
#OPT += -DFOF_OUTPUT_PID_ARRAY

# 2^type for the primary particle type to be used in the FOF calculation
#OPT += -DFOF_PRIMARY_LINK_TYPES=2

# 2^type for the types linked to nearest primaries
#OPT += -DFOF_SECONDARY_LINK_TYPES=49


#============================================================================
# AGN implementation
#============================================================================

# Enable black holes
#OPT += -DBLACK_HOLES

# Enable energy feedback from black holes
#OPT += -DBH_THERMALFEEDBACK

# Allow BHs to swallow gas particles when M_BH > M_PART
#OPT += -DSWALLOWGAS

# Enable HUGE BH logfiles:
#OPT += -DVERBOSE_LOGFILES

# Enable extra debugging variables for BH
#OPT += -DBH_DEBUG

# Enable repositioning low mass BHs on the potential minimum
#OPT += -DREPOSITION_ON_POTMIN


#============================================================================
# Things that are always recommended
#============================================================================

# This is a tuning option. When set, the code will bring the particles
# after each domain decomposition into Peano-Hilbert order. This
# improves cache utilizationand performance.
OPT += -DPEANOHILBERT

# If set, a wallclock timer is used by the code to measure internal
# time consumption (see cpu-log file).  Otherwise a timer that
# measures consumed processorticks is used.
OPT += -DWALLCLOCK

OPT += -DCPUSPEEDADJUSTMENT

# This replaces qsort() with Volker's sorting routine (much faster
# especially on AIX).
OPT += -DMYSORT


#============================================================================
# Things for special behaviour
#============================================================================

# This enables an internal memory handler (most recently allocated
# block needs to be freed first).
OPT += -DPEDANTIC_MEMORY_HANDLER

# This should be set to define the memory ceiling in MByte.
#OPT += -DPEDANTIC_MEMORY_CEILING=1800

# This is an alternative to PEDANTIC_MEMORY_CEILING and limits only
# the maximum memory used in the domain decomposition. Ordinary
# malloc() is used for all allocations, realloc() is not used.
#OPT += -DDOMAIN_MEMORY_CEILING=1600

#OPT += -DFORCE_DOMAIN_DECOMPOSITION #<------------

#OPT += -DMEMDEBUG
OPT += -DNOGRAVITY
#OPT += -DNOACCEL
#OPT += -DNOISMPRESSURE
OPT += -DNOVISCOSITYLIMITER

# If this is not set, the tree construction will succeed even when
# there are a few particles at identical locations. This is done by
# `rerouting' particles once the node-size has fallen below 1.0e-3 of
# the softening length. When this option is activated, this will be
# surpressed and the tree construction will always fail if there are
# particles at extremely close coordinates.
#OPT += -DNOTREERND

# If this is activated, the code will not terminate when the timestep
# falls below the value of MinSizeTimestep specified in the
# parameterfile. This is useful for runs where one wants to enforce a
# constant timestep for all particles. This can be done by activating
# this option, and by setting Min- and MaxSizeTimestep to an equal
# value.
#OPT += -DNOSTOP_WHEN_BELOW_MINTIMESTEP

# When this is set, the long-range timestep for the PM force
# computation is always determined by MaxSizeTimeStep. Otherwise, it
# is set to the minimum of MaxSizeTimeStep and the timestep obtained
# for the maximum long-range force with an effective softening scale
# equal to thePM smoothing-scale.
#OPT += -DNOPMSTEPADJUSTMENT

#OPT += -DNO_TREEDATA_IN_RESTART

# Add potential of an isothermal sphere.
#OPT += -DISOTHERM=200

#OPT += -DCOMPUTE_POTENTIAL_ENERGY
#OPT += -DALLOWEXTRAPARAMS
#OPT += -DLONGIDS

# If the gas distribution is spatially very different from collisionless
# particles, this helps to avoid problems in the domain decomposition
#OPT += -DINHOMOG_GASDISTR_HINT

# These options can be used together with PERIODIC and NOGRAVITY only.
# When set, the options define numerical factors that can be used to
# distorts the periodic simulation cube into a parallelepiped of
# arbitrary aspect ratio. This can be useful for idealized SPH tests.
#OPT += -DLONG_X=60
#OPT += -DLONG_Y=5
#OPT += -DLONG_Z=0.2

# This effectively switches of one dimension in SPH, i.e. the code
# follows only 2d hydrodynamics in the xy-, yz-, or xz-plane. This
# only works with NOGRAVITY, and if all coordinates of the third axis
# are exactly equal. Can be useful for idealizedSPH tests.
OPT += -DTWODIMS

# If this is set, particles with a particle-ID equal to zero do not
# receive any SPH acceleration. This can be useful for idealized SPH
# tests, where these particles represent fixed "walls".
#OPT += -DSPH_BND_PARTICLES

# Reads hsml from IC file.
#OPT += -DREAD_HSML

# Allows variable softening length for gas particles. This option
# require UNEQUALSOFTENINGLENGTH to be set.
#OPT += -DADAPTIVE_GRAVSOFT_FORGAS

# In case variable gravitational softening for gas is activated, this
# option couples the softening length to the SPH smoothing length.
#OPT += -DADAPTIVE_GRAVSOFT_FORGAS_HSML

# Splits type=1 particles in cosmological runs into gas and dark matter.
#OPT += -DGENERATE_GAS_IN_ICS

# Use the "Relative Pressure SPH" discretization of
# Abel (2010; arXiv: 1003.0937
#OPT += -DRPSPH


#============================================================================
# Time integration options
#============================================================================

#OPT += -DALTERNATIVE_VISCOUS_TIMESTEP #<------------


#============================================================================
# Output options
#============================================================================

# This will force the code to compute gravitational potentials for all
# particles each time a snapshot file is generated. This values are
# then included in the snapshot file. Note that the computation of the
# values of the potential costs additional time.
#OPT += -DOUTPUTPOTENTIAL

# Update potential every output even if EVALPOTENTIAL is set.
#OPT += -DRECOMPUTE_POTENTIAL_ON_OUTPUT

# This will include the physical acceleration of each particle in
# snapshot files.
#OPT += -DOUTPUTACCELERATION

# This will include the rate of change of entropy of gas particles in
# snapshot files.
#OPT += -DOUTPUTCHANGEOFENTROPY

# This will include an output of the timesteps actually taken by each
# particle.
#OPT += -DOUTPUTTIMESTEP

# Needed when HDF5 I/O support is desired.
OPT += -DHAVE_HDF5
OPT += -DH5_USE_16_API

# Needed when HDF5 v1.6.4 is used. It fixes a memory 'free' bug.
#OPT += -DUSE_HDF5_FIX #<------------

# Enables on-the-fly output of lines-of-sight.
#OPT += -DOUTPUTLINEOFSIGHT
#OPT += -DMyMaxLosBuffer=4096
#OPT += -DOUTPUTLINEOFSIGHT_SPECTRUM
#OPT += -DOUTPUTLINEOFSIGHT_PARTICLES

# Enable adaptive output of particles
#OPT += -DADAPTIVE_OUTPUT

#============================================================================
# Testing and Debugging options
#============================================================================

# This can be set to check the force accuracy of the code. The option
# needs to be set to a number between 0 and 1 (e.g. 0.01), which is
# taken to specify a random fraction of particles for which at each
# timestep forces by direct summation are computed. The normal
# tree-forces and the "correct" direct summation forces are collected
# in a file. Note that the simulation itself is unaffected by this
# option, but it will of course run much(!) slower if
# FORCETEST*NumPart*NumPart >> NumPart. NOTE: ParticleIDs must be set
# to numbers >=1 for this to work.
#OPT += -DFORCETEST=0.1

# Enables core-dumps and FPU exceptions.
#OPT += -DDEBUG

# Auxiliary communication of IDs.
#OPT += -DPARTICLE_DEBUG


#============================================================================
# Static NFW Potential
#============================================================================

#OPT += -DSTATICNFW
#OPT += -DNFW_C=9
#OPT += -DNFW_M200=100.0
#OPT += -DNFW_Eps=0.025
#OPT += -DNFW_DARKFRACTION=0.9


#============================================================================
# Static Hernquist Potential
#============================================================================

#OPT += -DSTATICHQ
#OPT += -DHQ_M200=1.0
#OPT += -DHQ_C=10
#OPT += -DHQ_DARKFRACTION=0.9


#============================================================================
# Dark energy
#============================================================================

# Enables Dark Energy.
#OPT += -DDARKENERGY

# Read w(z) from a file.
#OPT += -DTIMEDEPDE

# Rescale v_ini in read_ic / read_ic_cluster.
#OPT += -DRESCALEVINI


#============================================================================
# SPH algorithm
#============================================================================

OPT += -DPRESSURE_ENTROPY_SPH

# Cubic B-spline kernel (suggested number of neighbours is 58).
#OPT += -DSPH_KERNEL_GADGET

# Cubic B-spline kernel (suggested number of neighbours is 58).
#OPT += -DSPH_KERNEL_CUBIC

# Quartic B-spline kernel (suggested number of neighbours is ?).
#OPT += -DSPH_KERNEL_QUARTIC

# Quintic B-spline kernel (suggested number of neighbours is 195).
#OPT += -DSPH_KERNEL_QUINTIC

# Wendland's C2 kernel (suggested number of neighbours is 100).
#OPT += -DSPH_KERNEL_C2

# Wendland's C4 kernel (suggested number of neighbours is 200).
#OPT += -DSPH_KERNEL_C4

# Wendland's C6 kernel (suggested number of neighbours is 400).
OPT += -DSPH_KERNEL_C6

OPT += -DUSER_GAMMA=1.4


#============================================================================
# SPH viscosity options
#============================================================================

# Enables the old viscosity.
#OPT += -DCONVENTIONAL_VISCOSITY     

# Enables time dependend viscosity.
#OPT += -DTIME_DEP_ART_VISC          

# Turns off the shear viscosity supression.
OPT += -DNO_SHEAR_VISCOSITY_LIMITER 

# Start with high rather than low viscosity.
#OPT += -DHIGH_ART_VISC_START        

# Enables alternative viscosity based on div(v).
#OPT += -DALTVISCOSITY               


#============================================================================
# Glass making
#============================================================================

#OPT += -DMAKEGLASS


#============================================================================
# Baryonic physics - basic options
#============================================================================

# Enables gas cooling (Wiersma et al., 2008).
#OPT += -DBG_COOLING

# Switch this on if you want to use the old (buggy) cooling tables
#OPT += -DBG_COOLING_OLD_TABLES

# Mimic self-shielding of UV background above some density threshold
#OPT += -DBG_COOLING_SHIELDING
#OPT += -DBG_INSTANTANEOUS_REIONIZATION

# Enables star formation (Schaye & Dalla Vecchia, 2008).
#OPT += -DBG_SFR


#============================================================================
# Baryonic physics - kinetic SNII feedback
#============================================================================

# Enables kinetic feedback (Dalla Vecchia & Schaye, 2008).
#OPT += -DBG_SNII_KINETIC_FEEDBACK

# Enables temporary switch off of star formation for wind particles.
#OPT += -DBG_SNII_KINETIC_FEEDBACK_SF_DECOUPLING


#============================================================================
# Baryonic physics - thermal SNII feedback
#============================================================================

# Enables thermal feedback (Dalla Vecchia & Schaye, 2009).
#OPT += -DBG_SNII_THERMAL_FEEDBACK

#OPT += -DBG_THERMAL_FEEDBACK_MAKE_PARTICLES_ACTIVE
#OPT += -DTIMESTEP_LIMITER


#============================================================================
# Baryonic physics - POPIII star formation and evolution
#============================================================================

# Enables POPIII star formation mode
#OPT += -DBG_POPIII

# Enables POPIII thermal feedback
#OPT += -DBG_POPIII_THERMAL_FEEDBACK

# Enables seeding of BHs from POPIII stars
#OPT += -DBG_POPIII_BH_SEEDS


#============================================================================
# Baryonic physics - chemodynamics
#============================================================================

# Enables stellar evolution (metal production, mass transfer, SNIa
# feedback) (Wiersma et al., 2009).
#OPT += -DBG_STELLAR_EVOLUTION

# Uses smoothed metal fractions and metallicities for cooling.
#OPT += -DBG_METALSMOOTHING

# Enables the use of a top-heavy IMF above some density threshold.
#OPT += -DBG_DOUBLE_IMF

# Tracks the Iron produced by SNIa.
#OPT += -DBG_SNIA_IRON

# Enable dust production
#OPT += -DBG_DUST
#OPT += -DBG_DUST_METAL_COOLING_CORRECTION
#OPT += -DBG_DUST_DESTRUCTION_SPUTTERING
#OPT += -DBG_DUST_DESTRUCTION_SUBLIMATION

#OPT += -DLW_BACKGROUND
#OPT += -DLW_LOCAL

# Computes redshifts weighted by metallicity.
#OPT += -DBG_Z_WEIGHTED_REDSHIFT

# Outputs a 3-dimensional (temperature, density, metallicity) grid of
# gas mass.
#OPT += -DBG_OUTPUT_GRID


#============================================================================
# Baryonic physics - molecules
#============================================================================

# Option for molecular network
#OPT += -DBG_MOL_NETWORK

# Use tabulated rate coefficients
#OPT += -DBG_MOL_NETWORK_TABULATED_RATES

# Option for molecular (H2 and HD) cooling
#OPT += -DBG_MOL_COOLING

# Option for molecular rates based time step
#OPT += -DBG_MOL_NETWORK_TIMESTEP
#OPT += -DBG_MOL_NETWORK_H_TIMESTEP
#OPT += -DBG_MOL_NETWORK_H2_TIMESTEP


#============================================================================
# Baryonic physics - extra arrays
#============================================================================

# Outputs maximum temperature and entropy, and when they are achieved.
#OPT += -DBG_EXTRA_ARRAYS


#============================================================================
# Baryonic physics - extra arrays
#============================================================================

# Saves a snapshot instead of a restart file. The snapshot directory
# is called restart_NNN where NNN is the number of the last
# snapshot. At restart the next snapshot will have number NNN+1.
#OPT += -DBG_SNAPSHOT_AS_RESTART_FILE


#============================================================================
# Baryonic physics - extra arrays
#============================================================================

# Output extra information - not sure it still works.
#OPT += -DBG_VERBOSE

# Output extra (debugging) log files - they may be huge!
#OPT += -DBG_VERBOSE_LOGFILES


#============================================================================
# Code testing
#============================================================================

# Tests the cooling routines. If set to 1 it will output interpolated
# cooling rate values between two redshifts (to be updated). If set to
# 2 it will evolve the temperature at the mean density including
# cooling and heating.
#OPT += -DBG_TEST_COOL=2


#============================================================================
# Select here the target architecture
#============================================================================

SYSTYPE="SFC_INTEL"
#SYSTYPE="SFC_GNU"
#SYSTYPE="Hydra"
#SYSTYPE="MacOSX"


#============================================================================
# Specific compilation flags
#============================================================================

CC = mpicc -std=c99
OPTIMIZE = -Wall -g

MPICHLIB = -lmpich
HDF5LIB = -lhdf5

ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(OPT)))
CC = mpiCC
OPTIMIZE = -g
OPT += -DX86FIX
endif

ifeq ($(SYSTYPE),"SFC_INTEL")
CC = mpiicc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC = mpiicpc
endif
OPTIMIZE = -m64 -parallel -O3 -ipo -funroll-loops -no-vec -Wall
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)
FFTW_INCL = -I$(FFTW_HOME)/include
FFTW_LIBS = -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB = -L$(HDF5_HOME)/lib -lhdf5
endif

ifeq ($(SYSTYPE),"SFC_GNU")
CC = mpigcc
ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC = mpigxx
endif
sOPTIMIZE = -m64 -O3 -Wall
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)
FFTW_INCL = -I$(FFTW_HOME)/include
FFTW_LIBS = -L$(FFTW_HOME)/lib
MPICHLIB =
HDF5INCL = -I$(HDF5_HOME)/include
HDF5LIB = -L$(HDF5_HOME)/lib -lhdf5
endif

ifeq ($(SYSTYPE),"Hydra")
CC        = mpiicc
OPTIMIZE  = -m64 -parallel -O3 -ipo -funroll-loops -no-vec -Wall
GSL_INCL  = -I$(GSL_HOME)/include
GSL_LIBS  = -L$(GSL_HOME)/lib
FFTW_INCL = -I$(FFTW_HOME)/include
FFTW_LIBS = -L$(FFTW_HOME)/lib
HDF5INCL  = -I$(HDF5_HOME)/include
HDF5LIB   = -L$(HDF5_HOME)/lib -lhdf5 -lz
MPICHLIB  =
endif

ifeq ($(SYSTYPE),"MacOSX")
CC        = mpicc
OPTIMIZE  = -m64 -O3 -funroll-loops -Wall
GSL_INCL  = -I/opt/local/include
GSL_LIBS  = -L/opt/local/lib
FFTW_INCL =
FFTW_LIBS =
HDF5INCL  =
HDF5LIB   = -lhdf5 -lz
MPICHLIB  =
endif

ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(OPT)))
HDF5INCL =
HDF5LIB =
endif


OPTIONS = $(OPTIMIZE) $(OPT)
 
EXEC = AnarchySPHCosmo

OBJS = accel.o adaptive_outputs.o allocate.o allvars.o begrun.o bg_bhfeedback.o \
       bg_cooling.o bg_cooling_old.o bg_cool_tables.o bg_cool_tables_old.o \
       bg_enrich.o bg_feedback.o bg_init_cool.o bg_init_cool_old.o \
       bg_init.o bg_init_yields.o bg_interpol.o bg_io.o bg_molecules.o \
       bg_popevolv.o bg_sfr.o bg_yield_tables.o binning.o \
       blackhole.o darkenergy.o density.o domain.o driftfac.o endrun.o fof.o \
       forcetree.o global.o gravtree.o hdf5_mymalloc.o hydra.o init.o io.o \
       kernel.o lineofsight.o longrange.o main.o mymalloc.o ngb.o parallel_sort.o \
       peano.o pm_nonperiodic.o pm_periodic.o potential.o predict.o \
       read_ic.o restart.o run.o sigvel.o system.o timeofday.o \
       timestep.o

# added source files: kernel.o
# removed source files: bg_dust.o, lw.o, bg_popiii.o

INCL = Makefile allvars.h bg_cooling.h bg_cooling_old.h bg_dust.h bg_proto.h bg_molecules.h \
       bg_vars.h bg_yields.h dd.h domain.h fof.h forcetree.h functions.h kernel.h proto.h tags.h

# added header files: functions.h kernel.h
# removed header files: 

CFLAGS = $(OPTIONS) $(FFTW_INCL) $(HDF5INCL) $(GSL_INCL)

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

LIBS = -lm $(HDF5LIB) $(FFTW_LIB) $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas


$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
