genfitCR
===============

genfit package from Sakamoto san

###1. Documents/Homepages

1. [Genfit](http://genfit.sourceforge.net/Main.html)
2. [NIM on Genfit](http://dx.doi.org/10.1016/j.nima.2010.03.136)
3. [ROOT VMC](http://root.cern.ch/drupal/content/vmc)
4. [GEANE](http://nuint.ps.uci.edu/dcasper/files/GEANE.pdf)

###2. Installation

Genfit (C/C++) uses GEANE (Geant3) via ROOT's Geant3 VMC (C++).
Following is a procedure in the terminal to compile Genfit and test code for Comet Phase 1.

    $ tar xzvf genfit_2012_04_09.tgz                       # Extract tar file
    $ ls 
    README GeaneTrackRep2.cxx GeaneTrackRep2.h GFAbsTrackRep.h media.C geant321+_vmc.1.12.tar.gz genfitV1.1.tar.gz
    $ tar xzvf genfitV1.1.tar.gz                           # Directory `genfit` is created
    $ tar xzvf geant321+_vmc.1.12.tar.gz -C genfit         # Update geant3 directory (ROOT's Geant3 VMC)
    $ cp media.C genfit/geometry/                          # Patch 1. Add meterial
    $ cp GFAbsTrackRep.h genfit/core/                      # Patch 2. Define a custom function to retrieve a volume name during tracking
    $ cp GeaneTrackRep2.h genfit/GeaneTrackRep2/           # Cont'd
    $ cp GeaneTrackRep2.cxx genfit/GeaneTrackRep2/         # Cont'd
    $ root-config --f77                                    # Check if your ROOT knows fortran compiler (If nothing shown, you need to install gfortran then recompile ROOT)
    gfortran
    $ cd genfit/geant3                                     # Move to geant3 direcotry
    $ make                                                 # Take a few minuites
    $ ls lib/tgt_linuxx8664gcc/                            # On Ubuntu (64bit)
    libgeant321.so
    $ ln -s lib/tgt_linuxx8664gcc lib/tgt_linux            # Make symbolic link 
    $ cd ../                                               # Move to top direcotry of Genfit
    $ export GENFIT=`pwd`                                  # Setup enrivonment
    $ export VMC=$GENFIT/geant3                            # Cont'd
    $ ./makeEnv.sh                                         # Make env.sh, needed at every time you compile/run GENFIT
    $ source env.sh                                        # Read environment for Genfit
    $ cmake .                                              # Generate Makfile for Genfit (If cmake is not install, do sudo apt-get install cmake)
    $ make
    $ ls lib/
    libgenfitGeane.so  libgenfitLSL.so  libgenfitRK.so  libgenfitRKXY.so  libgenfitSlTrackRep.so  libgenfit.so
    $ tar xzvf ../CometPhase1_2012_04_09.tgz -C test       # Extract test code for Comet phase 1
    $ cd test/CometPhase1
    $ cmake .                                              # Generate Makefile
    $ make
    $ ./main -h                                            # Options and example case is displayed
    $ ./main <options>                                     # Run program (~ 20k events/hour with Intel(R) Core(TM) i7 CPU 940  @ 2.93GHz)

###3. Issues/Limits

1.  Deposit energy in materials are not recorded in ROOT.
2.  With option `-j RKTrackRep|GeaneTrackRep2`, hit positions are calculated by following a track in a step size as

        R <= target's radius: step_size = half_of_target_thickness
        R > target's radius:  step_size = 0.1 cm
        
    ==> Call Geant3 user routine (GUKINE etc.) instead of RKTrackRep or GeaneTrackRep2 for making true hits.
    More efficient way than current. (and also energy deposit can be recorded)  

3.  In tracking, hit positions after second turn are not handled by Genfit since the limit of TMAXFD = 20 degree  
    geant3/gphys/gphysi.F:725  

        IF (TMAXFD.LE.0..OR. (IGAUTO.NE.0.AND.TMAXFD.GT.20.)) THEN  
            TMAXFD=20.  

    ==> No idea at the moment
