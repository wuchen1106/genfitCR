genfitCR
===============

genfit tookit for Comet Reconstruction

###1. Documents/Homepages

1. [Genfit](http://genfit.sourceforge.net/Main.html)
2. [NIM on Genfit](http://dx.doi.org/10.1016/j.nima.2010.03.136)
3. [ROOT VMC](http://root.cern.ch/drupal/content/vmc)
4. [GEANE](http://nuint.ps.uci.edu/dcasper/files/GEANE.pdf)

###2. Installation

Genfit (C/C++) uses GEANE (Geant3) via ROOT's Geant3 VMC (C++).
Following is a procedure in the terminal to compile Genfit and test code for Comet Phase 1.  

    $ root-config --f77                         # Check if your ROOT knows fortran compiler 
                                                  (If nothing shown, you need to install gfortran then recompile ROOT)
                                                  ( for ArchLinux: sudo pacman -Sy gcc-fortran)
      gfortran
    $ cd genfit/geant3                          # Move to geant3 direcotry
    $ make                                      # Take a few minuites
    $ ls lib/tgt_linuxx8664gcc/                 # On 64bit Linux
      libgeant321.so
    $ ln -s lib/tgt_linuxx8664gcc lib/tgt_linux # Make symbolic link 
    $ cd ../                                    # Move to top direcotry of Genfit
    $ export GENFIT=`pwd`                       # Setup enrivonment
    $ export VMC=$GENFIT/geant3                 # Cont'd
    $ ./makeEnv.sh                              # Make env.sh, needed at every time you compile/run GENFIT
    $ source env.sh                             # Read environment for Genfit
    $ cmake .                                   # Generate Makfile for Genfit
                                                  ( you should install cmake first )
                                                  ( for Ubuntu: sudo apt-get install cmake )
                                                  ( for ArchLinux: sudo pacman -Sy cmake)
    $ make
    $ ls lib/
      libgenfitGeane.so  libgenfitLSL.so  libgenfitRK.so  libgenfitRKXY.so  libgenfitSlTrackRep.so  libgenfit.so
    $ cd ../CometPhase1
    $ cmake .                                   # Generate Makefile
    $ make
    $ ./main -h                                 # Options and example case is displayed
    $ ./main <options>                          # Run program (~ 20k events/hour with Intel(R) Core(TM) i7 CPU 940  @ 2.93GHz)

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
