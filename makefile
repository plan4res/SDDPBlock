##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SDDPBlock                                                    #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by SDDPBlock. #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the libSMS++ library itself                     #
#           $(StcBlkH)     = the .h files to include for StochasticBlock     #
#           $(StcBlkINC)   = the -I$( StochasticBlock source directory )     #
#           $(libStOptINC) = the -I$(include directories) for libStOpt       #
#           $(SDDPBkSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(SDDPBkOBJ)   = the final object(s) / library                   #
#           $(SDDPBkH)     = the .h files to include                         #
#           $(SDDPBkINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

SDDPBkOBJ = $(SDDPBkSDR)obj/SDDPBlock.o $(SDDPBkSDR)obj/SDDPSolver.o \
	$(SDDPBkSDR)obj/SDDPGreedySolver.o \
	$(SDDPBkSDR)obj/ParallelSDDPSolver.o

SDDPBkINC = -I$(SDDPBkSDR)include/

SDDPBkH   = $(SDDPBkSDR)include/SDDPBlock.h \
	$(SDDPBkSDR)include/ScenarioSet.h \
	$(SDDPBkSDR)include/ScenarioSimulator.h \
	$(SDDPBkSDR)include/SDDPSolver.h \
	$(SDDPBkSDR)include/SDDPGreedySolver.h \
	$(SDDPBkSDR)include/ParallelSDDPSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SDDPBkOBJ) $(MILPBSDR)*~

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(SDDPBkSDR)obj/SDDPBlock.o: $(SDDPBkSDR)src/SDDPBlock.cpp \
	$(SDDPBkSDR)include/SDDPBlock.h \
	$(SDDPBkSDR)include/ScenarioSet.h \
	$(SDDPBkSDR)include/ScenarioSimulator.h $(StcBlkH) $(SMS++OBJ)
	$(CC) -c $(SDDPBkSDR)src/SDDPBlock.cpp -o $@ $(SDDPBkINC) \
	$(StcBlkINC) $(SMS++INC) $(SW)

$(SDDPBkSDR)obj/SDDPSolver.o: $(SDDPBkSDR)src/SDDPSolver.cpp \
	$(SDDPBkSDR)include/SDDPSolver.h $(SDDPBkSDR)include/SDDPBlock.h \
	$(SDDPBkSDR)include/ScenarioSet.h \
	$(SDDPBkSDR)include/ScenarioSimulator.h $(StcBlkH) $(SMS++OBJ)
	$(CC) -c $(SDDPBkSDR)src/SDDPSolver.cpp -o $@ $(SDDPBkINC) \
	$(StcBlkINC) $(libStOptINC) $(SMS++INC) $(SW)

$(SDDPBkSDR)obj/SDDPGreedySolver.o: $(SDDPBkSDR)src/SDDPGreedySolver.cpp \
	$(SDDPBkSDR)include/SDDPGreedySolver.h \
	$(SDDPBkSDR)include/SDDPBlock.h \
	$(SDDPBkSDR)include/ScenarioSet.h \
	$(SDDPBkSDR)include/ScenarioSimulator.h $(StcBlkH) $(SMS++OBJ)
	$(CC) -c $(SDDPBkSDR)src/SDDPGreedySolver.cpp -o $@ $(SDDPBkINC) \
	$(StcBlkINC) $(SMS++INC) $(SW)

$(SDDPBkSDR)obj/ParallelSDDPSolver.o: $(SDDPBkSDR)src/ParallelSDDPSolver.cpp \
	$(SDDPBkSDR)include/ParallelSDDPSolver.h \
	$(SDDPBkSDR)include/SDDPSolver.h $(SDDPBkSDR)include/SDDPBlock.h \
	$(SDDPBkSDR)include/ScenarioSet.h \
	$(SDDPBkSDR)include/ScenarioSimulator.h $(StcBlkH) $(SMS++OBJ)
	$(CC) -c $(SDDPBkSDR)src/ParallelSDDPSolver.cpp -o $@ $(SDDPBkINC) \
	$(StcBlkINC) $(libStOptINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
