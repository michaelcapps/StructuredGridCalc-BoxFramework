#include "LBLevel.H"
#include "Stopwatch.H"
#include <chrono>

/******************************************************************************/
/**
 * \file latticeBoltzmann.cpp
 *
 * \brief Lattice-Boltzmann solution to a problem
 *
 *//*+*************************************************************************/

int main(int argc, const char* argv[])
{
	//Test code
	Box domain(IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(63, 31, 31)));
 	 // 2 boxes in x direction

#ifdef USE_MPI
 	//std::ostringstream ost;
	const bool verbose = 0;

	DisjointBoxLayout::initMPI(argc, argv);
  	int numProc = DisjointBoxLayout::numProc();
  	int procID = DisjointBoxLayout::procID();
	const bool masterProc = (procID == 0);

  	if (masterProc)
    	{
      		if (verbose) std::cout << "Using " << numProc << " processors\n";
    	}
#endif

	Stopwatch<std::chrono::steady_clock> stopwatch;
	stopwatch.start();
  	DisjointBoxLayout dbl(domain, 16*IntVect::Unit);
	LBLevel lblvl(dbl); //constructor with dbl

	for(int k = 0; k<4001; ++k)
	{
		//iterate
		if(k%200==0)
		{
			std::cout << "Writing during iteration "<< k << std::endl;
			lblvl.writePlotFile(k);
		}
		
		//std::cout << lblvl.computeTotalMass() << std::endl;
		lblvl.advance();
	}
	
	stopwatch.stop();
	std::cout <<stopwatch.time() << std::endl;
	DisjointBoxLayout::finalizeMPI();				
	
  // Setup input parameters
  // Setup data structures
  // Write the initial plot file
  // Iterate
  // Write the final plot file
  return 0;
}
