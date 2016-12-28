#ifndef _LBLEVEL_H
#define _LBLEVEL_H

#include "LBParameters.H"
#include "BaseFabMacros.H"
#include "LevelData.H"
#ifdef USE_MPI
#include "pcgnslib.h"
#else
#include "cgnslib.h"
#endif

class LBLevel
{
  using LevelSolData = LevelData<BaseFab<Real> >;

public:
	//Constructors
	LBLevel();//default
	LBLevel(DisjointBoxLayout &a_dbl); //construct with dbl
	LBLevel(const DisjointBoxLayout &a_dbl); //const construct with dbl
	//Destructors
 	//~LBLevel();//default

public: //member functions
	void initialData();
	void advance();
  	int writePlotFile(int iter) const;
  	Real computeTotalMass() const;

protected: //data members
	DisjointBoxLayout m_dbl;
	LevelSolData m_curr;
	LevelSolData m_prev;
	LevelSolData m_macro_comps;

	//Store some IntVects to make filling ghost cells simple
	IntVect e6  = IntVect(0,0,1); // +z
	IntVect e14 = IntVect(1,0,1); // +x,+z
	IntVect e13 = IntVect(-1,0,1); // -x,+z
	IntVect e18 = IntVect(0,1,1); // +y,+z
	IntVect e17 = IntVect(0,-1,1); // -y, +z
	IntVect e5  = IntVect(0,0,-1); // -z
	IntVect e12 = IntVect(1,0,-1); // +x,-z
	IntVect e11 = IntVect(-1,0,-1); // -x,-z
	IntVect e16 = IntVect(0,1,-1); // +y,-z
	IntVect e15 = IntVect(0,-1,-1); // -y,-z
};


/*********CONSTRUCTORS**********/
//Default
inline LBLevel::LBLevel()
:
m_dbl(),
m_curr(),
m_prev(),
m_macro_comps()
{}

//Construction with dbl
inline LBLevel::LBLevel(DisjointBoxLayout &a_dbl)
:
m_dbl(a_dbl),
m_curr(a_dbl,LBParameters::g_numVelDir,LBParameters::g_numGhost),
m_prev(a_dbl,LBParameters::g_numVelDir,LBParameters::g_numGhost),
m_macro_comps(a_dbl,4,LBParameters::g_numGhost)
{
	initialData();
}

//Construction with const dbl
inline LBLevel::LBLevel(const DisjointBoxLayout &a_dbl)
:m_dbl(a_dbl),
m_curr(a_dbl,LBParameters::g_numVelDir,LBParameters::g_numGhost),
m_prev(a_dbl,LBParameters::g_numVelDir,LBParameters::g_numGhost),
m_macro_comps(a_dbl,4,LBParameters::g_numGhost)
{
	initialData();
}

/********MEMBER FUNCTIONS*********/
inline void LBLevel::initialData()
{
	for(int k = 0; k<LBParameters::g_numVelDir;++k)
	{
		m_curr.setVal(k,LBParameters::g_weight[k]);
		m_prev.setVal(k,LBParameters::g_weight[k]);
	}
	m_macro_comps.setVal(0,1.0);
	for(int k = 1;k<4;++k)
	{
		m_macro_comps.setVal(k,0.0);
	}
}

inline int LBLevel::writePlotFile(int iter) const
{
	const bool verbose = 0;
	IntVect origin = IntVect::Zero;
	Real dx = 1;
	int cgerr;
	//Open the CGNS file
	int indexFile;
  	char  fileName[33]={0};
	strcat(fileName,"plot/solution");
	char tmp_str[33];
	sprintf(tmp_str,"%04d",iter);
	strcat(fileName,tmp_str);
	strcat(fileName,".cgns");
	
	#ifdef USE_MPI
  		cgerr = cgp_open(fileName, CG_MODE_WRITE, &indexFile);
	#else
  		cgerr = cg_open(fileName, CG_MODE_WRITE, &indexFile);
	#endif
  	if (cgerr)
    	{
      		if (verbose) cg_error_print();
      		return cgerr;
    	}
	
	//Create the base node
	int indexBase;
  	int iCellDim = g_SpaceDim;
  	int iPhysDim = g_SpaceDim;
  	cgerr = cg_base_write(indexFile, "Base", iCellDim, iPhysDim, &indexBase);
  	if (cgerr)
    	{
      	if (verbose) cg_error_print();
      		return cgerr;
    	}	

	// Write the zones and coordinates (indexZoneOffset needs to be determined
	//   // so we can use 'indexCGNSzone = globalBoxIndex + indexZoneOffset'.  This
	//     // value is almost certainly 1).	
	
	if (verbose)
    	{
      		std::cout << "Writing Grid Coordinates\n";
    	}
	int indexZoneOffset;  	// The difference between CGNS indexZone and the
                        	// globalBoxIndex.
	cgerr = m_dbl.writeCGNSZoneGrid(indexFile,
                                  indexBase,
                                  indexZoneOffset,
                                  origin,
                                  dx);
 	cg_error_print();
  	if (cgerr)
    	{
      		if (verbose) cg_error_print();
      		return cgerr;
    	}
	if (verbose)
    	{
      		std::cout << "Writing Solution Data\n";
    	}
  	cgerr = m_macro_comps.writeCGNSSolData(indexFile,
                             indexBase,
                             indexZoneOffset,
                             LBParameters::stateNames());
  	if (cgerr)
    	{
      	if (verbose) cg_error_print();
      		return cgerr;
    	}
		
	//Close the CGNS file
	#ifdef USE_MPI
  		cgerr = cgp_close(indexFile);
	#else
  		cgerr = cg_close(indexFile);
	#endif
  	if (cgerr)
    	{
      	if (verbose) cg_error_print();
      		return cgerr;
    	}

  	if (verbose)
    	{
      		std::cout << "Successfully wrote " << fileName << std::endl;
    	}
	return 0;
		
}//end writePlotFile

#endif  //header guard
