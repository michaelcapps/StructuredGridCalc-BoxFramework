#ifndef _LBPATCH_H_
#define _LBPATCH_H_
#include "LBPhysics.H"
#include "LevelData.H"
#include "BaseFab.H"

namespace LBPatch
{
using SolFab = BaseFab<Real>;
void macroscopic(DisjointBoxLayout& a_dbl,LevelData<SolFab>& curr,LevelData<SolFab>& macro)
{
	IntVect temp;
	for(DataIterator dit(a_dbl);dit.ok();++dit)
	{
		MD_BOXLOOP_OMP(a_dbl[dit],i)
		{
			temp = IntVect(i0,i1,i2);
			LBPhysics::macroscopic(macro[dit],curr[dit],temp);
			//std::cout << macro[dit](temp,0) << std::endl;
		}	
	}

}

//Collision function
void collision(LevelData<SolFab> &curr, LevelData<SolFab>& macro,DisjointBoxLayout &a_dbl)
{
	IntVect temp;
	Real u[3];
	Real rho;
	for(DataIterator dit(a_dbl);dit.ok();++dit)
	{
		MD_BOXLOOP_OMP(a_dbl[dit],i)
		{
			temp = IntVect(i0,i1,i2);
			rho = (macro[dit])(temp,0);
			u[0] = (macro[dit])(temp,1);
			u[1] = (macro[dit])(temp,2);
			u[2] = (macro[dit])(temp,3);		
			for(int k = 0;k<LBParameters::g_numVelDir;++k)	
			{
				LBPhysics::collision(k,curr[dit](temp,k),u,rho);
			}
		}	
	}	
}


void stream(DisjointBoxLayout& a_dbl, LevelData<SolFab>& m_curr, LevelData<SolFab>& m_prev)
{
	Box src_box;
	Box dst_box;
	IntVect shift_dir;
	for(DataIterator dit(a_dbl);dit.ok();++dit)
	{
		for(int k = 0;k<LBParameters::g_numVelDir;++k)
		{
			shift_dir=LBParameters::latticeVelocity(LBParameters::oppositeVelDir(k));		
			dst_box = a_dbl[dit];
			src_box = dst_box;
			src_box = src_box.shift(shift_dir);
			m_prev[dit].copy(dst_box,k,m_curr[dit],src_box,k,1); 
		}
	}

	//Move m_prev to m_curr and m_curr to m_prev
	LevelData<SolFab> temp;
	temp = std::move(m_curr);
	m_curr = std::move(m_prev);
	m_prev = std::move(temp);
}//end stream

}//end namespace LBPatch


#endif //header guard
