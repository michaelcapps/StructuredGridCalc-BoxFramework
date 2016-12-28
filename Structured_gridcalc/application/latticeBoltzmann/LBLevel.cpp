#include "LBLevel.H"
#include "LBPatch.H"
#include "LevelData.H"

//Advance a time step
void LBLevel::advance()
{
	//Collision
	LBPatch::collision(m_curr,m_macro_comps,m_dbl);
	
	//Exchange
	Copier copier;
	copier.defineExchangeLD(m_curr,PeriodicX | PeriodicY,TrimCorner);	
	m_curr.exchange(copier);

	//Fill ghost cells on top/bottom boundary and fill x,y ghost cells using periodic conditions
	//Using non-slip condictions
	Box temp_box;
	IntVect temp_lo;
	IntVect temp_hi;
	IntVect curr_vect;
	for(DataIterator dit(m_dbl);dit.ok();++dit)
	{
		/* NOTE: This stuff was to do periodic conditions "manually", but it would definitely fail if not on 1 or 2 processors
		if(m_dbl[dit].loVect(0)==(m_dbl.problemDomain()).loVect(0))
		{//on left y boundary
			temp_lo = m_dbl[dit].loVect();
			temp_hi = m_dbl[dit].hiVect();
			temp_hi[0] = m_dbl[dit].loVect(0);
			temp_box = Box(temp_lo,temp_hi);
			//std::cout << temp_box << std::endl;
			temp_box = temp_box.shift(IntVect(-1,0,0));//shift so dst box are ghost cells
			DataIterator temp_dit = dit;
			for(int ii = 0;ii<3;++ii){++temp_dit;}
			IntVect period_hi = m_dbl[temp_dit].hiVect();//Should fix the +3 to be a stride
			IntVect period_lo = m_dbl[temp_dit].loVect();
			period_lo[0] = period_hi[0];
			Box periodic_box = Box(period_lo,period_hi);//right boundary
			//std::cout << periodic_box << std::endl;
			m_curr[dit].copy(temp_box,0,m_curr[temp_dit],periodic_box,0,LBParameters::g_numVelDir);
			temp_box = temp_box.shift(IntVect(1,0,0));
			periodic_box = periodic_box.shift(IntVect(1,0,0));
			m_curr[temp_dit].copy(periodic_box,0,m_curr[dit],temp_box,0,LBParameters::g_numVelDir);
		}		
		if(m_dbl[dit].loVect(1)==(m_dbl.problemDomain()).loVect(1))
		{//on low x boundary
			temp_lo = m_dbl[dit].loVect();
			temp_hi = m_dbl[dit].hiVect();
			temp_hi[1] = temp_lo[1];
			temp_box = Box(temp_lo,temp_hi);
			temp_box = temp_box.shift(IntVect(0,-1,0));
			//std::cout << temp_box << std::endl;
			DataIterator temp_dit = dit;
			for(int ii = 0;ii<4;++ii){++temp_dit;}
			IntVect period_hi = m_dbl[temp_dit].hiVect();
			IntVect period_lo = m_dbl[temp_dit].loVect();
			period_lo[1] = period_hi[1];
			Box periodic_box = Box(period_lo,period_hi);
			//std::cout << periodic_box << std::endl;
			m_curr[dit].copy(temp_box,0,m_curr[temp_dit],periodic_box,0,LBParameters::g_numVelDir);
			temp_box = temp_box.shift(IntVect(0,1,0));
			periodic_box = periodic_box.shift(IntVect(0,1,0));
   			m_curr[temp_dit].copy(periodic_box,0,m_curr[dit],temp_box,0,LBParameters::g_numVelDir);
		}		
		*/
		if(m_dbl[dit].hiVect(2) == (m_dbl.problemDomain()).hiVect(2))
                {//on top of domain     
                        temp_lo = IntVect(m_dbl[dit].loVect(0),m_dbl[dit].loVect(1),m_dbl[dit].hiVect(2));
                        temp_hi = m_dbl[dit].hiVect();
                        temp_box = Box(temp_lo,temp_hi);
                        MD_BOXLOOP_OMP(temp_box,i)
                        {
                                curr_vect = IntVect(i0,i1,i2);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(6),5)   = m_curr[dit](curr_vect,6);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(13),12) = m_curr[dit](curr_vect,13);
				m_curr[dit](curr_vect+LBParameters::latticeVelocity(14),11) = m_curr[dit](curr_vect,14);
				m_curr[dit](curr_vect+LBParameters::latticeVelocity(17),16) = m_curr[dit](curr_vect,17);
				m_curr[dit](curr_vect+LBParameters::latticeVelocity(18),15) = m_curr[dit](curr_vect,18);
				
				/*
                                m_curr[dit](curr_vect+e6,5)   = m_curr[dit](curr_vect,6);
                                m_curr[dit](curr_vect+e13,12) = m_curr[dit](curr_vect,13);
                                m_curr[dit](curr_vect+e14,11) = m_curr[dit](curr_vect,14);
                                m_curr[dit](curr_vect+e17,16) = m_curr[dit](curr_vect,17);
                                m_curr[dit](curr_vect+e18,15) = m_curr[dit](curr_vect,18);
				*/
                        }
                }
                else if(m_dbl[dit].loVect(2) == (m_dbl.problemDomain()).loVect(2))
                {//on bottom of domain
                        temp_lo = m_dbl[dit].loVect();
                        temp_hi = m_dbl[dit].hiVect();
                        temp_hi[2] = m_dbl[dit].loVect(2);
                        temp_box = Box(temp_lo,temp_hi);
                        MD_BOXLOOP_OMP(temp_box,i)
                        {
                                curr_vect = IntVect(i0,i1,i2);
                         
				m_curr[dit](curr_vect+LBParameters::latticeVelocity(5),6)   = m_curr[dit](curr_vect,5);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(12),13) = m_curr[dit](curr_vect,12);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(11),14) = m_curr[dit](curr_vect,11);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(16),17) = m_curr[dit](curr_vect,16);
                                m_curr[dit](curr_vect+LBParameters::latticeVelocity(15),18) = m_curr[dit](curr_vect,15);
				/*
			       m_curr[dit](curr_vect+e5,6)   = m_curr[dit](curr_vect,5);
                                m_curr[dit](curr_vect+e12,13) = m_curr[dit](curr_vect,12);
                                m_curr[dit](curr_vect+e11,14) = m_curr[dit](curr_vect,11);
                                m_curr[dit](curr_vect+e16,17) = m_curr[dit](curr_vect,16);
                                m_curr[dit](curr_vect+e15,18) = m_curr[dit](curr_vect,15);
				*/
                        }
                }
	}	
	
	
	//Stream
	LBPatch::stream(m_dbl,m_curr,m_prev);
	
	//Macroscopic
	LBPatch::macroscopic(m_dbl,m_curr,m_macro_comps);
}

/*--------------------------------------------------------------------*/
//  Compute mass in domain
/** Just a sum of fi()
 *  \return             Total mass in domain in process 0
 *//*-----------------------------------------------------------------*/

Real LBLevel::computeTotalMass() const
{
  Real localDomainMass = 0.;
  for (DataIterator dit(m_dbl); dit.ok(); ++dit)  //**FIX m_boxes
    {
      const Box& box = m_dbl[dit];
      MD_ARRAY_RESTRICT(arrfi, m_curr[dit]);  //**FIX method to access current fi
      for (int iVel = 0; iVel != LBParameters::g_numVelDir; ++iVel)
        {
          MD_BOXLOOP_OMP(box, i)
            {
              localDomainMass += arrfi[MD_IX(i, iVel)];
            }
        }
    }
  Real globalDomainMass = localDomainMass;
#ifdef USE_MPI
  // Accumulated the result into process 0
  MPI_Reduce(&localDomainMass, &globalDomainMass, 1, BXFR_MPI_REAL, MPI_SUM, 0,
             MPI_COMM_WORLD);
#endif
  return globalDomainMass;

}
