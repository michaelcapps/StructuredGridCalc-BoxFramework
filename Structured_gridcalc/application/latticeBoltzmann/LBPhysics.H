#ifndef _LBPHYSICS_H_
#define _LBPHYSICS_H_
namespace LBPhysics
{
void collision(int i, Real &fi, Real* u, Real rho)
{
	Real ei_dot_u = u[0]*LBParameters::latticeVelocityP(i)[0]+u[1]*LBParameters::latticeVelocityP(i)[1]+u[2]*LBParameters::latticeVelocityP(i)[2];

	Real fi_eq = LBParameters::g_weight[i]*rho*(1+ei_dot_u/LBParameters::g_cs2 + 
		ei_dot_u*ei_dot_u/(2*LBParameters::g_cs2*LBParameters::g_cs2)-
		(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])/(2*LBParameters::g_cs2));

	fi = fi + (fi_eq - fi)/LBParameters::g_tau + 3*LBParameters::g_weight[i]*LBParameters::latticeVelocityP(i)[0]*LBParameters::g_bodyForce;
}//end collision


void macroscopic(BaseFab<Real>& macro, BaseFab<Real>& curr,IntVect& a_cell)
{
	for(int k = 0; k<4; ++k){macro(a_cell,k)=0;}//clear rho and u for new computations

	for(int k = 0; k<LBParameters::g_numVelDir; ++k)
	{
		macro(a_cell,0) += curr(a_cell,k);
		macro(a_cell,1) += curr(a_cell,k)*LBParameters::latticeVelocityP(k)[0];
		macro(a_cell,2) += curr(a_cell,k)*LBParameters::latticeVelocityP(k)[1];
		macro(a_cell,3) += curr(a_cell,k)*LBParameters::latticeVelocityP(k)[2];
	}

	//Divide by rho
	macro(a_cell,1) = macro(a_cell,1)/macro(a_cell,0);
	macro(a_cell,2) = macro(a_cell,2)/macro(a_cell,0);
	macro(a_cell,3) = macro(a_cell,3)/macro(a_cell,0);
}//end macroscopic

}//end LBPhysics

#endif  //header guard
