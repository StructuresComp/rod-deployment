#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper,
	double m_viscosity)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;

	Id3<<1,0,0,
		0,1,0,
        0,0,1;

}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i=0; i < rod->nv; i++)
	{
		u = rod->getVelocity(i);
		f = - viscosity * u * rod->voronoiLen(i);
	}
}


void dampingForce::computeJd()
{
	// Remember that dF/dx = 1/dt * dF/dv

	for (int i = 0; i < rod->nv; i++)
	{
		jac = - viscosity * Id3 * rod->voronoiLen(i)/rod->dt;
		for (int kx = 0; kx < 3; kx++)
		{
			indx = 4 * i + kx;
			for (int ky = 0; ky < 3; ky++)
			{
				indy = 4 * i + ky;
				stepper->addJacobian(indx, indy, - jac(kx,ky)); // subtracting external force
			}
		}
	}
}
