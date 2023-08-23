#include "world.h"
#include <sstream>
#include <iomanip>
#include <limits>

bool isNegativeNaN(double x) {
    // Create a bit pattern of -nan
    double negativeNaN = -std::numeric_limits<double>::quiet_NaN();

    // Type punning to access the underlying bits
    unsigned long long xBits = *reinterpret_cast<unsigned long long*>(&x);
    unsigned long long negativeNaNBits = *reinterpret_cast<unsigned long long*>(&negativeNaN);

    // Check if the bit patterns are equal
    return xBits == negativeNaNBits;
}


world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean
	gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
	maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	dl_b = m_inputData.GetScalarOpt("dl_b"); // rope mesh size

	hangLength = m_inputData.GetScalarOpt("hangLength"); // length of the suspended part
	kappaT = m_inputData.GetScalarOpt("kappaT"); // target kappa
	kappaI = m_inputData.GetScalarOpt("kappaI"); // initial guess kappa
	Poisson = m_inputData.GetScalarOpt("Poisson");       // Poisson's ratio
	normKs = m_inputData.GetScalarOpt("normKs"); // normalized stretching stiffness

	xTop = m_inputData.GetScalarOpt("xTop"); // robot grasp x
	yTop = m_inputData.GetScalarOpt("yTop"); // robot grasp y
	zTop = m_inputData.GetScalarOpt("zTop"); // robot grasp z
	// scaling all quantities
	youngM = 1.8e5;
	rodRadius = 1.6e-3;
	density = 1180;
	Lgb = pow(rodRadius * rodRadius * youngM/(8 * density * abs(gVector(2))), 1.0/3.0);

	// compute density
	gVector = gVector * density/1.180;
	density = 1.180;
	cout << youngM/density << endl;


  Lgb = pow(rodRadius * rodRadius * youngM/(8 * density * abs(gVector(2))), 1.0/3.0);
	rodRadius = sqrt(Lgb * Lgb * 4/normKs);
	youngM = pow(Lgb, 3) * (8 * density * abs(gVector(2)))/pow(rodRadius, 2);

	density = youngM * density/1.8e5;
	gVector(2) = -10 * 1180/density;

	Lgb = pow(rodRadius * rodRadius * youngM/(8 * density * abs(gVector(2))), 1.0/3.0);

	// cout << rodRadius << " "<<youngM << endl;
	// cout << youngM/density <<" "<<Lgb<<	 endl;
	// 	exit(0);
  // for numeric stability
  finiteDiff = 1e-9;
  if (hangLength < 3)
  {
    finiteDiff = 1e-7;
  }
  else if (hangLength < 7)
  {
    finiteDiff = 1e-8;
  }


	hangLength = hangLength * Lgb;

	kappaI = kappaI/Lgb;
	kappaT = kappaT/Lgb;

	deltaKappa = 0.01/Lgb;
	xTop = xTop * Lgb;
	yTop = yTop * Lgb;
	zTop = zTop * Lgb;
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus

	savetemp = false;
	innerIter = 0;
  addingInertial = false;

  addStep = 0;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile, string filename)
{
	if (saveData==false) return;

	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;

	name << std::setprecision(4);
	// name << "datafiles/"<< filename;
	// name << "_hangLength_"<< hangLength/Lgb;
	// name << "_ks_"<<normKs;
	// name << "_kappaT_"<<kappaT * Lgb;
	// name << "_Poisson_"<<Poisson;
	// name << "_dl_"<<dl_b;
	// name <<".txt";
	name << "datafiles/simDER.txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false)
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false)
		return;

  Matrix3d R0, Rt;
	R0.col(0) = rod->m1.row(0);
	R0.col(1) = rod->m2.row(0);
	R0.col(2) = rod->tangent.row(0);
	Rt.col(0) = m1Local;
	Rt.col(2) = tLocal;
	Rt.col(1) = tLocal.cross(m1Local);

	Rt = Rt * R0.transpose();



	outfile << xTop/Lgb <<" "<<yTop/Lgb <<" "<<zTop/Lgb <<" "<<F0(0) <<" "<<F0(1)<<" "<<F0(2)<<" "<< yForce/rod->EI*dl_b * Lgb <<
	" "<<Rt(0, 0)<<" "<< Rt(0, 1) << " "<<Rt(0, 2)<<" "<<Rt(1, 0) <<" "<<Rt(1, 1) <<" "<< Rt(1, 2) <<" "<< Rt(2, 0) <<" "<<Rt(2, 1) <<" "<<Rt(2, 2)<<endl;

	if (abs(yForce/rod->EI*dl_b - kappaT)/(kappaT+1e-4) < 0.01)
	{
		exit(0);
	}
	yForce = yForce + deltaKappa * rod->EI/dl_b;
	yForce = yForce < kappaT * rod->EI/dl_b ? yForce : kappaT * rod->EI/dl_b;
}


void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();

	// Create the rod
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, theta);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;

	// cout << forceTol << endl;
	// exit(0);

	// Set up boundary condition
	rodBoundaryCondition();
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup

	yForce = kappaI * rod->EI/dl_b; // set up y force

	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();

	saveData = false;
	searchStep = 1;
	gStep = 0;

	timeStep = 0;
	currentTime = 0.0;
}

void world::rodGeometry()
{

	//read data for rodGeometry
	numVertices = hangLength/dl_b;
	dl_b = hangLength/numVertices;
	RodLength = hangLength + dl_b;
	numVertices = RodLength/dl_b  + 1;

	// cout << numVertices <<" "<<RodLength << endl;
	// exit(0);
	vertices = MatrixXd(numVertices, 3);


	for (int i = 0; i<numVertices; i++)
	{
		vertices(i, 0) = -dl_b+ i*dl_b;
		vertices(i, 1) = 0;
		vertices(i, 2) = 0;
	}

	theta = VectorXd::Zero(numVertices-1);
}



void world::rodBoundaryCondition()
{

	rod->setVertexBoundaryCondition(rod->getVertex(0),0);
  rod->setThetaBoundaryCondition(rod->getTheta(0),0);
  rod->setVertexBoundaryCondition(rod->getVertex(1),1);

  rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1),numVertices-1);
}


void world::updateBoundary()
{
	// Steps
	// 1.move to startpoint
  // 2.move along vertical direction

	saveData = false;
	getGoal = false;
	switch(searchStep){
		case 1:{
			move2NextP(); // xtop, ytop, ztop
			break;
		}
		case 2:{
			sampleFinteDifferenceX();
			break;
		}
		case 3:{
			sampleFinteDifferenceY();
			break;
		}
		case 4:{
			sampleFinteDifferenceZ();
			break;
		}
		case 5:{
			if (run) getGoal = move2goal(Vector3d(xTop - dx(0) * rate, yTop - dx(1) * rate, zTop - dx(2) * rate));
			break;
		}
	}
}


void world::move2NextP()
{
	getGoal = move2goal(Vector3d(xTop, yTop, zTop));

  if (addStep == 0 && getGoal)
  {
    // check all nodes above the ground
    bool underground = false;
    for (int i = 2; i < rod->nv; i++)
    {
      Vector3d xLocal = rod->getVertex(i);
      if (xLocal(2) < 0)
      {
        underground = true;
        break;
      }
    }
    if (underground)
    {
      getGoal = false;
      zTop += 0.01 * hangLength;
    }
  }

  if (getGoal)
  {
    addStep++;
  }
}

void world::sampleFinteDifferenceX()
{
	getGoal = move2goal(Vector3d(xTop+finiteDiff, yTop, zTop));
}

void world::sampleFinteDifferenceY()
{
	getGoal = move2goal(Vector3d(xTop, yTop + finiteDiff, zTop));
}

void world::sampleFinteDifferenceZ()
{
	getGoal = move2goal(Vector3d(xTop, yTop, zTop+finiteDiff));
}

bool world::move2goal(Vector3d goal)
{
	Vector3d temp = rod->getVertex(numVertices-1);
  temp = goal - temp;
	innerIter++;
	if (temp.norm() <= 0.01 * deltaTime)
	{
		rod->setVertexBoundaryCondition(goal , numVertices-1);
		if (rod->u.norm()<finiteDiff)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		temp = temp/temp.norm();
		rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1)+0.01 * temp*deltaTime
			 , numVertices-1);
	}
	return false;
}

void world::updateTimeStep()
{
	// cout <<"begin exploring the optimal"<<endl;

  addingInertial = false;
	bool solved = false;
	while (!solved)
	{
		if (currentTime > 1)
		{
			updateBoundary();
		}
		rod->updateGuess();
		newtonMethod(solved); //get configuration when no collisoins
		if (!solved)
		{
			cout <<"shrink Time step" <<endl;
			deltaTime = 0.1 * deltaTime;
			rod->dt = deltaTime;
		}

		if (deltaTime < 1e-9)
		{
      if (~addingInertial)
      {
        addingInertial = true;
        deltaTime = 0.1;
      }
      else{
        break;
      }
		}
	}

	cout <<"curvature: " << yForce/rod->EI * dl_b * Lgb<<" normalized ls: "<<hangLength/Lgb <<" "
	<<"normalized ks: "<<normKs<<" "<<rod->u.norm()<< " "<< deltaTime <<" "<<searchStep<<endl;


	rod->updateTimeStep();

  // record force
	switch (searchStep)
	{
		case 1:{
			if (rod->u.norm()<finiteDiff && getGoal)
			{
				// record F0;
				F0 = stepper->force.segment(4, 3);
				m1Local = rod->m1.row(numVertices-2);
				tLocal = rod->tangent.row(numVertices-2);
				searchStep = 2;
				getGoal = false;
				innerIter = 0;
			}
			break;
		}
		case 2:{
			if ( rod->u.norm()<finiteDiff &&getGoal)
			{
				// record F0;
				Fx = stepper->force.segment(4, 3);
				searchStep = 3;
				getGoal = false;
				innerIter = 0;
			}
			break;
		}
		case 3:{
			if (rod->u.norm()< finiteDiff && getGoal)
			{
				// record F0;
				Fy = stepper->force.segment(4, 3);
				searchStep = 4;
				getGoal = false;
				innerIter = 0;
			}
			break;
		}
		case 4:{
			if (rod->u.norm()<finiteDiff && getGoal)
			{
				// record F0;
				Fz = stepper->force.segment(4, 3);
				rate = 0.1;
				run = false;
				searchStep = 5;
				innerIter = 0;
				getGoal = false;
				dx = Vector3d(0, 0, 0);
			}
			break;
		}
		case 5:{
			calSearchDirection();
			break;
		}
	}
	// cout << searchStep << " "<<currentTime << " " <<rod->u.norm()<< endl;


	currentTime += deltaTime;

	timeStep++;

	if (deltaTime < 1e-1 && iter <= 5)
	{
		deltaTime = 10 * deltaTime;
		rod->dt = deltaTime;
	}

	if (solved == false)
	{
		currentTime = totalTime; // we are exiting
		cout <<"solved false issue" << endl;
	}
}


void world::calSearchDirection()
{
	run = false;
	Vector3d FRef = F0;
	FRef(1) = FRef(1) - yForce;
	if (dx.norm() == 0)
	{
		Jc.col(0) = (Fx - F0)/finiteDiff;
		Jc.col(1) = (Fy - F0)/finiteDiff;
		Jc.col(2) = (Fz - F0)/finiteDiff;

		dx = Jc.colPivHouseholderQr().solve(FRef);
	}

	cout << rate * dx.norm() <<" "<<FRef.norm()<<" "<<rate*dx.norm()<<endl;

	if ((rate * dx.norm() < 1e-10 && FRef.norm()<1e-6) || FRef.norm() < 1e-10 || rate * dx.norm() < 1e-20 )
	{
		// cout <<"deltaKappa "<<deltaKappa << endl;
		saveData = true;
		searchStep = 1;
		gStep = 0;
		innerIter = 0;
	}

	if (rate* dx.norm() < 1e-2)
	{
		run = true;
	}
	else
	{
		rate *= 0.5;
	}

	if (getGoal)
	{
		Vector3d temp = stepper->force.segment(4, 3);
		temp(1) = (temp(1) - yForce);
		if (temp.norm() >= FRef.norm())
		{
			rate *= 0.5;
		}
		else
		{
			xTop = xTop - dx(0) * rate;
			yTop = yTop - dx(1) * rate;
			zTop = zTop - dx(2) * rate;

			gStep = gStep + 1;
			searchStep = 1;
		}
	}

}

void world::newtonMethod(bool &solved)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	iter = 0;
	while (solved == false)
	{
		rod->prepareForIteration();

		stepper->setZero();

		// Compute the forces and the jacobians
    if (addingInertial || addStep == 0)
    {
      m_inertialForce->computeFi();
      m_inertialForce->computeJi();
    }

    // m_inertialForce->computeFi();
    // m_inertialForce->computeJi();

		m_stretchForce->computeFs();
		m_stretchForce->computeJs();

		m_bendingForce->computeFb();
		m_bendingForce->computeJb();

		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		// m_dampingForce->computeFd();
		// m_dampingForce->computeJd();

		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);
		if (isNegativeNaN(normf) || isnan(normf))
		{
			break;
		}
		if (iter == 0)
		{
			normf0 = normf;
		}

		if (normf <= forceTol )
		{
			solved = true;
		}
		else if(iter > 0 && (normf <= normf0 * stol || normf <= 1e-8))
		{
			solved = true;
		}

		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
      cout << normf << endl;
			break;
		}
	}
}

int world::simulationRunning()
{
	if (currentTime<totalTime)
		return 1;
	else
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] /RodLength ;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}
