#include "Limiter.h"

Limiter::Limiter(Type limiter, std::shared_ptr<Zone> zone)
{
	_limiter = limiter;
	_polyOrder = zone->getPolyOrder();
	_num_cell = zone->getGrid()->getNumCell();
	_size_cell = zone->getGrid()->getSizeX();
}

Limiter::~Limiter()
{

}

void Limiter::hMLP_Limiter(std::shared_ptr<Zone> zone)
{
	// No limiter if PO
	if ((_polyOrder == 0) || (_limiter == "none")) return;

	// Temporary variables for hMLP
	// Current projected degree
	std::vector<int_t> projectDegree(_num_cell, _polyOrder);

	// hMLP limiting process
	for (int_t step = 0; step < _polyOrder; ++step)
	{
		// Troubled-cell marker
		std::vector<bool> marker(_num_cell, true);

		// Marking troubled-cell
		for (int_t icell = GHOST; icell < _num_cell - GHOST; ++icell)
		{
			marker[icell] = troubleCellMarker(zone, projectDegree[icell], icell);
		}

		// Project troubled-cell
		troubleCellProject(zone, projectDegree, marker);
	}

	// Update zone
	zone->calSolution();
}

void Limiter::troubleCellProject
(std::shared_ptr<Zone> zone, std::vector<int_t>& degree, const std::vector<bool>& marker)
{
	// Temporary DOF and Zone
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Projection
	for (int_t icell = 0; icell < _num_cell; ++icell)
	{
		if (marker[icell] == false)
		{
			if (degree[icell] > 2)
			{
				temp_DOF[degree[icell]][icell] = 0.0;
				degree[icell]--;
			}
			else if (degree[icell] == 2)
			{
				temp_DOF[2][icell] = 0.0;
				degree[icell]--;
			}
			else if (degree[icell] == 1)
			{
				temp_DOF[1][icell] = MLP_limit_ftn(zone, icell)*temp_DOF[1][icell];
			}
			else ERROR("step degree");
		}
	}

	// Update zone
	zone->setDOF(temp_DOF);
	zone->calSolution();
}

bool Limiter::troubleCellMarker(std::shared_ptr<Zone> zone, int_t degree, int_t icell) const
{
	bool marker = augmentMLPmarker(zone, icell);
	if ((marker == false) && (degree > 1)) marker = extremaDetector(zone, icell);
	
	return marker;
}

bool Limiter::augmentMLPmarker(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Augmented MLP condition marker
	bool marker = true;

	// Temporary DOF
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();

	// Variables
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;

	// Variables to MLP condition
	real_t max_avgQ; /// maximum averaged Q
	real_t min_avgQ; /// minimum averaged Q
	real_t max_appQ; /// maximum approximated Q
	real_t min_appQ; /// minimum approximated Q
	real_t appQ; /// approximated vertex Q

	// Left cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_left);
	max_appQ = std::max(zone->getPolySolution(icell - 1, coord_x_left), zone->getPolySolution(icell, coord_x_left));
	min_appQ = std::min(zone->getPolySolution(icell - 1, coord_x_left), zone->getPolySolution(icell, coord_x_left));
	max_avgQ = std::max(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]);
	min_avgQ = std::min(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]);
	if (!((max_avgQ > max_appQ) && (min_appQ > min_avgQ)))
		marker = false;

	// Right cell MLP condition
	appQ = zone->getPolySolution(icell, coord_x_right);
	max_appQ = std::max(zone->getPolySolution(icell, coord_x_right), zone->getPolySolution(icell + 1, coord_x_right));
	min_appQ = std::min(zone->getPolySolution(icell, coord_x_right), zone->getPolySolution(icell + 1, coord_x_right));
	max_avgQ = std::max(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]);
	min_avgQ = std::min(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]);
	if (!((max_avgQ > max_appQ) && (min_appQ > min_avgQ)))
		marker = false;

	return marker;
}

bool Limiter::extremaDetector(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Decomposing the DG-Pn approximation
	real_t avgQ = zone->getDOF()[0][icell];
	real_t max_avgQ;
	real_t min_avgQ;
	real_t P1_projected;
	real_t Pn_projected_slope;
	real_t P1_filtered_Pn;
	
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;

	real_t leftQ = zone->getPolySolution(icell, coord_x_left);
	real_t rightQ = zone->getPolySolution(icell, coord_x_right);

	// Deactivation threshold
	real_t threshold = std::max(0.001*avgQ, _size_cell);
	if ((abs(leftQ - avgQ) <= threshold) && (abs(rightQ - avgQ) <= threshold))
		return true; /// deactivation

	// Left vertex
	bool leftMarker = false;
	P1_projected = projectionTo(1, zone, icell, coord_x_left);
	Pn_projected_slope = P1_projected - avgQ;
	P1_filtered_Pn = zone->getPolySolution(icell, coord_x_left) - P1_projected;
	max_avgQ = std::max(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]);
	min_avgQ = std::min(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]);
	
	// Marking
	if ((Pn_projected_slope > 0.0) && (P1_filtered_Pn < 0.0) && (leftQ > min_avgQ)) leftMarker = true;
	if ((Pn_projected_slope < 0.0) && (P1_filtered_Pn > 0.0) && (leftQ < max_avgQ)) leftMarker = true;

	// Right vertex
	bool rightMarker = false;
	P1_projected = projectionTo(1, zone, icell, coord_x_right);
	Pn_projected_slope = P1_projected - avgQ;
	P1_filtered_Pn = zone->getPolySolution(icell, coord_x_right) - P1_projected;
	max_avgQ = std::max(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]);
	min_avgQ = std::min(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]);

	// Marking
	if ((Pn_projected_slope > 0.0) && (P1_filtered_Pn < 0.0) && (rightQ > min_avgQ)) rightMarker = true;
	if ((Pn_projected_slope < 0.0) && (P1_filtered_Pn > 0.0) && (rightQ < max_avgQ)) rightMarker = true;

	// Smooth extrema detect
	return (leftMarker && rightMarker);
}

real_t Limiter::MLP_limit_ftn(std::shared_ptr<Zone> zone, int_t icell) const
{
	// Variables
	real_t limit_ftn_left;
	real_t limit_ftn_right;
	real_t coord_x_left = zone->getGrid()->getCell()[icell]->getPosX() - 0.5*_size_cell;
	real_t coord_x_right = zone->getGrid()->getCell()[icell]->getPosX() + 0.5*_size_cell;
	real_t avgQ = zone->getDOF()[0][icell];
	real_t del_m = projectionTo(1, zone, icell, coord_x_right) - avgQ;

	// Compute MLP function
	if (del_m > epsilon)
	{
		limit_ftn_right = limit_PI(std::max(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]) - avgQ, del_m);
		limit_ftn_left = limit_PI(std::min(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]) - avgQ, -del_m);
	}
	else if (del_m < -epsilon)
	{
		limit_ftn_right = limit_PI(std::min(zone->getDOF()[0][icell], zone->getDOF()[0][icell + 1]) - avgQ, del_m);
		limit_ftn_left = limit_PI(std::max(zone->getDOF()[0][icell - 1], zone->getDOF()[0][icell]) - avgQ, -del_m);
	}
	else return 1.0;

	return std::min(limit_ftn_right, limit_ftn_left);
}

real_t Limiter::limit_PI(real_t del_p, real_t del_m) const
{
	if (_limiter == "MLP-u1") return std::min(1.0, del_p / del_m);
	if (_limiter == "MLP-u2") return MLP_u2(del_p, del_m);
	ERROR("cannot find limiter");
	return 0;
}

real_t Limiter::MLP_u2(real_t del_p, real_t del_m) const
{
	real_t ep = pow(CONST_K*_size_cell, 1.5);
	real_t num = (pow(del_p, 2.0) + pow(ep, 2.0))*del_m + 2.0*pow(del_m, 2.0)*del_p;
	real_t den = del_m*(pow(del_p, 2.0) + 2.0*pow(del_m, 2.0) + del_m*del_p + pow(ep, 2.0));

	return num / den;
}

real_t Limiter::projectionTo(int_t degree, std::shared_ptr<Zone> zone, int_t icell, real_t coord_x) const
{
	// Temporary object
	std::shared_ptr<Zone> temp_zone = std::make_shared<Zone>(*zone);
	std::vector<std::vector<real_t> > temp_DOF = zone->getDOF();
	std::vector<real_t> zero_DOF(_num_cell, 0.0);
	
	// Projection to n degree
	for (int_t idegree = _polyOrder; idegree > degree; --idegree)
	{
		temp_DOF[idegree] = zero_DOF;
	}

	temp_zone->setDOF(temp_DOF);

	return temp_zone->getPolySolution(icell, coord_x);
}