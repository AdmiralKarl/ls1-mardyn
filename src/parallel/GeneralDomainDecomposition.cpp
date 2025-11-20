/**
 * @file GeneralDomainDecomposition.cpp
 * @author georg
 * @date 6.11.25
 */

#include "GeneralDomainDecomposition.h"
#ifdef ENABLE_ALLLBL
#include "ALLLoadBalancer.h"
#endif

#include "Domain.h"
#include "NeighborAcquirer.h"
#include "NeighbourCommunicationScheme.h"

#include "utils/String_utils.h"
#include "utils/mardyn_assert.h"

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <sstream>




GeneralDomainDecomposition::GeneralDomainDecomposition(double interactionLength, Domain* domain) : GeneralDomainDecomposition(interactionLength, domain, MPI_COMM_WORLD) {}

GeneralDomainDecomposition::GeneralDomainDecomposition(double interactionLength, Domain* domain, MPI_Comm comm) : 
DomainDecompMPIBase(comm),
_interactionLength{interactionLength},
_domainLength{domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2)},
_gridSize({0,0,0}), 
_coords{0} {
	initMPIGridDims();
}

void GeneralDomainDecomposition::initMPIGridDims() {
	mardyn_assert(DIMgeom == 3);
	int period[DIMgeom] = {1, 1, 1}; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder = 1; // 1(true) if the ranking may be reordered by MPI_Cart_create

	MPI_CHECK(MPI_Dims_create( _numProcs, DIMgeom, _gridSize.data()));
	MPI_CHECK(MPI_Cart_create(_comm, DIMgeom, _gridSize.data(), period, reorder, &_comm));

	Log::global_log->info() << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1] << ", " << _gridSize[2] << std::endl;
	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));
	MPI_CHECK(MPI_Cart_coords(_comm, _rank, DIMgeom, _coords.data()));
	Log::global_log->info() << "MPI coordinate of current process: " << _coords[0] << ", " << _coords[1] << ", " << _coords[2] << std::endl;

	std::tie(_boxMin, _boxMax) = initializeRegularGrid(_domainLength, _gridSize, _coords);
}

GeneralDomainDecomposition::~GeneralDomainDecomposition() {
	MPI_Comm_free(&_comm);
}


void GeneralDomainDecomposition::initializeALLLoadBalancer() {
	Log::global_log->info() << "initializing ALL load balancer..." << std::endl;
#ifdef ENABLE_ALLLBL
	const std::vector<double> minimalPartitionSize = {_interactionLength, _interactionLength, _interactionLength};
	_loadBalancer = std::make_unique<ALLLoadBalancer>(_boxMin, _boxMax, 4 /*gamma*/, this->getCommunicator(), _gridSize,  minimalPartitionSize);
#else
	std::ostringstream error_message;
	error_message << "ALL load balancing library not enabled. Aborting." << std::endl;
	MARDYN_EXIT(error_message.str());
#endif
	Log::global_log->info() << "GeneralDomainDecomposition initial box: [" << _boxMin[0] << ", " << _boxMax[0] << "] x ["
					   << _boxMin[1] << ", " << _boxMax[1] << "] x [" << _boxMin[2] << ", " << _boxMax[2] << "]"
					   << std::endl;
}

void GeneralDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	// Ensures that the readXML() call to DomainDecompMPIBase forces the direct-pp communication scheme.
	_forceDirectPP = true;

	DomainDecompMPIBase::readXML(xmlconfig);

#ifdef MARDYN_AUTOPAS
	Log::global_log->info() << "AutoPas only supports FS, so setting it." << std::endl;
	setCommunicationScheme("direct-pp", "fs");
#endif

	xmlconfig.getNodeValue("updateFrequency", _rebuildFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition update frequency: " << _rebuildFrequency << std::endl;

	xmlconfig.getNodeValue("initialPhaseTime", _initPhase);
	Log::global_log->info() << "GeneralDomainDecomposition time for initial rebalancing phase: " << _initPhase << std::endl;

	xmlconfig.getNodeValue("initialPhaseFrequency", _initFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition frequency for initial rebalancing phase: " << _initFrequency
					   << std::endl;

	if (xmlconfig.changecurrentnode("loadBalancer")) {
		std::string loadBalancerString = "None";
		xmlconfig.getNodeValue("@type", loadBalancerString);
		Log::global_log->info() << "Chosen Load Balancer: " << loadBalancerString << std::endl;

		std::transform(loadBalancerString.begin(), loadBalancerString.end(), loadBalancerString.begin(), ::tolower);

		if (loadBalancerString.find("all") != std::string::npos) {
			initializeALLLoadBalancer();
		} else {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: Unknown load balancer " << loadBalancerString
								<< ". Aborting! Please select a valid option! Valid options: ALL";
			MARDYN_EXIT(error_message.str());
		}
		_loadBalancer->readXML(xmlconfig);
	} else {
		std::ostringstream error_message;
		error_message << "loadBalancer section missing! Aborting!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	xmlconfig.changecurrentnode("..");
}

double GeneralDomainDecomposition::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

double GeneralDomainDecomposition::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }


void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {
    std::ostringstream error_message;
	error_message << "Not Implementet! Aborting!" << std::endl;
	MARDYN_EXIT(error_message.str());

}



std::tuple<std::array<double, DIMgeom>, std::array<double, DIMgeom>> GeneralDomainDecomposition::initializeRegularGrid(const std::array<double, DIMgeom>& domainLength, 
	const std::array<int, DIMgeom>& gridSize, const std::array<int, DIMgeom>& gridCoords) {
	std::array<double, DIMgeom> boxMin{0.};
	std::array<double, DIMgeom> boxMax{0.};
	// initialize it as regular grid!
	for (int dim = 0; dim < DIMgeom; ++dim) {
		boxMin[dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		boxMax[dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			boxMax[dim] = domainLength[dim];
		}
	}
	return std::make_tuple(boxMin, boxMax);
}
