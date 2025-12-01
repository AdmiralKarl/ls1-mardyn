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

bool GeneralDomainDecomposition::checkRebalancing(size_t step) {
	return step <= _initPhase ? step % _initFrequency == 0 : step % _rebuildFrequency == 0;
}

void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {							
	const bool doRebalance = checkRebalancing(_steps) or forceRebalancing;
	
	if (_steps == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		initCommunicationPartners(domain, moleculeContainer);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		++_steps;
		return;
	}

	if (doRebalance) {
		rebalance(lastTraversalTime, moleculeContainer, domain);
	} else {
		if (sendLeavingWithCopies()) {
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving and Halos." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
		} else {
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
			#ifndef MARDYN_AUTOPAS
				moleculeContainer->deleteOuterParticles();
			#endif
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Halos." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		}
	}
	_boundaryHandler.setLocalRegion(_boxMin.data(),_boxMax.data());
	_boundaryHandler.updateGlobalWallLookupTable();
	++_steps;		
}

void GeneralDomainDecomposition::initCommunicationPartners(Domain* domain, ParticleContainer* moleculeContainer) { 
	auto coversWholeDomain = _loadBalancer->getCoversWholeDomain();
	for (int d = 0; d < DIMgeom; ++d) {
		_neighbourCommunicationScheme->setCoverWholeDomain(d, coversWholeDomain[d]);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(moleculeContainer->getCutoff(), domain, this,
															 moleculeContainer);
}


void GeneralDomainDecomposition::rebalance(double lastTraversalTime, ParticleContainer* moleculeContainer, Domain* domain) {
	Log::global_log->info() << "GeneralDomainDecomposition: rebalancing..." << std::endl;
	Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving." << std::endl;
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

	moleculeContainer->deleteOuterParticles();

	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "work:" << lastTraversalTime << std::endl;
	Log::global_log->set_mpi_output_root(0);
	auto [newBoxMin, newBoxMax] = _loadBalancer->rebalance(lastTraversalTime);
	
																	
	Log::global_log->debug() << "migrating particles" << std::endl;
	migrateParticles(domain, moleculeContainer, newBoxMin, newBoxMax);

	#ifndef MARDYN_AUTOPAS
			moleculeContainer->update();
	#endif

	_boxMin = newBoxMin;
	_boxMax = newBoxMax;

	Log::global_log->debug() << "updating communication partners" << std::endl;
	initCommunicationPartners(domain, moleculeContainer);
	Log::global_log->debug() << "rebalancing finished" << std::endl;

	Log::global_log->debug() << "GeneralDomainDecomposition: Sending Halos." << std::endl;
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
}




void GeneralDomainDecomposition::migrateParticles(Domain* domain, ParticleContainer* particleContainer,
												  std::array<double, 3> newMin, std::array<double, 3> newMax) {
	std::array<double, 3> oldBoxMin{particleContainer->getBoundingBoxMin(0), particleContainer->getBoundingBoxMin(1),
									particleContainer->getBoundingBoxMin(2)};
	std::array<double, 3> oldBoxMax{particleContainer->getBoundingBoxMax(0), particleContainer->getBoundingBoxMax(1),
									particleContainer->getBoundingBoxMax(2)};

	HaloRegion ownDomain{}, newDomain{};
	for (size_t i = 0; i < 3; ++i) {
		ownDomain.rmin[i] = oldBoxMin[i];
		newDomain.rmin[i] = newMin[i];
		ownDomain.rmax[i] = oldBoxMax[i];
		newDomain.rmax[i] = newMax[i];
		ownDomain.offset[i] = 0;
		newDomain.offset[i] = 0;
	}
	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "migrating from"
						<< " [" << oldBoxMin[0] << ", " << oldBoxMax[0] << "] x"
						<< " [" << oldBoxMin[1] << ", " << oldBoxMax[1] << "] x"
						<< " [" << oldBoxMin[2] << ", " << oldBoxMax[2] << "] " << std::endl;
	Log::global_log->debug() << "to"
						<< " [" << newMin[0] << ", " << newMax[0] << "] x"
						<< " [" << newMin[1] << ", " << newMax[1] << "] x"
						<< " [" << newMin[2] << ", " << newMax[2] << "]." << std::endl;
	Log::global_log->set_mpi_output_root(0);
	std::vector<HaloRegion> desiredDomain{newDomain};
	std::vector<CommunicationPartner> sendNeighbors{}, recvNeighbors{};

	std::tie(recvNeighbors, sendNeighbors) =
		NeighborAcquirer::acquireNeighbors(_domainLength, &ownDomain, desiredDomain, _comm);

	std::vector<Molecule> dummy;
	for (auto& sender : sendNeighbors) {
		sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
						false /*don't use invalid particles*/, true /*do halo position change*/,
						true /*removeFromContainer*/);
	}
	// TODO: Using the autopas resizeBox functions would be more efficient.
	std::vector<Molecule> ownMolecules{};
	ownMolecules.reserve(particleContainer->getNumberOfParticles());
	for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
		ownMolecules.push_back(*iter);
		if(debugMode){
			if (not iter->inBox(newMin.data(), newMax.data())) {
				std::ostringstream error_message;
				error_message
					<< "Particle still in domain that should have been migrated."
					<< "BoxMin: "
					<< particleContainer->getBoundingBoxMin(0) << ", "
					<< particleContainer->getBoundingBoxMin(1) << ", "
					<< particleContainer->getBoundingBoxMin(2) << "\n"
					<< "BoxMax: "
					<< particleContainer->getBoundingBoxMax(0) << ", "
					<< particleContainer->getBoundingBoxMax(1) << ", "
					<< particleContainer->getBoundingBoxMax(2) << "\n"
					<< "Particle: \n" << *iter
					<< std::endl;
				MARDYN_EXIT(error_message.str());
			}
		} 
	}
	particleContainer->clear();
	particleContainer->rebuild(newMin.data(), newMax.data());
	particleContainer->addParticles(ownMolecules);
	bool allDone = false;
	double waitCounter = 30.0;
	double deadlockTimeOut = 360.0;
	double startTime = MPI_Wtime();
	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (auto& sender : sendNeighbors) {
			allDone &= sender.testSend();
		}

		// unpack molecules
		for (auto& recv : recvNeighbors) {
			allDone &= recv.iprobeCount(this->getCommunicator(), this->getMPIParticleType());
			allDone &= recv.testRecv(particleContainer, false);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			Log::global_log->warning() << "GeneralDomainDecomposition::migrateParticles: Deadlock warning: Rank " << _rank
								  << " is waiting for more than " << waitCounter << " seconds" << std::endl;
			waitCounter += 1.0;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			Log::global_log->error() << "GeneralDomainDecomposition::migrateParticles: Deadlock error: Rank " << _rank
								<< " is waiting for more than " << deadlockTimeOut << " seconds" << std::endl;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
			break;
		}
	}
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
