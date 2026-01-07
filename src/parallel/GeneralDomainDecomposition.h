/**
 * @file GeneralDomainDecomposition.h
 * @author georg
 * @date 6.11.25
 */

#pragma once

#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "DomainDecompMPIBase.h"
#include "LoadBalancer.h"



/**
 * This decomposition is meant to be able to call arbitrary load balancers.
 */
class GeneralDomainDecomposition : public DomainDecompMPIBase {
public:
    /**
	 * Constructor for the GeneralDomainDecomposition.
	 * @param interactionLength
	 * @param domain
	 * @param forceGrid
	 */
	//GeneralDomainDecomposition(double interactionLength, Domain* domain, bool forceGrid);

	GeneralDomainDecomposition(double interactionLength, Domain* domain);

	GeneralDomainDecomposition(double interactionLength, Domain* domain, MPI_Comm comm);

	// documentation see father class (DomainDecompBase.h)
	~GeneralDomainDecomposition() override;

	// Read in XML configuration for GeneralDomainDecomposition and all its included objects.
	void readXML(XMLfileUnits& xmlconfig) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer,
							Domain* domain) override;

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanks() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
											double etime) override {
		throw std::runtime_error(
			"GeneralDomainDecomposition::queryBalanceAndExchangeNonBlocking() not yet implemented");
	}

	std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion,
																  double cutoff) override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighboursFromHaloRegion() not yet implemented");
	}
private:
    /**
	 * Method that initializes the ALLLoadBalancer
	 */
	void initializeALLLoadBalancer();
	
	/**
	Creates a new MPI communicator with topology information added.
	*/
	void initMPIGridDims();

    std::tuple<std::array<double, 3>, std::array<double, 3>> initializeRegularGrid(const std::array<double, DIMgeom>& domainLength, const std::array<int, DIMgeom>& gridSize,
	const std::array<int, DIMgeom>& gridCoords);
	

	/**
	 * Initializes communication partners
	 * @param moleculeContainer
	 * @param domain
	 */
	void initCommunicationPartners(Domain* domain, ParticleContainer* moleculeContainer);

	/**
	 * Calculate new distribution on process and migrate accordingly.
	 * @param lastTraversalTime time of last calculation
	 * @param domain
	 * @param particleContainer
	 */
	void rebalance(double lastTraversalTime, ParticleContainer* moleculeContainer, Domain* domain);

	/**
	 * Exchange the particles, s.t., particles are withing the particleContainer of the process they belong to.
	 * This function will rebuild the particleContainer.
	 * @param domain
	 * @param particleContainer
	 * @param newMin new minimum of the own subdomain
	 * @param newMax new maximum of the own subdomain
	 */
	void migrateParticles(Domain* domain, ParticleContainer* particleContainer, std::array<double, 3> newMin,
						  std::array<double, 3> newMax);

	/**
	 * Check whether a rebalancing is necessary.
	 * @param step current step of the simulation
	 */
	bool checkRebalancing(size_t step);

    // variables
	const bool debugMode = false;
	bool _cartCommunicatorCreated = false; // Indicates whether a communicator with topology information has already been created.

	std::array<double, DIMgeom> _boxMin;
	std::array<double, DIMgeom> _boxMax;

	std::array<double, 3> _domainLength;
	double _interactionLength;

	size_t _steps{0};
	size_t _rebuildFrequency{10000};

	size_t _initPhase{0};
	size_t _initFrequency{500};
	
	// the LoadBalancer used
	std::unique_ptr<LoadBalancer> _loadBalancer{nullptr};
	
	 // Number of processes in each dimension of the MPI process grid
	std::array<int, DIMgeom> _gridSize;
	
	// Coordinate of the process in the MPI process grid
	std::array<int, DIMgeom> _coords;
};