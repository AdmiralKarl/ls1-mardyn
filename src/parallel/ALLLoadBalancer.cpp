/**
 * @file ALLLoadBalancer.cpp
 * @author seckler
 * @date 04.06.19
 */

#include "ALLLoadBalancer.h"
#include "parallel/DomainDecompMPIBase.h"

ALLLoadBalancer::ALLLoadBalancer(std::array<double, DIMgeom> localBoxMin, std::array<double, DIMgeom> localBoxMax, double gamma,
								 MPI_Comm comm, std::array<int, 3> globalSize, std::vector<double> minimalPartitionSize) {
	
	_localBoxMin = localBoxMin;
	_localBoxMax = localBoxMax;
	
	_coversWholeDomain = {globalSize[0] == 1, globalSize[1] == 1, globalSize[2] == 1};;

	_all = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, DIMgeom, gamma);
	_all->setCommunicator(comm);
	
	_all->setMinDomainSize(minimalPartitionSize);
    _all->setup();

}
std::tuple<std::array<double, DIMgeom>, std::array<double, DIMgeom>> ALLLoadBalancer::rebalance(double work) {
	std::vector<ALL::Point<double>> domain(2, ALL::Point<double>(DIMgeom));

	for (int i = 0; i < DIMgeom; ++i) {
		domain[0][i] = _localBoxMin[i];
		domain[1][i] = _localBoxMax[i];
	}

	_all->setVertices(domain);
	_all->setWork(work);
	_all->balance();

	std::vector<ALL::Point<double>> updatedVertices = _all->getVertices();

	for (int i = 0; i < DIMgeom; ++i) {
		_localBoxMin[i] = updatedVertices[0][i];
		_localBoxMax[i] = updatedVertices[1][i];
	}

	return std::make_tuple(_localBoxMin, _localBoxMax);
}
