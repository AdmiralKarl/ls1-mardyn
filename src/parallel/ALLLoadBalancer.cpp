/**
 * @file ALLLoadBalancer.cpp
 * @author seckler
 * @date 04.06.19
 */

#include "ALLLoadBalancer.h"
#include <string>
#include "ALL.hpp"
#include "parallel/DomainDecompMPIBase.h"

ALLLoadBalancer::ALLLoadBalancer(std::array<double, DIMgeom> localBoxMin, std::array<double, DIMgeom> localBoxMax, double gamma,
								 MPI_Comm comm, std::array<int, 3> globalSize, std::vector<double> minimalPartitionSize) {
	
	_localBoxMin = localBoxMin;
	_localBoxMax = localBoxMax;
	_comm = comm;
	_gamma = gamma;
	_minimalPartitionSize = minimalPartitionSize;
	
	_coversWholeDomain = {globalSize[0] == 1, globalSize[1] == 1, globalSize[2] == 1};;
}

void ALLLoadBalancer::readXML(XMLfileUnits& xmlconfig){
	ALL::LB_t mode = ALL::LB_t::UNIMPLEMENTED;

	std::string loadBalancer("TENSOR");
	xmlconfig.getNodeValue("mode", loadBalancer);
	
	if (loadBalancer == "STAGGERED") {
		mode = ALL::LB_t::STAGGERED;
	} else if (loadBalancer == "TENSOR") {
		mode = ALL::LB_t::TENSOR;
	} else if (loadBalancer == "FORCEBASED") {
		mode = ALL::LB_t::FORCEBASED; 
	} else if (loadBalancer == "ALL_VORONOI_ACTIVE") {
		#ifdef ALL_VORONOI_ACTIVE
			mode = ALL::LB_t::VORONOI;
		#else
			std::ostringstream error_message;
			error_message << "ALLLoadBalancer: ALL libery has VORONOI not active. Aborting! Please select a valid option!";
			MARDYN_EXIT(error_message.str());
		#endif
	} else if (loadBalancer == "HISTOGRAM") {
		mode = ALL::LB_t::HISTOGRAM;
	} else if (loadBalancer == "TENSOR_MAX") {
		mode = ALL::LB_t::TENSOR_MAX;
	} else {
		std::ostringstream error_message;
		error_message << "ALLLoadBalancer: Unsupported load balancer " << loadBalancer << " was selected. Aborting! Please select a valid option!";
		MARDYN_EXIT(error_message.str());
	}

	Log::global_log->info() << "ALLLoadBalancer: using the " << loadBalancer << " load balancer" << std::endl;

	_all = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, DIMgeom, _gamma);
	_all->setCommunicator(_comm);
	_all->setMinDomainSize(_minimalPartitionSize);
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
