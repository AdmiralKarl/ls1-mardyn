/**
 * @file ALLLoadBalancer.h
 * @author seckler
 * @date 04.06.19
 */

#pragma once
#ifdef ENABLE_ALLLBL
#include <ALL.hpp>
#include "LoadBalancer.h"
#include "parallel/DomainDecompMPIBase.h"

#include <tuple>
class ALLLoadBalancer : public LoadBalancer {
public:
	ALLLoadBalancer(std::array<double, DIMgeom> localBoxMin, std::array<double, DIMgeom> localBoxMax, double gamma,
								 MPI_Comm comm, std::array<int, 3> globalSize, std::vector<double> minimalPartitionSize);

	~ALLLoadBalancer() override = default;
	std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) override;
	void readXML(XMLfileUnits& xmlconfig) override {
		// nothing yet.
	}

	std::array<bool, 3> getCoversWholeDomain() override { return _coversWholeDomain; }

private:
	std::unique_ptr<ALL::ALL<double, double>> _all;
	std::array<double, 3> _localBoxMin;
	std::array<double, 3> _localBoxMax;

	std::array<double, 3> _minimalPartitionSize{};
	std::array<bool, 3> _coversWholeDomain{};
};
#endif
