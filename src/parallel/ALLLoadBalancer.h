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
	void readXML(XMLfileUnits& xmlconfig) override;

	std::array<bool, 3> getCoversWholeDomain() override { return _coversWholeDomain; }

	std::tuple<std::array<double, 3>, std::array<double, 3>> getlocalDomain() override;

	void setlocalDomain(std::array<double, 3> newBoxMin, std::array<double, 3> newBoxMax) override;

private:
	std::unique_ptr<ALL::ALL<double, double>> _all;
	MPI_Comm _comm;
	double _gamma;

	std::array<double, 3> _localBoxMin;
	std::array<double, 3> _localBoxMax;

	std::vector<double> _minimalPartitionSize{};
	std::array<bool, 3> _coversWholeDomain{};
};
#endif
