#ifndef COMPONENT_H
#define COMPONENT_H

#include <vector>
using std::vector;

#include <Eigen/Sparse>
using Eigen::Triplet;

typedef Triplet<double> T;

class Component {
public:
	~Component() {};
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) = 0;
	virtual int numVoltagesToAdd(vector<int> &nodes) = 0;
	virtual int numCurrentsToAdd() = 0;
};
#endif