#ifndef VCCS_H
#define VCCS_H

#include "Component.h"
#include "Netlist.h"

class VCCS : public Component {
public:
	VCCS(string id, int nodeOne, int nodeTwo, int nodeThree, int nodeFour, double gain);
	~VCCS() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
private:
	string id;
	int nodeOne;
	int nodeTwo;
	int nodeThree;
	int nodeFour;
	double gain;
};
#endif