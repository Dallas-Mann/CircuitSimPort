#ifndef VSTEP_H
#define VSTEP_H

#include "Component.h"
#include "Netlist.h"

class VStep : public Component {
public:
	VStep(string id, int nodeOne, int nodeTwo, double amplitude, double riseTime, double minTime, double maxTime, int numSteps);
	~VStep() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNewIndex(int val);

	string id;
	int nodeOne;
	int nodeTwo;
	int newIndex;
	double amplitude;
	double riseTime;
	double minTime;
	double maxTime;
	int numSteps;
};
#endif