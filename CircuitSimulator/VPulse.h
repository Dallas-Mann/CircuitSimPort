#ifndef VPULSE_H
#define VPULSE_H

#include "Component.h"
#include "Netlist.h"

class VPulse : public Component {
public:
	VPulse(string id, int nodeOne, int nodeTwo, double amplitude, double riseTime, double pulseWidth, double minTime, double maxTime, int numSteps);
	~VPulse() = default;
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
	double pulseWidth;
	double minTime;
	double maxTime;
	int numSteps;
};
#endif