#ifndef VAC_H
#define VAC_H

#include "Component.h"
#include "Netlist.h"

class VAC : public Component {
public:
	VAC(string id, int nodeOne, int nodeTwo, double amplitude, double frequency, double minFrequency, double maxFrequency, int numSteps);
	~VAC() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNewIndex(int val);

	string id;
	int nodeOne;
	int nodeTwo;
	int newIndex;
	double amplitude;
	double frequency;
	double minFrequency;
	double maxFrequency;
	int numSteps;
};

#endif