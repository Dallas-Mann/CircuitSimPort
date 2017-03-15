#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "Component.h"
#include "Netlist.h"

class Inductor : public Component {
public:
	Inductor(string id, int nodeOne, int nodeTwo, double capacitance);
	~Inductor() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNewIndex(int val);
private:
	string id;
	int nodeOne;
	int nodeTwo;
	double inductance;
	int newIndex;
};
#endif