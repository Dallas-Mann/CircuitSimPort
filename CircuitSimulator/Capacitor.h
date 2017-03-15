#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "Component.h"
#include "Netlist.h"

class Capacitor : public Component {
public:
	Capacitor(string id, int nodeOne, int nodeTwo, double capacitance);
	~Capacitor() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
private:
	string id;
	int nodeOne;
	int nodeTwo;
	double capacitance;
};
#endif