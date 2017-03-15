#ifndef RESISTOR_H
#define RESISTOR_H

#include "Component.h"
#include "Netlist.h"

class Resistor : public Component {
public:
	Resistor(string id, int nodeOne, int nodeTwo, double resistance);
	~Resistor() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
private:
	string id;
	int nodeOne;
	int nodeTwo;
	double resistance;
};
#endif