#ifndef INDVOLTAGESOURCE_H
#define INDVOLTAGESOURCE_H

#include "Component.h"
#include "Netlist.h"

class IndVoltageSource : public Component {
public:
	IndVoltageSource(string id, int nodeOne, int nodeTwo, double voltage);
	~IndVoltageSource() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNewIndex(int val);
private:
	string id;
	int nodeOne;
	int nodeTwo;
	double voltage;
	int newIndex;
};
#endif