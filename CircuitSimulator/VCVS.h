#ifndef VCVS_H
#define VCVS_H

#include "Component.h"
#include "Netlist.h"

class VCVS : public Component {
public:
	VCVS(string id, int nodeOne, int nodeTwo, int nodeThree, int nodeFour, double gain);
	~VCVS() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNewIndex(int val);
private:
	string id;
	int nodeOne;
	int nodeTwo;
	int nodeThree;
	int nodeFour;
	int newIndex;
	double gain;
};
#endif