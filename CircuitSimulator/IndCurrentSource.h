#ifndef INDCURRENTSOURCE_H
#define INDCURRENTSOURCE_H

#include "Component.h"
#include "Netlist.h"

class IndCurrentSource : public Component {
public:
	IndCurrentSource(string id, int nodeOne, int nodeTwo, double current);
	~IndCurrentSource() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
private:
	string id;
	int nodeOne;
	int nodeTwo;
	double current;
};
#endif