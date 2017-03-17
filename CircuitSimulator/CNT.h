#ifndef CNT_H
#define CNT_H

#include "Component.h"
#include "Netlist.h"

class CNT : public Component {
public:
	CNT(string id, int nodeOne, int nodeTwo, string sourceFile);
	~CNT() = default;
	virtual void insertStamp(vector<T> &GVals, vector<T> &XVals, vector<T> &CVals, vector<T> &BVals) override;
	virtual int numVoltagesToAdd(vector<int> &nodes) override;
	virtual int numCurrentsToAdd() override;
	void setNodeOffset(int offset);
	vector<Component*> subComponents;
	void readNetlist(string inputFilename);
	int modifyNodeValues();
	string sourceFile;

private:
	string id;
	int nodeOne;
	int nodeTwo;
	int nodeOffset;
};
#endif