#include "VPulse.h"

VPulse::VPulse(string id, int nodeOne, int nodeTwo, double amplitude, double riseTime, double pulseWidth, double minTime, double maxTime, int numSteps) :
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, amplitude{ amplitude }, riseTime{ riseTime }, 
	pulseWidth{ pulseWidth }, minTime{ minTime }, maxTime{ maxTime }, numSteps{ numSteps } {
}

void VPulse::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	if (!nodeOne == 0) {
		GVals.push_back({ indexOne, newIndex, 1 });
		GVals.push_back({ newIndex, indexOne, 1 });
	}
	if (!nodeTwo == 0) {
		GVals.push_back({ indexTwo, newIndex, -1 });
		GVals.push_back({ newIndex, indexTwo, -1 });
	}
	BVals.push_back({ newIndex, 0, amplitude });
}

int VPulse::numVoltagesToAdd(vector<int> &nodes) {
	int val = 0;
	if (std::find(nodes.begin(), nodes.end(), nodeOne) == nodes.end()) {
		nodes.push_back(nodeOne);
		val++;
	}
	if (std::find(nodes.begin(), nodes.end(), nodeTwo) == nodes.end()) {
		nodes.push_back(nodeTwo);
		val++;
	}
	return val;
}

int VPulse::numCurrentsToAdd() {
	return 1;
}

void VPulse::setNewIndex(int val)
{
	newIndex = val;
}
