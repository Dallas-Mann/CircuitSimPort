#include "IndCurrentSource.h"

IndCurrentSource::IndCurrentSource(string id, int nodeOne, int nodeTwo, double current) :
	id{ id }, nodeOne{ nodeOne }, nodeTwo{ nodeTwo }, current{ current } {
}

void IndCurrentSource::insertStamp(vector<T>& GVals, vector<T>& XVals, vector<T>& CVals, vector<T>& BVals) {
	// 0th node is ground node, and thus not implemented in our matrices
	// because of this we need to offset all the matrix indices by -1
	int indexOne = nodeOne - 1;
	int indexTwo = nodeTwo - 1;
	if (!nodeOne == 0) {
		BVals.push_back({ indexOne, 0, -current });
	}
	if (!nodeTwo == 0) {
		BVals.push_back({ indexTwo, 0, current });
	}
}

int IndCurrentSource::numVoltagesToAdd(vector<int> &nodes) {
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

int IndCurrentSource::numCurrentsToAdd() {
	return 0;
}
