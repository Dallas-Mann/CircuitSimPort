#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <string>
using std::string;

// Print the correct usage in case of user syntax error.
static void Usage(char* arg, string message) {
	cerr << endl << "Message: " << message << endl;
	cerr << "Usage:   " << arg << " <netlist.txt> <soultion.txt>" << endl;
	exit(-1);
}

static void Error(string message) {
	cerr << endl << "Message: " << message << endl;
	exit(-1);
}
#endif