#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include "Component.h"

class Utilities {
public:
	static void Usage(char* arg, string message);
	static void Error(string message);
	static vector<string> parseLine(string& nextLine);
	static Component* parseComponent(vector<string>& tokens);
	static double convert(string& token);
	static bool isANumber(string& token);
	static string* splitString(string& token);
	static bool areSame(double a, double b);
};
#endif