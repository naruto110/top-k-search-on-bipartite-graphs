#ifndef TOOLS_H
#define TOOLS_H

#include "global.h"

void print_vector(vector<int> v) {
	if (v.size() > 0) {
		cout << "{";
		for (int i = 0; i<int(v.size()); i++) {
			cout << v[i] << ",";
		}
		cout << "}";
	}
	else {
		cout << "{}";
	}
	cout << endl;
}

#endif

