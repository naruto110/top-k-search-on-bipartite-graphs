#ifndef TASK_H
#define TASK_H

#include "global.h"

class Task{

public:

	size_v task_id;
	size_v in_degree; //induced degree
	vector<size_v> in_nbs;  // induced neighbors

	Task(size_v task_id) 
	{
		this->task_id = task_id;
		in_degree = 0;
	}

	~Task() {}



};


#endif

