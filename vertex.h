#ifndef VERTEX_H
#define VERTEX_H

#include "global.h"
#include "task.h"

class Vertex {

public:
	typedef hash_map<size_v, Task*> TaskMap;
	typedef typename TaskMap::iterator TaskMapIter;

	size_v _type;  // 0:U set, 1:V set
	size_v id;
	size_v degree;

	vector<size_v> nbs; // neighbor list
	TaskMap tMap;

	Vertex()
	{

	}

	~Vertex() {}

	Task* get_task(size_v task_id)   // designed for parallel processing
	{
		TaskMapIter tIter = tMap.find(task_id);
		if (tIter != tMap.end()) // the task existing
		{
			return tIter->second;
		}
		else
		{
			//Task* task = new Task;
			//tMap[task_id] = task;
			return NULL;
		}
	}

	void free(size_v task_id)
	{
		TaskMapIter iter = tMap.find(task_id);
		if (iter != tMap.end())
		{
			delete iter->second;
			tMap.erase(iter);
		}
	}
};

#endif