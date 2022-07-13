#ifndef GRAPH_H
#define GRAPH_H

#include "global.h"
#include "vertex.h"


class Graph {

public:
	typedef hash_map<size_v, size_v> Map;
	typedef typename Map::iterator MapIter;

	vector<Vertex*> vertexesU;
	vector<Vertex*> vertexesV;

	Map posU;
	Map posV;

	size_v edge_num;

	Graph()
	{
		edge_num = 0;
	}

	~Graph(){}

	void initU(vector<Vertex*> vertexesU)
	{
		int n = vertexesU.size();
		for (int i = 0; i < n; i++) {
			Vertex* v = vertexesU[i];
			posU[v->id] = i;
		}
	}

	int get_vposU(size_v vertex_id)
	{
		//return -1 if the vertex with the specified id is not found
		MapIter it = posU.find(vertex_id);
		if (it == posU.end()) return -1;
		else return it->second;
	}

	void initV(vector<Vertex*> vertexesV)
	{
		int n = vertexesV.size();
		for (int i = 0; i < n; i++) {
			Vertex* v = vertexesV[i];
			posV[v->id] = i;
		}
	}

	int get_vposV(size_v vertex_id)
	{
		//return -1 if the vertex with the specified id is not found
		MapIter it = posV.find(vertex_id);
		if (it == posV.end()) return -1;
		else return it->second;
	}

	void toVertexU(char* line)
	{
		char* pch;
		pch = strtok(line, " ");
		size_v id = atoi(pch);
		size_v vec_capcity = 0;

		if (global_max_uid < id)
		{
			global_max_uid = id;
		}

		Vertex* u = new Vertex;
		u->_type = 0;
		u->id = id;
		pch = strtok(NULL, " ");
		vec_capcity = atoi(pch);
		u->degree = atoi(pch);  
		edge_num += u->degree;

		u->nbs.reserve(vec_capcity);
		while (pch = strtok(NULL, " "))
		{
			size_v nb = atoi(pch);
			u->nbs.push_back(nb);
		}

		vertexesU.push_back(u);
	}

	void toVertexV(char* line)
	{
		char* pch;
		pch = strtok(line, " ");
		size_v id = atoi(pch);
		size_v vec_capcity = 0;

		if (global_max_vid < id)
		{
			global_max_vid = id;
		}

		Vertex* v = new Vertex;
		v->_type = 0;
		v->id = id;
		pch = strtok(NULL, " ");
		vec_capcity = atoi(pch);
		v->degree = atoi(pch);
		v->nbs.reserve(vec_capcity);
		while (pch = strtok(NULL, " "))
		{
			size_v nb = atoi(pch);
			v->nbs.push_back(nb);
		}

		vertexesV.push_back(v);
	}

	void load_graphU(char* line)
	{
		toVertexU(line);
	}
	void load_graphV(char* line)
	{
		toVertexV(line);
	}

	void load_graph(const WorkerParams& params)
	{
	#ifdef _WIN32
	#else
		ResetTimer(WORKER_TIMER);
	#endif
		
		string file_path;
		string file_pathU;
		string file_pathV;

		file_pathU = params.input_pathU;
		file_pathV = params.input_pathV;

		// load U vertices
		cout << "file_path: " << file_pathU << endl;

		ifstream ori_fileU(file_pathU.c_str());
		string line;
		if (ori_fileU)
		{
			while (getline(ori_fileU, line)) 
			{
				if (line.size() != 0)
				{
					load_graphU(const_cast<char*>(line.c_str()));
				}
			}
			ori_fileU.close();
		}
		else
		{
			cout << params.input_pathU.c_str() << endl;
			cout << "no such file" << endl;
			exit(-1);
		}

		// load V vertices
		ifstream ori_fileV(file_pathV.c_str());
		string lineV;
		if (ori_fileV)
		{
			while (getline(ori_fileV, lineV)) 
			{
				if (lineV.size() != 0)
				{
					load_graphV(const_cast<char*>(lineV.c_str()));
					//	free(line);
				}
			}
			ori_fileV.close();
		}
		else 
		{
			cout << params.input_pathV.c_str() << endl;
			cout << "no such file" << endl;
			exit(-1);
		}

		initU(vertexesU);
		initV(vertexesV);

		cout << "Total vertexesU number: " << vertexesU.size() << endl;
		cout << "Total vertexesV number: " << vertexesV.size() << endl;
		cout << "Total edges number: " << edge_num << endl;

	#ifdef _WIN32
	#else
		StopTimer(WORKER_TIMER);
		PrintTimer("Load Graph Time", WORKER_TIMER);
	#endif

	}
};

#endif


