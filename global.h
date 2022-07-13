#ifndef GLOBAL_H
#define GLOBAL_H

using namespace std;

#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 

#include <iostream>
#include <vector>
#include "fstream"
#include <sstream>
#include <iostream>
#include <map>
#include <cmath>
#include <omp.h>

#ifdef _WIN32
	#include <hash_set>
	#include <hash_map>
#else
	#include <ext/hash_set>
	#include <ext/hash_map>
	#define hash_map __gnu_cxx::hash_map
	#define hash_set __gnu_cxx::hash_set
	#include "unistd.h"
	#include "time.h"
#endif

#include <string.h>
#include <time.h>
#include<algorithm>

#define size_v int 
#define size_b float

struct WorkerParams {
	string input_pathU;
	string input_pathV;
	string task_path;
	size_v k;

	WorkerParams()
	{

	}
};

typedef struct edge
{
	size_v edge_id;
	size_v src;
	size_v des;
	size_v src_degree;
	size_v des_degree;
	size_v ex_edge_num;
	size_v global_idx;  // index in global_edge_vec
	size_b betweenness;
	size_b up_bound;
	bool processed;

	edge()
	{
		ex_edge_num = 0;
		processed = 0;
	}

	void set_up_bound()
	{
		size_v dmax = max(src_degree, des_degree) - 1;
		size_v dmin = min(src_degree, des_degree) - 1;
		size_v D = ex_edge_num / dmax;
		size_v d = ex_edge_num % dmax;

		if (D > dmin)
		{
			cout << "error in bound update !!!" << endl;
		}

		up_bound = ((size_b)src_degree - 1) * ((size_b)des_degree - 1) + 1
			- (size_b)D * dmax - ((size_b)dmin - D) * (1 - (size_b)1 / ((size_b)1 + (size_b)D)) * dmax
			- d - ((size_b)dmax - d) * ((size_b)1 / ((size_b)1 + (size_b)D) - (size_b)1 / ((size_b)1 + (size_b)D + d))
			- ((size_b)dmin - D - 1) * ((size_b)1 / ((size_b)1 + (size_b)D) - (size_b)1 / ((size_b)1 + (size_b)D + 1)) * d;

		//cout << "ex_edge_num: " << ex_edge_num << " D: " << D << " d: " << d << endl;
		//cout << "uid: " << src << " vid: " << des << " upper bound: " << up_bound << endl;
	}

} Edge;

bool edge_cmp(Edge* x, Edge* y)
{
	return x->up_bound > y->up_bound;
}

bool result_cmp(Edge* x, Edge* y)
{
	return x->betweenness > y->betweenness;
}

typedef struct edgeidx
{
	size_v src;
	size_v des;

	bool operator< (const edgeidx& x) const
	{
		if (src < x.src)
		{
			return true;
		}
		else if (src == x.src)
		{
			return des < x.des;
		}
		return src < x.src;
	}
} EdgeIdx;

typedef map<EdgeIdx, Edge*> EdgeMap;
typedef typename EdgeMap::iterator EdgeIter;

typedef map<size_b, size_v> ResultIdxMap;
typedef typename ResultIdxMap::iterator ResultIdxIter;

size_v global_dynamic_ctr = 0;
size_v global_priority_ctr = 0;
size_v global_processed_edge = 0;
size_v global_total_threads = 1;

double global_egogen_time = 0.0;
double global_betcal_time = 0.0;

size_v global_max_vid = 0;
size_v global_max_uid = 0;

EdgeMap global_edge_map;

size_v global_edge_idx = 0; 

vector<Edge*> global_edge_vec;
size_v global_cur_idx = 0;   // the edge processed currently

void update_edge_map(vector<Edge*>& global_edge_vec, Edge* edge)
{
	size_v start_idx = edge->global_idx;
	size_v end_idx = start_idx + 1;
	for (size_v edgeidx = start_idx + 1; edgeidx  < global_edge_vec.size(); edgeidx ++)
	{
		Edge* tmp_edge = global_edge_vec[edgeidx];
		if (tmp_edge->up_bound > edge->up_bound)
		{
			continue;
		}
		else
		{
			end_idx = edgeidx;
			break;
		}
	}
	// update global_edge_vec
	if (end_idx == start_idx + 1)
	{
		return;
	}
	else
	{
		for (size_v edgeidx = start_idx; edgeidx < end_idx - 1; edgeidx++)
		{
			global_edge_vec[edgeidx] = global_edge_vec[edgeidx + 1];
			global_edge_vec[edgeidx]->global_idx = edgeidx;
		}
		global_edge_vec[end_idx - 1] = edge;
		edge->global_idx = end_idx - 1;
	}
}

#endif

