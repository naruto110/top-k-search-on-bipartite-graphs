#ifndef EGONETWORK_H
#define EGONETWORK_H

#include "global.h"
#include "graph.h"
#include <omp.h>
#include "./parallel_hashmap/phmap.h"
using phmap::flat_hash_map;
//#include "oneapi/tbb/concurrent_hash_map.h"
//#include "oneapi/tbb/concurrent_vector.h"


class Ego{
public:
	typedef hash_map<size_v, vector<size_v> > NeiMap;  // induced neighbors
	typedef typename NeiMap::iterator NeiIter;

	//typedef hash_map<size_v, size_v> BetweenMap;  // induced neighbors
	//typedef typename BetweenMap::iterator BetweenIter;

	//typedef hash_map<size_v, size_v> BetweenMap;  // induced neighbors
	typedef flat_hash_map<size_v, size_v> BetweenMap;  // induced neighbors

	//typedef typename BetweenMap::iterator BetweenIter;

	//typedef oneapi::tbb::concurrent_hash_map<size_v, size_v> BetweenMap;
	//typedef typename BetweenMap::iterator BetweenIter;

	//typedef hash_map<size_v, vector<size_v> > EgoMap;  // induced neighbors
	typedef flat_hash_map<size_v, vector<size_v> > EgoMap;  // induced neighbors
	typedef typename EgoMap::iterator EgoIter;

	size_v uid;  // src
	size_v vid;  // des
	size_v task_id;
	size_v in_edge_num;
	size_v global_idx;
	size_b betweenness;

	size_v max_vid;
	size_v max_uid;

	vector<size_v> uset;
	vector<size_v> vset;

	// thread safe
	EgoMap ego_map_u;
	EgoMap ego_map_v;

	size_v* uset_arr;
	size_v* vset_arr;


	Ego(size_v src, size_v des, size_v task_id)
	{
		this->uid = src;
		this->vid = des;
		this->task_id = task_id;
		this->in_edge_num = 0;
		this->betweenness = 0.0;
	}

	Ego(Edge* edge)
	{
		this->uid = edge->src;
		this->vid = edge->des;
		this->task_id = edge->edge_id;
		this->in_edge_num = 0;
		this->betweenness = 0.0;
		this->global_idx = edge->global_idx;
	}

	~Ego() {
	
		free(uset_arr);
		free(vset_arr);
	}

	void set_sets(Graph graph)
	{
		size_v upos = graph.get_vposU(uid);
		size_v vpos = graph.get_vposV(vid);

		vset = graph.vertexesU[upos]->nbs;
		uset = graph.vertexesV[vpos]->nbs;

		max_uid = uset[0];
		max_vid = vset[0];
		for (size_v uidx = 0; uidx < uset.size(); uidx++)
		{
			if (uset[uidx] > max_uid)
			{
				max_uid = uset[uidx];
			}
		}
		for (size_v vidx = 0; vidx < vset.size(); vidx++)
		{
			if (vset[vidx] > max_vid)
			{
				max_vid = vset[vidx];
			}
		}

	}

	void generate_ego(Graph &graph)
	{
	#ifdef _WIN32
	#else
		ResetTimer(EGOGENERATE_TIMER);
	#endif
		// step 1: initialize induced vertices
		for (size_v uidx = 0; uidx < uset.size(); uidx++)
		{
			size_v uid = uset[uidx];
			size_v upos = graph.get_vposU(uid);
			Vertex* u = graph.vertexesU[upos];
			Task* t = new Task(task_id);
			u->tMap[task_id] = t;
		}

		for (size_v vidx = 0; vidx < vset.size(); vidx++)
		{
			size_v vid = vset[vidx];
			size_v vpos = graph.get_vposV(vid);
			Vertex* v = graph.vertexesV[vpos];
			Task* t = new Task(task_id);
			v->tMap[task_id] = t;
		}

		// step 2: construct induced edges
		for (size_v uidx = 0; uidx < uset.size(); uidx++)
		{
			size_v uid = uset[uidx];
			size_v upos = graph.get_vposU(uid);
			Vertex* u = graph.vertexesU[upos];
			Task* tmp_task_u = u->get_task(task_id);
			if (tmp_task_u == NULL)
			{
				cout << "error in ego-network construction !!!" << endl;
				exit;
			}
			for (size_v neiidx = 0; neiidx < u->nbs.size(); neiidx++)
			{
				size_v nei_vid = u->nbs[neiidx];
				size_v nei_vpos = graph.get_vposV(nei_vid);  // u's neighbors are all from V set
				Vertex* neiv = graph.vertexesV[nei_vpos];
				Task* tmp_task = neiv->get_task(task_id);
				if (tmp_task != NULL)  // nei_vid is an induced neighbor
				{
					tmp_task->in_nbs.push_back(uid);
					tmp_task->in_degree += 1;

					tmp_task_u->in_nbs.push_back(nei_vid);
					tmp_task_u->in_degree += 1;

					this->in_edge_num++;
				}
			}
		}

		// subgraph construct
		for (size_v uidx = 0; uidx < uset.size(); uidx++)
		{
			size_v uid = uset[uidx];
			size_v upos = graph.get_vposU(uid);
			Vertex* u = graph.vertexesU[upos];
			vector<size_v> tmp_nei = u->tMap[task_id]->in_nbs;
			ego_map_u[uid] = tmp_nei;
			u->free(task_id);
		}

		for (size_v vidx = 0; vidx < vset.size(); vidx++)
		{
			size_v vid = vset[vidx];
			size_v vpos = graph.get_vposV(vid);
			Vertex* v = graph.vertexesV[vpos];
			vector<size_v> tmp_nei = v->tMap[task_id]->in_nbs;
			ego_map_v[vid] = tmp_nei;
			v->free(task_id);
		}

		//uset_arr = (size_v*)malloc((max_uid + 1) * sizeof(size_v));

	#ifdef _WIN32
	#else
		StopTimer(EGOGENERATE_TIMER);
		//global_egogen_time += get_timer((EGOGENERATE_TIMER));
		//PrintTimer("Ego generate time: ", EGOGENERATE_TIMER);
	#endif
		cout << "Edge (" << uid << ", " << vid << ") ego generate ok !!! (" << uset.size() << ", " << vset.size() << ") number edges: " << in_edge_num << endl;
	}

	size_b cal_egobetweenness_opt(Graph graph)
	{
#ifdef _WIN32
#else
		ResetTimer(BETCAL_TIMER);
#endif
		if (uset.size() < vset.size())// start from vset  opt: <  baseline: >
		{
			#pragma omp parallel for num_threads(global_total_threads)
			//	#pragma omp parallel for num_threads(1)
			for (size_v vidx = 0; vidx < vset.size(); vidx++)
			{
				//cout << omp_get_thread_num() << " to " << vidx << endl;
				size_v tmpvid = vset[vidx];
				vector<size_v> v_nbs = ego_map_v[tmpvid];

				size_b maxBetweenness = uset.size();
				BetweenMap bMap;

				for (size_v uidx = 0; uidx < uset.size(); uidx++)
				{
					bMap[uset[uidx]] = 0;
				}
				for (size_v nbidx = 0; nbidx < v_nbs.size(); nbidx++)
				{
					bMap[v_nbs[nbidx]] = -1; // in the end, plus 1 to maxBetweenness, for the original edge
				}

				//#pragma omp parallel for
				for (size_v nbidx = 0; nbidx < v_nbs.size(); nbidx++)
				{
					// test openmp
				//	cout << omp_get_thread_num() << " to " << nbidx << endl;

					size_v _1hop_nei = v_nbs[nbidx];  // u

					vector<size_v> _2hop_nbs_vec = ego_map_u[_1hop_nei];
					//#pragma omp parallel for
					for (size_v _2hop_idx = 0; _2hop_idx < _2hop_nbs_vec.size(); _2hop_idx++)
					{
						size_v _2hop_nei = _2hop_nbs_vec[_2hop_idx]; // v

						if (_2hop_nei == tmpvid)
						{
							continue;
						}

						vector<size_v> _3hop_nbs_vec = ego_map_v[_2hop_nei];
						//#pragma omp parallel for
						for (size_v _3hop_idx = 0; _3hop_idx < _3hop_nbs_vec.size(); _3hop_idx++)
						{
							size_v _3hop_nei = _3hop_nbs_vec[_3hop_idx]; //u

							if (_3hop_nei == _1hop_nei)
							{
								continue;
							}

							if (bMap[_3hop_nei] != -1)
							{
								size_v tmp_path_num = bMap[_3hop_nei];
								tmp_path_num += 1;
								//#pragma omp critical
								bMap[_3hop_nei] = tmp_path_num;
							}
						}
					}
				}

				for (size_v uidx = 0; uidx < uset.size(); uidx++)
				{
					size_v tmpval = bMap[uset[uidx]];
					if (bMap[uset[uidx]] == -1)
					{
						maxBetweenness -= 1.0;
					}
					else
					{
						if (bMap[uset[uidx]] <= 0)
						{
							cout << "error in vset processing !!!" << endl;
							cout << "bMap[uset[uidx]]: " << bMap[uset[uidx]] << endl;
							while (1)
							{
								;
							}
						}
						size_b tmp_val = 1.0 - 1.0 / (size_b)bMap[uset[uidx]];
						maxBetweenness -= tmp_val;
					}
				}

				//maxBetweenness += 1;
#pragma omp atomic
				betweenness += maxBetweenness;
				//	cout << "contributer: " << uid << " donation: " << maxBetweenness << endl;

			}
			betweenness += 1;
			//cout << "from vset" << endl;
		}
		else  // start from uset
		{
            #pragma omp parallel for num_threads(global_total_threads)
			//	#pragma omp parallel for num_threads(1)
			for (size_v uidx = 0; uidx < uset.size(); uidx++)
			{
				size_v tmpuid = uset[uidx];
				vector<size_v> u_nbs = ego_map_u[tmpuid];

				size_b maxBetweenness = vset.size();
				BetweenMap bMap;

				for (size_v vidx = 0; vidx < vset.size(); vidx++)
				{
					bMap[vset[vidx]] = 0;
				}
				for (size_v nbidx = 0; nbidx < u_nbs.size(); nbidx++)
				{
					bMap[u_nbs[nbidx]] = -1; // in the end, plus 1 to maxBetweenness, for the original edge

					// dynamic upper bound update
					if (global_dynamic_ctr == 1)
					{
						if (tmpuid != uid && u_nbs[nbidx] != vid)
						{
							EdgeIdx eidx;
							eidx.src = tmpuid;
							eidx.des = u_nbs[nbidx];
							Edge* tmpedge = global_edge_map[eidx];
							if (tmpedge->processed == 0 && tmpedge->global_idx < global_idx + 100)  // unprocessed
							{
								tmpedge->ex_edge_num += 1;
								tmpedge->set_up_bound();
								update_edge_map(global_edge_vec, tmpedge);
							}
						}
					}
				}

				for (size_v nbidx = 0; nbidx < u_nbs.size(); nbidx++)
				{
					size_v _1hop_nei = u_nbs[nbidx];  // v
					vector<size_v> _2hop_nbs_vec = ego_map_v[_1hop_nei];
					//#pragma omp parallel for
					for (size_v _2hop_idx = 0; _2hop_idx < _2hop_nbs_vec.size(); _2hop_idx++)
					{
						size_v _2hop_nei = _2hop_nbs_vec[_2hop_idx]; // u

						if (_2hop_nei == tmpuid)
						{
							continue;
						}

						vector<size_v> _3hop_nbs_vec = ego_map_u[_2hop_nei];
						for (size_v _3hop_idx = 0; _3hop_idx < _3hop_nbs_vec.size(); _3hop_idx++)
						{
							size_v _3hop_nei = _3hop_nbs_vec[_3hop_idx]; //v

							if (_3hop_nei == _1hop_nei)
							{
								continue;
							}

							if (bMap[_3hop_nei] != -1)
							{
								size_v tmp_path_num = bMap[_3hop_nei];
								tmp_path_num += 1;
								bMap[_3hop_nei] = tmp_path_num;
							}
						}
					}
				}

				for (size_v vidx = 0; vidx < vset.size(); vidx++)
				{
					if (bMap[vset[vidx]] == -1)
					{
						maxBetweenness -= 1.0;
					}
					else
					{
						if (bMap[vset[vidx]] <= 0)
						{
							cout << "error in vset processing !!!" << endl;
							cout << "bMap[vset[vidx]]: " << bMap[vset[vidx]] << endl;
							while (1)
							{
								;
							}
						}

						size_b tmp_val = 1.0 - 1.0 / (size_b)bMap[vset[vidx]];
						maxBetweenness -= tmp_val;
					}
				}

				//maxBetweenness += 1;
#pragma omp atomic
				betweenness += maxBetweenness;
				//	cout << "contributer: " << uid << " donation: " << maxBetweenness << endl;
			}
			betweenness += 1;
			//cout << "from uset" << endl;
		}
		//	cout << "Edge (" << uid << ", " << vid << ") " << "ego-betweenness: " << betweenness << endl;
#ifdef _WIN32
#else
		StopTimer(BETCAL_TIMER);
		//global_betcal_time += get_timer((BETCAL_TIMER));
		//PrintTimer("Betweenness calculate time: ", BETCAL_TIMER);
#endif
		return betweenness;
	}

};


#endif

