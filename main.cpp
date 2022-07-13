#include "global.h"
#include "graph.h"
#include "tools.h"
#include "egonetwork.h"
#include <omp.h>

#ifdef _WIN32
#include <iostream>
#include <fstream> // c++文件操作
#include <iomanip> // 设置输出格式
#else
//#include "./parallel_hashmap/phmap.h"
#endif


vector<string> in_path_U;
vector<string> in_path_V;
vector<string> task_path;


void set_path(vector<string>& in_path_U, vector<string>& in_path_V, vector<string>& task_path)
{
	in_path_U.push_back("./data/demoU.txt");  // 0: demo
	in_path_V.push_back("./data/demoV.txt");

	//in_path_U.push_back("./data/case/caseU.txt");  // 1: demo
	//in_path_V.push_back("./data/case/caseV.txt");

	in_path_U.push_back("../../data/dblp-author/dblp-authorU.txt");  // 1: dblp  ok  8191, 3090144
	in_path_V.push_back("../../data/dblp-author/dblp-authorV.txt");
	//task_path.push_back("../tipdec/data/dblp-author/dblp-author_insert.txt");

	in_path_U.push_back("../../data/amazon/amazonU.txt");  // 2: amazon  ok  10662, 302
	in_path_V.push_back("../../data/amazon/amazonV.txt");

	in_path_U.push_back("../../data/livejournal/livejournalU.txt");  // 3: livejournal converted  ok  54974, 6691
	in_path_V.push_back("../../data/livejournal/livejournalV.txt");

	in_path_U.push_back("../../data/google/googleU.txt");  // 4: google  ok   10231071, 111
	in_path_V.push_back("../../data/google/googleV.txt");

	in_path_U.push_back("../../data/IMDB/IMDBU.txt");  // 5: IMDB   ok  2280, 98733
	in_path_V.push_back("../../data/IMDB/IMDBV.txt");

	in_path_U.push_back("../../data/twitter/twitterU.txt");  // 6: twitter
	in_path_V.push_back("../../data/twitter/twitterV.txt");

	in_path_U.push_back("../../data/tracker/trackerU.txt");  // 7: tracker
	in_path_V.push_back("../../data/tracker/trackerV.txt");

	in_path_U.push_back("../../data/synthetic/RMAT-20U.txt");  // 8: RMAT-20
	in_path_V.push_back("../../data/synthetic/RMAT-20V.txt");

	in_path_U.push_back("../../data/synthetic/RMAT-21U.txt");  // 9: RMAT-21
	in_path_V.push_back("../../data/synthetic/RMAT-21V.txt");

	in_path_U.push_back("../../data/synthetic/RMAT-22U.txt");  // 10: RMAT-22
	in_path_V.push_back("../../data/synthetic/RMAT-22V.txt");

	in_path_U.push_back("../../data/synthetic/RMAT-23U.txt");  // 11: RMAT-23
	in_path_V.push_back("../../data/synthetic/RMAT-23V.txt");

	in_path_U.push_back("../../data/synthetic/RMAT-24U.txt");  // 12: RMAT-24
	in_path_V.push_back("../../data/synthetic/RMAT-24V.txt");
}

vector<Edge*> top_k_edges_cal(const WorkerParams& params, Graph& graph, vector<Edge*>& global_edge_vec)
{
	size_v k = params.k;
	vector<Edge*> result_edge_vec;
	ResultIdxMap result_idx_map;
	size_v min_idx = 0;
	size_b min_result = 0.0;

	// process k edges in parallel
	//size_v start_edge_idx = size_v(global_edge_vec.size() * 0.1);
	size_v start_edge_idx = 0;

	for (size_v edgeIdx = start_edge_idx; edgeIdx < global_edge_vec.size(); edgeIdx++)
	{

		//if (edgeIdx == 20)  // used to compare the time cost of ego generate and betweenness calculate
		//{
		//	return result_edge_vec;
		//}

		global_cur_idx = edgeIdx;
		if (result_edge_vec.size() == k && min_result >= global_edge_vec[edgeIdx]->up_bound)
		{
			return result_edge_vec;
		}
		else
		{
			Edge* edge = global_edge_vec[edgeIdx];
		//	cout << "Edge (" << edge->src << ", " << edge->des << ") degree: (" << edge->src_degree << ", " << edge->des_degree << ") up-bound: " << edge->up_bound << endl;

			Ego* ego = new Ego(edge);
			ego->set_sets(graph);
			ego->generate_ego(graph);

		//	cout << "generate ego ok !!!" << endl;
			if (global_dynamic_ctr == 2 && result_edge_vec.size() == k)
			{
				edge->ex_edge_num = ego->in_edge_num - edge->src_degree - edge->des_degree + 1;
				edge->set_up_bound();

				if (edge->up_bound <= min_result)
				{
					edge->processed = 1;

					continue;
				}
			}

		//	edge->betweenness = ego->cal_egobetweenness(graph);
			edge->betweenness = ego->cal_egobetweenness_opt(graph);
			
			edge->processed = 1;  // been processed

			global_processed_edge++;
			cout << global_processed_edge << ": Edge (" << edge->src << ", " << edge->des << ") up-bound: " << edge->up_bound << " ego-betweenness: " << edge->betweenness << endl;
			cout << endl;

			if (result_edge_vec.size() < k)
			{
				result_edge_vec.push_back(edge);
				if (result_edge_vec.size() == k)
				{
					min_result = result_edge_vec[0]->betweenness;
					for (size_v resultIdx = 0; resultIdx < result_edge_vec.size(); resultIdx++)
					{
						Edge* tmpedge = result_edge_vec[resultIdx];
						if (tmpedge->betweenness < min_result)
						{
							min_result = tmpedge->betweenness;
							min_idx = resultIdx;
						//	result_idx_map[tmpedge->betweenness] = resultIdx; // used to find the edge with minimal betweenness
						}
					}
				}
			}
			else
			{
				if (edge->betweenness > min_result)
				{
					//size_v edge_pos = result_idx_map[min_result];
					//result_edge_vec[edge_pos] = edge;
					//result_idx_map[edge->betweenness] = edge_pos;

					result_edge_vec[min_idx] = edge;

					min_result = edge->betweenness;
					// update min_result
					for (size_v resultIdx = 0; resultIdx < result_edge_vec.size(); resultIdx++)
					{
						Edge* tmpedge = result_edge_vec[resultIdx];
						if (tmpedge->betweenness < min_result)
						{
							min_result = tmpedge->betweenness;
							min_idx = resultIdx;
						}
					}
				}
			}
		}
	}
	return result_edge_vec;
}

vector<Edge*> top_k_edges_cal_parallel(const WorkerParams& params, Graph& graph, vector<Edge*>& global_edge_vec)
{
	size_v k = params.k;
	vector<Edge*> result_edge_vec;
	ResultIdxMap result_idx_map;
	size_v min_idx = 0;
	size_b min_result = 0.0;
	size_v start_edge_idx = size_v(global_edge_vec.size() * 0.2);
	start_edge_idx = 0;
	// process k edges in parallel
	int num_thread = min(k, global_total_threads);
#pragma omp parallel for num_threads(num_thread)
	for (size_v edgeIdx = start_edge_idx; edgeIdx < start_edge_idx + k; edgeIdx++)
	{
		//global_cur_idx = edgeIdx;
		Edge* edge = global_edge_vec[edgeIdx];

		Ego* ego = new Ego(edge);
		ego->set_sets(graph);
		ego->generate_ego(graph);

		edge->betweenness = ego->cal_egobetweenness_opt(graph);

		edge->processed = 1;  // been processed

	#pragma omp atomic
		global_processed_edge++;
		cout << global_processed_edge << ": Edge (" << edge->src << ", " << edge->des << ") up-bound: " << edge->up_bound << " ego-betweenness: " << edge->betweenness << endl;
		
	#pragma omp critical
		result_edge_vec.push_back(edge);
		if (result_edge_vec.size() == k)
		{
			min_result = result_edge_vec[0]->betweenness;
			for (size_v resultIdx = 0; resultIdx < result_edge_vec.size(); resultIdx++)
			{
				Edge* tmpedge = result_edge_vec[resultIdx];
				if (tmpedge->betweenness < min_result)
				{
					min_result = tmpedge->betweenness;
					min_idx = resultIdx;
					//	result_idx_map[tmpedge->betweenness] = resultIdx; // used to find the edge with minimal betweenness
				}
			}
		}
	}

	size_v terminal_flag = 0;
	size_v stride = k/2;
	size_v start_idx = k + start_edge_idx;
	while (terminal_flag == 0)
	{
		if (start_idx >= global_edge_vec.size())
		{
			break;
		}
	#pragma omp parallel for num_threads(stride)
		for (size_v edgeIdx = start_idx; edgeIdx < global_edge_vec.size(); edgeIdx++)
		{
			//global_cur_idx = edgeIdx;
			if (min_result >= global_edge_vec[edgeIdx]->up_bound)
			{
				terminal_flag = 1;
				//return result_edge_vec;
			}
			else
			{
				Edge* edge = global_edge_vec[edgeIdx];
				//	cout << "Edge (" << edge->src << ", " << edge->des << ") degree: (" << edge->src_degree << ", " << edge->des_degree << ") up-bound: " << edge->up_bound << endl;

				Ego* ego = new Ego(edge);
				ego->set_sets(graph);
				ego->generate_ego(graph);

				//	cout << "generate ego ok !!!" << endl;
				if (global_dynamic_ctr == 2)
				{
					edge->ex_edge_num = ego->in_edge_num - edge->src_degree - edge->des_degree + 1;
					edge->set_up_bound();

					if (edge->up_bound <= min_result)
					{
						edge->processed = 1;
					//	cout << global_processed_edge << ": Edge (" << edge->src << ", " << edge->des << ") up-bound: " << edge->up_bound << " ego-betweenness: " << edge->betweenness << endl;

						continue;
					}
				}

			//	edge->betweenness = ego->cal_egobetweenness(graph);
				edge->betweenness = ego->cal_egobetweenness_opt(graph);

				edge->processed = 1;  // been processed

			#pragma omp atomic
				global_processed_edge++;
				cout << global_processed_edge << ": Edge (" << edge->src << ", " << edge->des << ") up-bound: " << edge->up_bound << " ego-betweenness: " << edge->betweenness << endl;

			#pragma omp critical
				if (edge->betweenness > min_result)
				{
					//size_v edge_pos = result_idx_map[min_result];
					//result_edge_vec[edge_pos] = edge;
					//result_idx_map[edge->betweenness] = edge_pos;

					result_edge_vec[min_idx] = edge;

					min_result = edge->betweenness;
					// update min_result
					for (size_v resultIdx = 0; resultIdx < result_edge_vec.size(); resultIdx++)
					{
						Edge* tmpedge = result_edge_vec[resultIdx];
						if (tmpedge->betweenness < min_result)
						{
							min_result = tmpedge->betweenness;
							min_idx = resultIdx;
						}
					}
				}
			}
		}
		start_idx += stride;  // stride
	}

	return result_edge_vec;
}

int main(int argc, char* argv[])
{
	WorkerParams param;
	set_path(in_path_U, in_path_V, task_path);


#ifdef _WIN32
	//flat_hash_map<size_v, size_v> mp;
	//mp[1] = 1;
	//cout << "mp[1]: " << mp[1] << endl;
#else

#endif



#ifdef _WIN32
	string in_pathU = in_path_U[1];
	string in_pathV = in_path_V[1];
#else
	string in_pathU = in_path_U[atoi(argv[1])];
	string in_pathV = in_path_V[atoi(argv[1])];
#endif

	param.input_pathU = in_pathU;
	param.input_pathV = in_pathV;
#ifdef _WIN32
	param.k = 2;
	global_dynamic_ctr = 0;
	global_priority_ctr = 0;
	global_total_threads = 1;
#else
	global_total_threads = atoi(argv[3]);
	param.k = atoi(argv[4]);
	global_dynamic_ctr = atoi(argv[2]);
	init_timers();
#endif

	Graph graph;
	graph.load_graph(param);

	// openmp test

//#pragma omp parallel for num_threads(4)
//	for (int i = 0; i < 10; i++)
//	{
//		cout << omp_get_thread_num() << " to " << i << endl;
//	}

	 //test 1082368, 5
	//Edge* edge = new Edge;
	//edge->src = 1082368;
	//edge->des = 5;
	//edge->ex_edge_num = 0;
	//edge->edge_id = global_edge_idx;

	//size_v upos = graph.get_vposU(edge->src);
	//edge->src_degree = graph.vertexesU[upos]->degree;

	//size_v vpos = graph.get_vposV(edge->des);
	//edge->des_degree = graph.vertexesV[vpos]->degree;

	//Ego* ego = new Ego(edge);
	//ego->set_sets(graph);
	//ego->generate_ego(graph);

	//edge->ex_edge_num = ego->in_edge_num - edge->src_degree - edge->des_degree + 1;;
	//edge->set_up_bound();

	//ego->cal_egobetweenness_opt(graph);

	//cout << "Edge (" << edge->src << ", " << edge->des << ") " << "up-bound: " << floor(edge->up_bound) << " ego-betweenness: " << ego->betweenness << endl;

	//return 1;

#ifdef _WIN32
#else
	ResetTimer(WORKER_TIMER);
#endif
	for (size_v uidx = 0; uidx < graph.vertexesU.size(); uidx++)
	{
		size_v uid = graph.vertexesU[uidx]->id;
		for (size_v neiidx = 0; neiidx < graph.vertexesU[uidx]->nbs.size(); neiidx++)
		{
			size_v vid = graph.vertexesU[uidx]->nbs[neiidx];
			Edge* edge = new Edge;
			edge->src = uid;
			edge->des = vid;
			edge->ex_edge_num = 0;
			edge->edge_id = global_edge_idx;

			size_v upos = graph.get_vposU(uid);
			edge->src_degree = graph.vertexesU[upos]->degree;

			size_v vpos = graph.get_vposV(vid);
			edge->des_degree = graph.vertexesV[vpos]->degree;

			edge->up_bound = (size_b)(edge->src_degree - 1) * (size_b)(edge->des_degree - 1) + 1.0;
			edge->betweenness = edge->up_bound;

			EdgeIdx edge_idx;
			edge_idx.src = uid;
			edge_idx.des = vid;
			global_edge_map[edge_idx] = edge;

			global_edge_vec.push_back(edge);

			//Ego* ego = new Ego(edge);
			//ego->set_sets(graph);
			//ego->generate_ego(graph);
			//ego->cal_egobetweenness(graph);

		//	cout << "Edge (" << uid << ", " << vid << ") " << "up-bound: " << edge->up_bound << " ego-betweenness: " << edge->betweenness << endl;
			global_edge_idx++;
		}
	}

	sort(global_edge_vec.begin(), global_edge_vec.end(), edge_cmp);

	for (size_v edgeidx = 0; edgeidx < global_edge_vec.size(); edgeidx++)
	{
		global_edge_vec[edgeidx]->global_idx = edgeidx;

		//if (edgeidx < 100)
		//{
		//	cout << "Edge (" << global_edge_vec[edgeidx]->src << ", " << global_edge_vec[edgeidx]->des << ") up-bound: " << global_edge_vec[edgeidx]->up_bound << endl;
		//}
	}

#ifdef _WIN32
#else
	StopTimer(WORKER_TIMER);
	PrintTimer("Pre-processing time: ", WORKER_TIMER);
#endif

#ifdef _WIN32
#else
	ResetTimer(WORKER_TIMER);
#endif
//	vector<Edge*> result_vec = top_k_edges_cal(param, graph, global_edge_vec);  // used to compare the time cost of ego generate and betweenness calculate
	vector<Edge*> result_vec = top_k_edges_cal_parallel(param, graph, global_edge_vec);
	
#ifdef _WIN32
#else
	StopTimer(WORKER_TIMER);
	PrintTimer("Top k search time: ", WORKER_TIMER);
	cout << "Total ego generate time: " << global_egogen_time << endl;
	cout << "Total betweenness calculate time: " << global_betcal_time << endl;
#endif

	cout << "global_processed_edge: " << global_processed_edge << endl;

	sort(result_vec.begin(), result_vec.end(), result_cmp);

	cout << "Result: " << endl;

#ifdef _WIN32
	ofstream caseresult("./case.txt", std::ios::app);
#else

#endif

	for (size_v resultIdx = 0; resultIdx < result_vec.size(); resultIdx++)
	{
		Edge* tmpedge = result_vec[resultIdx];
		cout << tmpedge->global_idx << " Edge (" << tmpedge->src << ", " << tmpedge->des << ") " << "upper bound: " << tmpedge->up_bound << " ego-betweenness: " << tmpedge->betweenness << endl;
	#ifdef _WIN32
		caseresult << tmpedge->up_bound << " " << tmpedge->betweenness << endl;
	#else

	#endif
	}

	//system("pause");
	return 1;
}