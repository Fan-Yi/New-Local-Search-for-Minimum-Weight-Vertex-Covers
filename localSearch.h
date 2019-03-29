#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <cstdlib>
#include  <stdio.h>
#include <float.h>
#include <ctime>// include this header 

#include "graph.h"
#include "myBijection.h"
#include "coverHash.h"
#include "ansUpdate.h"
#include "hugeInt.h"
#include "constants.h"
#include "config.h"
#include "operandSets.h"
#include "weightBuckets.h"


class StochasticLocalSearch : private Graph{

private:
	string solver_name;
	string file_name;

	int seed;

	bool* in_cover;
	LL weighted_cover_size;
	LL cover_weight;

	ConfchangedWeightAgePickSet* ptr_to_no_loss_vertices;

	Bijection* ptr_to_complementary_set_without_fixed;

	Bijection* ptr_to_uncov_edges;
	AnsChangeSet* ptr_to_moved_v;

	LL *score;
	HugeInt *time_stamp;

	bool *confChange;

	HugeInt step;
	LL time_limit;

	LL start, stop;

	bool *best_in_cover;
	LL best_weighted_cover_size;
	LL best_cover_weight;
	LL first_cover_weight;
	double best_cmp_time;
	HugeInt best_solve_step;

#ifdef init_reduction_mode
	bool* is_fixed;
#endif

#ifdef detect_local_optimum_mode
	bool last_step_improved;
#endif
	
	WeightBucket *ptr_to_weight_bucket;
		
#ifdef cover_hash_mode
	CoverHash *ptr_to_hashed_cover;
#endif

#ifdef strategy_analysis_mode
	int ans_update_times;
	int restart_times;
#endif

#ifdef debug_mode
	bool in_local_search_phase;
#endif

public:
	StochasticLocalSearch(char* solv_name, char *fname, int sd, LL cut_off) : Graph(fname)
	{
#ifdef debug_mode
cout << "the problem instance has been constructed" << endl;
#endif
#ifdef debug_mode
		in_local_search_phase = false;
#endif
		solver_name = string(solv_name);
		file_name = string(fname);

		seed = sd;
		step = 0;
		time_limit = cut_off;

		start = clock();
		best_cmp_time = DBL_MAX;

		in_cover = new bool[v_num + 2];
		memset(in_cover, 0, sizeof(bool) * (v_num + 2));
		weighted_cover_size = 0;
		cover_weight = 0;
		first_cover_weight = 0;

		ptr_to_complementary_set_without_fixed = new Bijection(v_num);
		for(LL v = 1; v <= v_num; v++)
		{
			ptr_to_complementary_set_without_fixed->insert_element(v);
		}

#ifdef debug_mode
cout << 1 << endl;
#endif

		ptr_to_no_loss_vertices = new ConfchangedWeightAgePickSet(v_num);

#ifdef debug_mode
cout << 2 << endl;
#endif

		ptr_to_uncov_edges = new Bijection(e_num);
		for(LL e = 0; e < e_num; e++)
		{
			ptr_to_uncov_edges->insert_element(e);
		}

		score = new LL[v_num + 2];
		for(LL v = 1; v <= v_num; v++)
		{
			score[v] = vertices[v].get_degree();
		}

		time_stamp = new HugeInt[v_num + 2];
#ifdef edge_time_stamp_mode
		edge_time_stamp = new HugeInt[e_num + 2];
		for(LL e = 0; e < e_num; e++)
			edge_time_stamp[e] = 0;
#endif
		confChange = new bool[v_num + 2];
		for(LL v = 1; v <= v_num; v++)
		{
#ifdef age_init_i_mode
			time_stamp[v] = -v_num - 1 + v;
#endif
			confChange[v] = 1;
		}
		ptr_to_moved_v = new AnsChangeSet(v_num);

		best_in_cover = new bool[v_num + 2];
		memset(best_in_cover, 0, sizeof(bool) * (v_num + 2));
		best_weighted_cover_size = 0;
		best_cover_weight = 0;
#ifdef detect_local_optimum_mode
		last_step_improved = false;
#endif

#ifdef debug_mode
cout << 3 << endl;
#endif

#ifdef init_reduction_mode
		is_fixed = new bool[v_num + 2];
		memset(is_fixed, 0, sizeof(bool) * (v_num + 2));
		for(LL v = 1; v <= v_num; v++)
		{
			if(!vertices[v].get_degree())
			{
				is_fixed[v] = true;
				ptr_to_complementary_set_without_fixed->delete_element(v);
			}
		}
#endif

#ifdef debug_mode
cout << 4 << endl;
#endif

		ptr_to_weight_bucket = new WeightBucket(v_num, max_vertex_weight, vertices);

#ifdef debug_mode
cout << 5 << endl;
#endif

#ifdef cover_hash_mode
		ptr_to_hashed_cover = new CoverHash(v_num);
#endif

#ifdef debug_mode
cout << 6 << endl;
#endif

		init_solution();

#ifdef debug_mode
cout << 7 << endl;
#endif

#ifdef strategy_analysis_mode
		ans_update_times = 0;
		restart_times = 0;
#endif

#ifdef debug_mode
show_state();
#endif

#ifdef debug_mode
ptr_to_weight_bucket->show_weight_bucket();
cout << "A localSearch object has been constructed" << endl;
#endif
	}

	~StochasticLocalSearch()
	{
		delete[] in_cover;
		delete ptr_to_uncov_edges;
		delete[] time_stamp;
#ifdef edge_time_stamp_mode
		delete[] edge_time_stamp;
#endif
		delete[] score;
		delete[] confChange;
		delete[] best_in_cover;
		delete ptr_to_moved_v;

#ifdef init_reduction_mode
		delete[] is_fixed;
#endif

		delete ptr_to_weight_bucket;

		delete ptr_to_complementary_set_without_fixed;
		delete ptr_to_no_loss_vertices;
		
#ifdef cover_hash_mode
		delete ptr_to_hashed_cover;
#endif
	}

private:

	bool check_solution()
	{
		LL v, e;
		LL verified_weight_sum = 0;
		LL verified_cover_size = 0;
		for(v = 1; v <= v_num; v++)
		{
			if(best_in_cover[v])
			{
				verified_weight_sum += vertices[v].get_weight();
				verified_cover_size++;
			}
		}
		if(verified_weight_sum != best_cover_weight || verified_cover_size != best_weighted_cover_size)
		{
			cout << "verified_weight_sum: " << verified_weight_sum << endl;
			cout << "best_cover_weight: " << best_cover_weight << endl;
			cout << "verified_cover_size: " << verified_cover_size << endl;
			cout << "best_weighted_cover_size: " << best_weighted_cover_size << endl;
			cout << "the total weight is computed incorrectly, or the cover size" << endl;
			return false;
		}

		for(e = 0; e < e_num; e++)
		{
			LL v1, v2;
			edges[e].get_vertices(v1, v2);
			if(!best_in_cover[v1] && !best_in_cover[v2])
			{
				cout << "edge " << e << " is uncovered" << endl;
				cout << "its endpoints are " << v1 << " and " << v2 << endl;
				return false;
			}
		}
		return true;
	}

	bool check_curr_cand_solution()
	{
		LL v, e;
		LL verified_weight_sum = 0;
		LL verified_cover_size = 0;

		for(v = 1; v <= v_num; v++)
		{
			if(in_cover[v])
			{
				verified_weight_sum += vertices[v].get_weight();
				verified_cover_size++;
			}
		}

		if(verified_weight_sum != cover_weight || verified_cover_size != weighted_cover_size)
		{
			cout << "verified_weight_sum: " << verified_weight_sum << endl;
			cout << "cover_weight: " << cover_weight << endl;
			cout << "verified_cover_size: " << verified_cover_size << endl;
			cout << "weighted_cover_size: " << weighted_cover_size << endl;
			cout << "the total weight is computed incorrectly, or the cover size" << endl;
			return false;
		}

		LL verified_uncov_edge_num = 0;
		for(e = 0; e < e_num; e++)
		{
			LL v1, v2;
			edges[e].get_vertices(v1, v2);
			if(!in_cover[v1] && !in_cover[v2])
			{
				//cout << "edge " << e << " is uncovered" << endl;
				//cout << "its endpoints are " << v1 << " and " << v2 << endl;
				if(!ptr_to_uncov_edges->element_in(e))
				{
					cout << e << " is uncovered, but not in uncov_edge set" << endl;
					return false;
				}
				verified_uncov_edge_num++;
			}
		}

		if(verified_uncov_edge_num != ptr_to_uncov_edges->size())
		{
			cout << "verified_uncov_edge_num: " << verified_uncov_edge_num << endl;
			cout << "curr_uncov_edge_num: " << ptr_to_uncov_edges->size() << endl;
			return false;
		}
		return true;
	}


#ifdef init_reduction_mode
	bool degreeTwoReductionTQ(LL &v, LL &nb1, LL &nb2, LL &u, LL &u3, LL &u4)
	//preconditions: score[v]==2
	//return true if TQ reduction works
	{
#ifdef debug_mode
//cout << "considering degree-2 vertex: " << v << endl;
#endif
		LL i;
		bool flag = false;
		LL *nbs = vertices[v].get_neighbors();
		LL dgr = vertices[v].get_degree();
		for(i = 0; i < dgr; i++)
		{
			LL w = nbs[i];
			if(!in_cover[w])
			{
				if(!flag)
				{
					nb1 = w;
					flag = true;
				}
				else
				{
					nb2 = w;
					break;
				}
			}
		}
		u = 0;
		//Tri rule
		if(isConnected(nb1, nb2))
		{
			if(vertices[nb1].get_weight() < vertices[v].get_weight() && vertices[nb2].get_weight() < vertices[v].get_weight()) 
				return true;
		}
		else
		{
			if(vertices[nb1].get_weight() + vertices[nb2].get_weight() < vertices[v].get_weight()) 
				return true;
		}
#ifdef only_tri_rule
		return false;
#endif
		//Qua rule
		LL m, n;
		if(vertices[nb1].get_degree() <= vertices[nb2].get_degree())
		{
			m = nb1; n = nb2;
		}
		else
		{
			m = nb2; n = nb1;
		}
		nbs = vertices[m].get_neighbors();
		dgr = vertices[m].get_degree();
		for(i = 0; i < dgr; i++)
		{
			LL w = nbs[i];
			if(score[w] == 2 && w != v && isConnected(w, n))
			{
				if(isConnected(nb1, nb2))
				{
					if(vertices[nb1].get_weight() < vertices[w].get_weight() + vertices[v].get_weight())
						if(vertices[nb2].get_weight() < vertices[w].get_weight() + vertices[v].get_weight())
						{
							u = w;
							return true;
						}
				}
				else
				{
					if(vertices[nb1].get_weight() + vertices[nb2].get_weight() < vertices[w].get_weight() + vertices[v].get_weight())
					{
						u = w;
						return true;
					}
				}
			}
		}
#ifdef only_tri_and_qua_rule
		return false;
#endif
		// Chain rule
		if(!isConnected(nb1, nb2) && score[nb1] == 2 && score[nb2] == 2)
		{
			LL nb3, nb4;
			for(i = 0; i < vertices[nb1].get_degree(); i++)
			{
				LL w = vertices[nb1].get_neighbors()[i];
				if(!in_cover[w] && w != v)
				{
					nb3 = w;
					break;
				}
			}
			for(i = 0; i < vertices[nb2].get_degree(); i++)
			{
				LL w = vertices[nb2].get_neighbors()[i];
				if(!in_cover[w] && w != v)
				{
					nb4 = w;
					break;
				}
			}
#ifdef debug_mode
cout << "for chain rule, having just obtained " << nb3 << " and " << nb4 << endl;
cout << "3 weights: " << vertices[v].get_weight() << ", " << vertices[nb3].get_weight() << ", " << vertices[nb4].get_weight() << endl;
cout << "2 weights: " << vertices[nb1].get_weight() << ", " << vertices[nb2].get_weight() << endl;
#endif
			if(nb3 != nb4 && vertices[v].get_weight() + vertices[nb3].get_weight() + vertices[nb4].get_weight() < vertices[nb1].get_weight() + vertices[nb2].get_weight())
			{
				u3 = nb3;
				u4 = nb4;
				return true;
			}	
		}
		return false;
	}

	void putVertexIntoCover(const LL v, LL* considered_d1_v_array, LL &size_of_considered_d1_v_array, LL* considered_d2_v_array, LL &size_of_considered_d2_v_array)
	{
#ifdef debug_mode
if(v <= 0 || v > v_num)
{
	cout << "the vertex to be put has an illegal id " << v << endl;
	exit(1);
}

if(in_cover[v])
{
	cout << v << " has already been in the cover, cannot be put in again" << endl;
	exit(1);
}
#endif

#ifdef debug_mode
if(!check_curr_cand_solution())
{
exit(1);
}
#endif

		score[v] = -score[v];
		in_cover[v] = true;//line 
		weighted_cover_size++;//line 
		cover_weight += vertices[v].get_weight();

		ptr_to_weight_bucket->placeInVertexToCover(v, vertices);

		LL *nbs = vertices[v].get_neighbors();
		LL *adj_edgs = vertices[v].get_adj_edges();
		LL dgr = vertices[v].get_degree();

			for(LL j = 0; j < dgr; ++j)//line 
			{

				LL n = nbs[j];//line 

				LL e = adj_edgs[j];
				if(in_cover[n])//line 
				{
					score[n]++; 
					if(!score[n]) // those fixed vertices can never have a score of 0, because they have at least neighbor outside
						ptr_to_no_loss_vertices->insert_element(n);
				}
				else
				{ 
					score[n]--;

					if(score[n] == 2)//
					{
						considered_d2_v_array[size_of_considered_d2_v_array++] = n;//line 
					}
					else if(score[n] == 1)
					{
						considered_d1_v_array[size_of_considered_d1_v_array++] = n;
					}

					ptr_to_uncov_edges->delete_element(e);

				}
			}



	#ifdef cover_hash_mode
		//cout << "put " << v << endl;
		ptr_to_hashed_cover->update_hash_wrt_add(v);
		//cout << "curr hash entry: " << ptr_to_hashed_cover->get_hash_entry() << endl;
	#endif

#ifdef debug_mode
if(!check_curr_cand_solution())
{
exit(1);
}
#endif
	}

#endif

	void init_solution()
	{
		LL i, j;
		LL v, e;
		LL* nbs;
		LL dgr;

#ifdef init_reduction_mode

		bool optimality_still_guaranteed = true;

		LL* considered_d1_v_array = new LL[v_num + 2];
		LL size_of_considered_d1_v_array = 0;
		LL* considered_d2_v_array = new LL[v_num + 2];
		LL size_of_considered_d2_v_array = 0;
		
		for(v = 1; v <= v_num; v++)
		{
			if(vertices[v].get_degree() == 1)
				considered_d1_v_array[size_of_considered_d1_v_array++] = v;
			else if(vertices[v].get_degree() == 2)
				considered_d2_v_array[size_of_considered_d2_v_array++] = v;
		}
#endif

		while(ptr_to_uncov_edges->size())
		{

#ifdef init_reduction_mode
			LL randPtr;
			LL v_to_add;//init to be an illegal vertex

			while(size_of_considered_d1_v_array)// exists 1-gain vertices, line 8
			{

				LL gain_one_v = considered_d1_v_array[--size_of_considered_d1_v_array];//line 10
				if(score[gain_one_v] != 1)
					continue;
				LL i;
				nbs = vertices[gain_one_v].get_neighbors();
				dgr = vertices[gain_one_v].get_degree();

				for(i = 0; i < dgr; i++)//line 10
				{
					if(!in_cover[nbs[i]])//search for the unique neighbor
						break;
				}

				v_to_add = nbs[i];//line 10

				if(vertices[v_to_add].get_weight() < vertices[gain_one_v].get_weight())
				{				
					putVertexIntoCover(v_to_add, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
					if(optimality_still_guaranteed)
					{
						is_fixed[v_to_add] = is_fixed[gain_one_v] = true;
						
						ptr_to_complementary_set_without_fixed->delete_element(gain_one_v);

#ifdef debug_mode
cout << "having applied degree-1 rule, put " << v_to_add << " inside, and keep " << gain_one_v << " outside" << endl;
//getchar();
#endif
					}

					ptr_to_complementary_set_without_fixed->delete_element(v_to_add);
	//cout << "added " << v_to_add << endl;
				}

			}


			while(size_of_considered_d2_v_array)
			{
				LL u = considered_d2_v_array[--size_of_considered_d2_v_array];
				LL nb1, nb2;
				LL w = 0, w3 = 0, w4 = 0;
				if(score[u] == 2 && degreeTwoReductionTQ(u, nb1, nb2, w, w3, w4))
				{
					if(w3 == 0) // and w4 == 0
					{
						putVertexIntoCover(nb1, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
						putVertexIntoCover(nb2, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
						if(optimality_still_guaranteed)
						{
							is_fixed[u] = is_fixed[nb1] = is_fixed[nb2] = true;
							
							ptr_to_complementary_set_without_fixed->delete_element(u);
							if(w != 0) // for Qua rule
							{
								is_fixed[w] = true;
								ptr_to_complementary_set_without_fixed->delete_element(w);
							}
#ifdef debug_mode
	cout << "having applied degree-2 TQ rule, put " << nb1 << " and " << nb2 << " inside, and keep " << u << " outside, ";
							if(w != 0)
							{
								cout << "also keep " << w << " outside.";
							}
	cout << endl;
	//getchar();
#endif						
						}

						ptr_to_complementary_set_without_fixed->delete_element(nb1);
						ptr_to_complementary_set_without_fixed->delete_element(nb2);
				
		//cout << "added " << nb1 << " and " << nb2 << endl;
					}
					else
					{
#ifdef debug_mode
cout << u << ", " << w3 << ", " << w4 << endl;
#endif
						putVertexIntoCover(u, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);						
						putVertexIntoCover(w3, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
						putVertexIntoCover(w4, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
#ifdef debug_mode
cout << "having just add 3 vertices into cover" << endl;
#endif
						if(optimality_still_guaranteed)
						{
							is_fixed[u] = is_fixed[w3] = is_fixed[w4] = true;
							is_fixed[nb1] = is_fixed[nb2] = true;
							ptr_to_complementary_set_without_fixed->delete_element(nb1);
							ptr_to_complementary_set_without_fixed->delete_element(nb2);
#ifdef debug_mode
	cout << "having applied degree-2 TQ rule, put " << w3 << " and " << u << " and " << w4 << " inside, ";
	cout << "also keep " << nb1 << " and " << nb2 << " outside" << endl;
	//getchar();
#endif	
						}
						ptr_to_complementary_set_without_fixed->delete_element(u);
						ptr_to_complementary_set_without_fixed->delete_element(w3);
						ptr_to_complementary_set_without_fixed->delete_element(w4);						
					}
				}
				
			}


			if(size_of_considered_d1_v_array)// and size_of_considered_d2_v_array != 0
				continue;
			if(!ptr_to_uncov_edges->size()) break;

			optimality_still_guaranteed = false;

#endif

				LL rand_uncov_e = ptr_to_uncov_edges->rand_element();
				LL w1, w2;
				LL add_v;
				edges[rand_uncov_e].get_vertices(w1, w2);
				if(vertices[w1].get_degree() > vertices[w2].get_degree())
					add_v = w1;
				else
					add_v = w2;
#ifdef init_reduction_mode
				putVertexIntoCover(add_v, considered_d1_v_array, size_of_considered_d1_v_array, considered_d2_v_array, size_of_considered_d2_v_array);
				ptr_to_complementary_set_without_fixed->delete_element(add_v);
#else
				add(add_v);
#endif

	//cout << "added " << v_to_add << endl;
		}//end while(ptr_to_uncov_edges->size() != 0)

#ifdef debug_mode
LL fixed_vertex_num = 0;
for(LL v = 1; v <= v_num; v++)	
{
	if(is_fixed[v]) fixed_vertex_num++;
}
cout << "fixed_vertex_num: " << fixed_vertex_num << endl;
getchar();
#endif

#ifdef debug_mode
	#ifdef init_reduction_mode
if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
{
	cout << "right after reduction, the weight bucket is incorrect" << endl;
	exit(1);
}
	#endif
#endif

#ifdef init_reduction_mode
	#ifdef debug_mode
cout << "before filtering out vertices" << endl;
ptr_to_weight_bucket->show_weight_bucket();
	#endif
		ptr_to_weight_bucket->filter_out_dead_vertices(is_fixed, v_num);
		delete[] considered_d2_v_array;
		delete[] considered_d1_v_array;
	#ifdef debug_mode
//cout << "after filtering out vertices" << endl;
//ptr_to_weight_bucket->show_weight_bucket();
	#endif
#endif

#ifdef debug_mode
cout << "about to remove redundant vertices" << endl;
#endif

		for(v = 1; v <= v_num; v++)
		{
			if(in_cover[v] && !score[v])// in our program, if score[v]==0, then v must not be fixed
				remove(v);
		}

#ifdef debug_mode
cout << "having just removed all redundant vertices" << endl;
	#ifdef init_reduction_mode
if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
{
	cout << "at the end of the construction function of localSearch, the weight bucket is incorrect" << endl;
	exit(1);
}
	#endif
#endif
		
		update_best_weighted_cover();

#ifdef ans_update_op_mode
		ptr_to_moved_v->efficiently_update_best_weighted_covering_vertices();
#endif

		first_cover_weight = cover_weight;

		stop = clock();
		best_cmp_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;
		best_solve_step = step;

#ifndef only_test_init_mode
	#ifdef init_reduction_mode
	if(optimality_still_guaranteed)
	#else
	if(false)
	#endif
#endif
	{
		if(check_solution() == 1)
			{
				cout << "o " << best_cover_weight << endl;
				cout << "c size " << best_weighted_cover_size << endl;
				cout << "c searchSteps " << best_solve_step <<endl;
				cout << "c solveTime " << best_cmp_time << endl;
			}
			else
			{
				cout << "the solution is wrong." << endl;
			}
			exit(0);
	}


	#ifdef individual_analysis_on_init_sls_mode
			stop = clock();
			init_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;
			init_time = round(init_time * 100) / 100.0;
	#endif
	}



	void update_best_weighted_cover()
	{
		for(LL v = 1; v <= v_num; v++)
		{
			best_in_cover[v] = in_cover[v];
		}
		best_weighted_cover_size = weighted_cover_size;
		best_cover_weight = cover_weight;
#ifdef strategy_analysis_mode
ans_update_times++;
#endif
	}

	void add(const LL v)
	{
#ifdef debug_mode
cout << "step: " << step << endl;
cout << "add " << v << ", its gain is " << score[v] << endl;
if(in_local_search_phase)
{
	if(score[v] != 0)
	{
		cout << "about to add a vertex whose gain is NOT 0" << endl;
		exit(1);
	}
}
#endif
		if(confChange[v])
		{
			confChange[v] = false;
		}
		in_cover[v] = true;
		ptr_to_complementary_set_without_fixed->delete_element(v);
#ifndef init_reduction_mode
		if(score[v] == 0) // in the init_solution() function, we use add(const LL) to add a vertex whose score is bigger than 0
#endif
			ptr_to_no_loss_vertices->insert_element(v);
		weighted_cover_size++;
		cover_weight += vertices[v].get_weight();
		//score[v] = -score[v];
#ifdef debug_mode
//cout << "before adding..." << endl;
//ptr_to_weight_bucket->show_weight_bucket();
/*
#ifdef init_reduction_mode
if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
{
	cout << "about to be added into cover, but not in proper out-bucket" << endl;
	exit(1);
}
#endif

if(!ptr_to_weight_bucket->is_in_vertex_in_bucket_of_weight(v, vertices[v].get_weight()))
{
	cout << "about to be added into cover, but not in proper out-bucket" << endl;
	exit(1);
}
*/
#endif

		ptr_to_weight_bucket->placeInVertexToCover(v, vertices);	

#ifdef debug_mode
//cout << "after adding..." << endl;
//ptr_to_weight_bucket->show_weight_bucket();

#ifdef init_reduction_mode
if(!ptr_to_weight_bucket->check_weight_bucket(vertices, v_num, is_fixed))
{
	cout << "about to be added into cover, but not in proper out-bucket" << endl;
	exit(1);
}
#endif

if(!ptr_to_weight_bucket->is_in_vertex_in_bucket_of_weight(v, vertices[v].get_weight()))
{
	exit(1);
}
#endif	
		LL* nbs = vertices[v].get_neighbors();
		LL* adj_edgs = vertices[v].get_adj_edges();
		LL dgr = vertices[v].get_degree();
		for(LL i = 0; i < dgr; i++)
		{
			LL n = nbs[i];
			LL e = adj_edgs[i];
			if(in_cover[n])
			{
				score[n]++;
				if(!score[n])
					ptr_to_no_loss_vertices->insert_element(n);
				if(!confChange[n])
				{
					confChange[n] = true;
				}
			}
			else
			{
				score[n]--;
				ptr_to_uncov_edges->delete_element(e);
#ifdef edge_time_stamp_mode
				edge_time_stamp[e] = step;
#endif
			}
		}
		time_stamp[v] = step;
#ifdef cover_hash_mode
		ptr_to_hashed_cover->update_hash_wrt_add(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif
		step++;
#ifdef debug_mode
if(!check_curr_cand_solution())
{
	exit(1);
}
#endif
	}

	void remove(const LL v)
	{
#ifdef debug_mode
cout << "step: " << step << endl;
cout << "remove " << v << ", the loss is " << -score[v] << endl;
if(score[v] != 0)
{
	cout << "about to remove a vertex whose loss is NOT 0" << endl;
	exit(1);
}
#endif
		if(!confChange[v])
		{
			confChange[v] = true;
		}

		in_cover[v] = false;
#ifdef debug_mode
//cout << 1 << endl;
//ptr_to_no_loss_vertices->show_elements();
#endif
		ptr_to_no_loss_vertices->delete_element(v);// when we call remove(const LL), the score of v is always 0
#ifdef debug_mode
//cout << 2 << endl;
#endif
		ptr_to_complementary_set_without_fixed->insert_element(v);
		weighted_cover_size--;
		cover_weight -= vertices[v].get_weight();
		ptr_to_weight_bucket->placeOutVertexFromCover(v, vertices);
		LL* nbs = vertices[v].get_neighbors();
		LL* adj_edgs = vertices[v].get_adj_edges();
		LL dgr = vertices[v].get_degree();
		for(LL i = 0; i < dgr; i++)
		{
			LL n = nbs[i];
			LL e = adj_edgs[i];
			if(in_cover[n])
			{
				if(!score[n])
					ptr_to_no_loss_vertices->delete_element(n);
				score[n]--;
				if(!confChange[n])
				{
					confChange[n] = true;
				}
			}
			else
			{
				score[n]++;
				ptr_to_uncov_edges->insert_element(e);
#ifdef edge_time_stamp_mode
				edge_time_stamp[e] = step;
#endif
			}
		}
		time_stamp[v] = step;
#ifdef cover_hash_mode
		ptr_to_hashed_cover->update_hash_wrt_remove(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif
		step++;

#ifdef debug_mode
if(!check_curr_cand_solution())
{
	exit(1);
}
#endif
	}

	LL randVertexInCover()
	{
		LL v;
		do{
				v = rand() % v_num + 1;
			}while(!in_cover[v]);
		return v;
	}

	LL rand_cov_edge()
	{
		LL e;
		do{
			e = rand() % e_num;
		}while(ptr_to_uncov_edges->element_in(e));
		return e;
	}

	LL rand_vertex_outside_cover()
	{
		LL v;
		do{
			v = rand() % v_num + 1;
		}while(in_cover[v]);
		return v;
	}

	LL rand_covering_vertex(LL e)// e must be a covered edge
	{
		LL v1, v2;
		edges[e].get_vertices(v1, v2);
		if(!in_cover[v1]) return v2;
		else if(!in_cover[v2]) return v1;
		else if(rand() % 2) return v2;
		else return v1;
	}

	void clear()
	{
		while(ptr_to_complementary_set_without_fixed->size())
		{
			LL v = ptr_to_complementary_set_without_fixed->at(ptr_to_complementary_set_without_fixed->begin());
			add(v);
		}
	
		while(ptr_to_no_loss_vertices->size())
		{
			LL v = ptr_to_no_loss_vertices->rand_element();
			remove(v);
		}

#ifdef detect_local_optimum_mode
		last_step_improved = true;
#endif

#ifdef scenario_hash_mode
		ptr_to_hashed_cover->forget_all_visits();
#endif
	}

	void local_move()
	{
		LL i;
		LL u, v;

		LL best_remove_v;
		
		do{
			best_remove_v = ptr_to_no_loss_vertices->best_element(confChange, vertices, time_stamp);
			if(best_remove_v == 0)
				break;
			remove(best_remove_v);
#ifdef detect_local_optimum_mode
			last_step_improved = true;
#endif
		}while(true);


#ifdef detect_local_optimum_mode
		if(last_step_improved) // local optimum
		{
			if(cover_weight < best_cover_weight)
			{
#ifdef ans_update_op_mode
				ptr_to_moved_v->efficiently_update_best_weighted_covering_vertices();
				best_weighted_cover_size = weighted_cover_size;
				best_cover_weight = cover_weight;		
#else
				update_best_weighted_cover();
#endif

#ifdef strategy_analysis_mode
ans_update_times++;
#endif

				stop = clock();
				best_cmp_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;		
				best_solve_step = step;
			}
#ifdef cover_hash_mode
			if(ptr_to_hashed_cover->curr_hash_entry_visited())
			{
				clear();
#ifdef strategy_analysis_mode
restart_times++;
#endif
				return;
			}
			else
			{
				ptr_to_hashed_cover->mark_hash_entry();
			}
#endif
		}
#endif

		LL best_add_v;

		best_add_v = ptr_to_weight_bucket->vertex_outside_cover_of_smallest_weight_then_time_stamp(time_stamp);
		if(best_add_v != 0)
		{
			add(best_add_v);
		}
		else
		{
			clear();
#ifdef strategy_analysis_mode
restart_times++;
#endif
			return;
		}
#ifdef detect_local_optimum_mode
		last_step_improved = false;
#endif

	}

public:

	void cover_sls()
	{
#ifdef debug_mode
cout << "beginning local move" << endl;
in_local_search_phase = true;
#endif

		while(true)
		{
				LL step_bound = 0;
				if(LL(step / TRY_STEP) > step_bound)
				{
					step_bound++;

					stop = clock();
					double elap_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;
					if(LL(elap_time) >= LL(time_limit))
					{
#ifdef ans_update_op_mode
						ptr_to_moved_v->compute_final_answer(v_num, in_cover, best_in_cover);
#endif
#ifdef debug_mode
cout << "best_cover_weight: " << best_cover_weight << endl;
cout << "best_weighted_cover_size: " << best_weighted_cover_size << endl;
#endif
						return;
					}
				}
				local_move();
//show_state();
		}
	}

	void show_results()
	{

		stop = clock();
		double elap_time = (stop - start) / double(CLOCKS_PER_SEC) * 1000;
		if(check_solution())
		{
			cout << "o " << best_cover_weight << endl;
			cout << "c size " << best_weighted_cover_size << endl;
			cout << "c solveTime " << best_cmp_time << endl;
			cout << "c searchSteps " << best_solve_step << endl;
			cout << "c stepSpeed(/ms) " << step / elap_time << endl;
			cout << "c solver_name " << solver_name << endl;
			cout << "c instance_name " << file_name << endl;
			cout << "c seed " << seed << endl;
			cout << "c user_cutoff " << time_limit << endl;
			cout << "c actual_cutoff " << elap_time << endl;
			cout << "c v_num " << v_num << endl;
			cout << "c e_num " << e_num << endl;
			cout << "c avg_degree " << round(double(2 * e_num) / double(v_num) * 100) / 100.0 << endl;
			cout << "c density " << round(double(2 * e_num) / double(v_num * (v_num - 1)) * 1000000) / 1000000.0 << endl;
			cout << "c total_step " << step << endl;
			cout << "c first_cover_weight " << first_cover_weight << endl;
#ifdef strategy_analysis_mode
cout << "c restart_times " << restart_times << endl;
cout << "the total number of steps: " << step << endl;
cout << "answer update times: " << ans_update_times << endl;
#endif
		}
		else
		{
			cout << "sorry, something is wrong" << endl;
		}
	}	

	void show_state()
	{
		LL v, e;
/*
		for(e = 0; e < e_num; e++)
		{
			cout << edges[e].get_v1() << " " << edges[e].get_v2() << endl;
		}

		for(v = 1; v <= v_num; v++)
		{
			cout << "vertex " << v << ": " << endl;
			vertices[v].show_neighbors();
		}
*/
		LL i, j;
		cout << "step: " << step << endl;
		cout << "the cover: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(in_cover[v])
				cout << v << '\t';
		}
		cout << endl;
#ifdef cover_hash_mode
//cout << "curr_hash_entry: " << ptr_to_hashed_cover->get_hash_entry() << endl;
#endif
/*
		cout << "loss: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(in_cover[v])
				cout << -score[v] << '\t';
		}
		cout << endl;

		cout << "uncov edges: " << endl;
		for(i = ptr_to_uncov_edges->begin(); i < ptr_to_uncov_edges->end(); i++)
		{
			int e = ptr_to_uncov_edges->at(i);
			cout << e;
			int v1, v2;
			edges[e].get_vertices(v1, v2);
			cout << ", " << "endpoints: " << v1 << " and " << v2 << endl;
		}
		cout << endl;

		cout << "cover weight: " << cover_weight << endl;
		cout << "weighted cover size: " << weighted_cover_size << endl;


		cout << endl;
*/		

cout << "**************************************" << endl;
//getchar();
	}
};
#endif
