#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <unordered_set>
#include <math.h>

#include "config.h"

using namespace std;

class Vertex
{
private:
	LL *neighbors;
	LL *adj_edges;
	LL degree;
	LL weight;

public:
	Vertex()
	{
		neighbors = NULL;
		adj_edges = NULL;
		degree = 0;
		weight = 0;
	}

	~Vertex()
	{
		free_neighborhood_space();
	}

	void allocate_neighborhood_space(int dgr)
	{
		neighbors = new LL[dgr];
		adj_edges = new LL[dgr];
	}

	void free_neighborhood_space()
	{
		delete[] neighbors;
		delete[] adj_edges;
	}

	void add_neighbor(LL name, LL index)
	{
		neighbors[index] = name;
	}

	void add_adj_edge(LL name, LL index)
	{
		adj_edges[index] = name;
	}

	LL *get_neighbors()
	{
		return neighbors;
	}

	LL *get_adj_edges()
	{
		return adj_edges;
	}

	void set_degree(LL d)
	{
		degree = d;
	}

	LL get_degree()
	{
		return degree;
	}

	LL get_weight()
	{
		return weight;
	}

	void set_weight(LL w)
	{
		weight = w;
	}

	void show_neighbors()
	{
		cout << "neighbors: ";
		for(LL i = 0; i < degree; i++)
		{
			cout << neighbors[i] << '\t';
		}
		cout << endl;
	}
};

class Edge
{
private:
	LL v1, v2;

public:
	Edge(){}
	void set_vertices(LL u1, LL u2)
	{
		v1 = u1;
		v2 = u2;
	}

	void get_vertices(LL &u1, LL &u2)
	{
		u1 = v1;
		u2 = v2;
	}

	~Edge(){}
};

class Graph
{
protected:
	Vertex *vertices;
	Edge	*edges;
	LL v_num;
	LL e_num;
	LL max_degree;
	LL max_vertex_weight;
	LL min_vertex_weight;
	unordered_set<LL> edge_hash_id_set;

private:
	void encode_pairID(LL &pairID, LL n1, LL n2)
	{
		pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
	}

	void decode_pairID(LL pairID, LL &n1, LL &n2)
	{
		LL w = LL((sqrt(double((pairID << 3) + 1)) - 1) / 2);
		LL t = (w * w + w) >> 1;
		n2 = pairID - t;
		n1 = w - n2; 
	}

	void insertEdgeHashIDToSet(LL n1, LL n2)
	{
		LL edge_hash_id;
		LL u, v;
		if(n1 < n2)
		{
			u = LL(n1); v = LL(n2);
		}
		else
		{
			u = LL(n2); v = LL(n1);
		}
		encode_pairID(edge_hash_id, u, v);
		edge_hash_id_set.insert(edge_hash_id);
	}

public:
	Graph(char *filename)
	{
		ifstream infile(filename);
		if(!infile)
		{
			cout << "File " << filename << " cannot be opened" << endl;
			exit(1);
		}

		char line[1024];
		infile.getline(line, 1024);


		while(line[0] != 'p')
			infile.getline(line, 1024);


		char tempstr1[1024], tempstr2[1024];
		sscanf(line, "%s %s %lld %lld", tempstr1, tempstr2, &v_num, &e_num);


		vertices = new Vertex[v_num + 2];
		edges = new Edge[e_num + 2];

		char ch_tmp;

		LL v1, v2;
		LL *v_degree_tmp = new LL[v_num + 2];
		memset(v_degree_tmp, 0, sizeof(LL) * (v_num + 2));

		LL v;
		for(v = 1; v <= v_num; v++)
		{
			LL v_index;
			LL v_weight;
			infile >> ch_tmp >> v_index >> v_weight;
			if(v != v_index)
			{
				cout << "format error occurs when reading v weight" << endl;
				exit(1);
			}
			vertices[v].set_weight(v_weight);
		}

		LL e;
		for(e = 0; e < e_num; e++)
		{
			infile >> ch_tmp >> v1 >> v2;
			edges[e].set_vertices(v1, v2);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
			insertEdgeHashIDToSet(v1, v2);
		}

		for(v = 1; v <= v_num; v++)
		{
			vertices[v].allocate_neighborhood_space(v_degree_tmp[v]);
			vertices[v].set_degree(v_degree_tmp[v]);
		}
		max_degree = v_degree_tmp[1];
		for(LL i = 2; i <= v_num; i++)
		{
			if(v_degree_tmp[i] > max_degree)
				max_degree = v_degree_tmp[i];
		}

		max_vertex_weight = vertices[1].get_weight();
		for(LL i = 2; i <= v_num; i++)
		{
			if(vertices[i].get_weight() > max_vertex_weight)
				max_vertex_weight = vertices[i].get_weight();
		}
		min_vertex_weight = vertices[1].get_weight();
		for(LL i = 2; i <= v_num; i++)
		{
			if(vertices[i].get_weight() < min_vertex_weight)
				min_vertex_weight = vertices[i].get_weight();
		}

		memset(v_degree_tmp, 0, sizeof(LL) * (v_num + 2));
		for(e = 0; e < e_num; e++)
		{
			edges[e].get_vertices(v1, v2);
//cout << v1 << " and " << v2 << endl;
			vertices[v1].add_neighbor(v2, v_degree_tmp[v1]);
			vertices[v2].add_neighbor(v1, v_degree_tmp[v2]);
			vertices[v1].add_adj_edge(e, v_degree_tmp[v1]);
			vertices[v2].add_adj_edge(e, v_degree_tmp[v2]);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
		}

		delete[] v_degree_tmp;
		infile.close();
	}

	Vertex* get_vertices()
	{
		return vertices;
	}

	LL get_max_degree()
	{
		return max_degree;
	}

	bool isConnected(LL n1, LL n2)
	{
		LL edge_hash_id;
		LL u, v;
		if(n1 < n2)
		{
			u = LL(n1); v = LL(n2);
		}
		else
		{
			u = LL(n2); v = LL(n1);
		}
		encode_pairID(edge_hash_id, u, v);
		return edge_hash_id_set.count(edge_hash_id);
	}

	~Graph()
	{
		delete[] vertices;
		delete[] edges;
	}
};

#endif
