#include "Mygraph.h"

MyGraph::MyGraph()
{
}

MyGraph::~MyGraph()
{
	
}

void MyGraph::AddEdge(int from, int to)
{
	this->edges.push_back(Edge(from, to));
	this->out_degree[from]++;
	this->in_degree[to]++;
	int m = edges.size();
	G[from].push_back(m - 1);   //´òÓ¡´ÎÐòÒÀÀµ
	G_2[to].push_back(m - 1);   //´òÓ¡´ÎÐòÒÀÀµ
	G_3[to].push_back(m - 1);   //´òÓ¡´ÎÐò+Ôö²ÄÅö×²ÒÀÀµ
}

void MyGraph::AddEdge_2(int from, int to)
{
	this->edges.push_back(Edge(from, to));
	this->out_degree[from]++;
	this->in_degree[to]++;
	int m = edges.size();
	G_3[to].push_back(m - 1);
	//G_2[to].push_back(m - 1);
}

void MyGraph::UpdateDegree(int node_id, int change)
{
	for (int i = 0; i < G[node_id].size(); i++) {
		if (edge_deleted[G[node_id][i]]) continue;
		int to = this->edges[G[node_id][i]].GetTo();
		this->in_degree[to] += change;
	}
	this->out_degree[node_id] += change * G[node_id].size();
}

void MyGraph::UpdateDegree_2(int node_id, int change)
{
	for (int i = 0; i < G_3[node_id].size(); i++) {
		if (edge_deleted[G_3[node_id][i]]) continue;
		int from = this->edges[G_3[node_id][i]].GetFrom();
		this->out_degree[from] += change;
	}
	this->in_degree[node_id] += change * G_3[node_id].size();
}
