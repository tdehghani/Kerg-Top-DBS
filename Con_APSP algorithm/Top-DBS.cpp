/*
*   Top-DBS algorithm for Beta-sheet Topology Prediction Problem
*   Copyright (C) 2018 ،Toktam Dehghani
* 	Knowledge Engineering Research Group (KERG), Department of  Computer Engineering,
*   Faculty of Engineering, Ferdowsi University of Mashhad, Iran. https://kerg.um.ac.ir/
*   Toktam Dehghani, email: dehghani.toktam@mail.um.ac.ir
*
* main.c - This file is part of Top-DBS
*
* Top-DBS is a free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundation; either version 3 of the License,
* or (at your option) any later version.
*
* Top-DBS is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
*/
#define _CRT_SECURE_NO_DEPRECATE
#include<conio.h>
#include <math.h>
#include <sstream>
#include "strand.h"
#include <algorithm> 
#include "find_max_beta_topology.h"
using namespace std;
double infinite = 99999;
vector<double> vertices_score;
vector<vector<int>> neighbors_list;
vector<int> indexes;
int last_infinite = 0;
vector<double> maximum_possibility;
vector<vector<bool>> bool_PATH_matrix;//// Matrix for storing  path in bool(0,1) between pairs in each step
struct max_clique_struct
{
	double score;
	vector<int> vertices;
};
vector<max_clique_struct> max_clique;
bool my_function(int i1, int i2)
{
	return vertices_score[i1] < vertices_score[i2];
}
bool my_function2(max_clique_struct clique1, max_clique_struct clique2)
{
	return clique1.score < clique2.score;
}
void func_P_or_AP_clique(int x, int y,int num_row, vector<vector<int>> &pair_strands, vector<char> &interaction_type)
{
	vector<int> new_pair;
	char new_interaction;
	if (x < (num_row / 2) && y < (num_row / 2))
	{
		new_pair.push_back(x);
		new_pair.push_back(y);
		new_interaction = 'P';
	}
	if (x >= (num_row / 2) && y >= (num_row / 2))
	{
		new_pair.push_back(x - (num_row / 2));
		new_pair.push_back(y - (num_row / 2));
		new_interaction = 'P';
	}
	if (x < (num_row / 2) && y >= (num_row / 2))
	{
		new_pair.push_back(x);
		new_pair.push_back(y - (num_row / 2));
		new_interaction = 'A';
	}
	if (x >= (num_row / 2) && y< (num_row / 2))
	{
		new_pair.push_back(x - (num_row / 2));
		new_pair.push_back(y);
		new_interaction = 'A';
	}
	interaction_type.push_back(new_interaction);
	pair_strands.push_back(new_pair);
}
void func_print_sheet_clique(vector<int> PATH,int max_row, int max_col, int step, int num_row, FILE* ofp_conformation)// , vector<vector<double>> A_matrix, vector<vector<int>> P_matrix, vector<vector<vector<int>>> Path_matrix, int num_row)// clique algorithm find maximum sheets )
{
	vector<vector<int>> pair_strands;
	vector<char> interaction_type;

	for (int k = 0; k < PATH.size()-1; k++)
		func_P_or_AP_clique(PATH[k], PATH[k+1], num_row, pair_strands, interaction_type);
	printf("  		<The predicted strand pairs and strand alignments by Top - DBS:\n");
	
	for (int i = 0; i < interaction_type.size(); i++)
		{
		    printf("  		%d--", pair_strands[i][0] + 1);
			fprintf(ofp_conformation, "%d--", pair_strands[i][0] + 1);
			printf("%d", pair_strands[i][1] + 1);
			fprintf(ofp_conformation, "%d:", pair_strands[i][1] + 1);
			printf(": %c\n", interaction_type[i]);
			fprintf(ofp_conformation, "%c\n", interaction_type[i]);
			if (i != interaction_type.size() - 1)
				fprintf(ofp_conformation, "\n");
		}
}
//clique_info make_sheets_graph(int k_parameter, vector<vector<vector<double>>> A_matrix, vector<vector<vector<vector<int>>>> PATH_matrix, int num_row)

void make_sheets_graph(int k_parameter,vector<vector<vector<double>>> A_matrix, vector<vector<vector<vector<int>>>> PATH_matrix,int num_row)
{
	indexes.clear();
	vertices_score.clear();
	neighbors_list.clear();
	vertices_score.resize(num_row * num_row*(num_row+1));// score of each sheet
	indexes.resize(num_row * num_row*(num_row+1));
	bool_PATH_matrix.resize(num_row * num_row*(num_row + 1),vector<bool>(num_row/2));
	
	int m = 0;
	double density = 0;
	int num_edge = 0, num_vertice = 0;
	int treewidth = 0;
	// read vertices scores
	int h = 0;
	//cout << "create nodes of graph \n";
	for (int i = 0; i < num_row; i++)
	{
		for (int j = 0; j < num_row; j++)
		{
			for (int k = 0; k < num_row+1; k++)
			{
				vertices_score[h] = A_matrix[i][j][k];
				if (A_matrix[i][j][k] != (-1 * infinite))
					num_vertice++;					
					indexes[h] = h;
					for (int w = 0; w < PATH_matrix[i][j][k].size(); w++)
					{
						if (PATH_matrix[i][j][k][w] != -1)
							bool_PATH_matrix[h][(PATH_matrix[i][j][k][w]) % (num_row / 2)] = 1;
					}
					h++;
			}
		}
	}
	//cout << "sort\n";
	//sort index of vertices in descending order
	sort(indexes.begin(), indexes.end(), my_function);
	reverse(indexes.begin(), indexes.end());
	//cout << "create neighbor list\n";
	int temp = 0;
	vector<int> new_row;	
	//read nighbor_list
	for (int i = 0; i < vertices_score.size(); i++)
	{
		if (vertices_score[i] != (infinite*-1))
		{
			temp = 0;
			new_row.clear();
			for (int j = 0; j < vertices_score.size(); j++)
			{
				if(i!=j)
				{ 
					if (vertices_score[j] != (infinite*-1))
					{
						int conflict = 0;
						for (int w = 0; (w < num_row / 2) && (conflict == 0); w++)
							if (bool_PATH_matrix[i][w] == true && bool_PATH_matrix[j][w] == true)
								conflict = 1;
						if (conflict == 0)
						{
							new_row.push_back(j);
							num_edge++;
							temp++;
						}
					}
				}
			}
		}
		if (temp > treewidth)
			treewidth = temp;
		if (new_row.size() == 0)
			new_row.push_back(-1);
		neighbors_list.push_back(new_row);
	}

	//sort neighbor list
	int max_num_neighbor = num_row*(k_parameter);
	for (int i = 0; i < vertices_score.size(); i++)
	{
		if (neighbors_list[i][0] != (-1))
		{
			vector<int> neighbors_indexes;
			neighbors_indexes = neighbors_list[i];
			sort(neighbors_indexes.begin(), neighbors_indexes.end(), my_function);
			reverse(neighbors_indexes.begin(), neighbors_indexes.end());
			neighbors_list[i] = neighbors_indexes;			
			if(neighbors_list[i].size()>max_num_neighbor)
				neighbors_list[i]._Pop_back_n(neighbors_list[i].size()- max_num_neighbor);
		
		}
	}
	int list_max_num_neighbor = 0;
	int sum_num_edge = 0;
	int sum_num_node = 0;
	int max_neighbor = 0;
	for (int i = 0; i <vertices_score.size(); i++)
	{
		if (vertices_score[i] == (infinite*-1))
		{
			maximum_possibility.push_back(vertices_score[i]);
		}

		if (vertices_score[i] != (infinite*-1))
		{
			sum_num_node++;
			if (list_max_num_neighbor < neighbors_list[i].size())
				list_max_num_neighbor = neighbors_list[i].size();
			double count = 0;
			if (max_neighbor < neighbors_list[i].size())
				max_neighbor = neighbors_list[i].size();

			for (int j = 0; (j < neighbors_list[i].size() && (neighbors_list[i][0] != (-1))); j++)
			{
				sum_num_edge++;
				count = count + vertices_score[neighbors_list[i][j]];
			}
			maximum_possibility.push_back(count);
		}
	}
	//clique_info my_clique;

	density = num_edge ;//include i to j and j to i , no need to multiple by two
	density = density / (vertices_score.size()*(vertices_score.size()- 1));// density= 2*E/(v(v-1)) in undirect graphs
									
    /*my_clique.number_edge = num_edge / 2;
	my_clique.number_node = num_vertice;
	my_clique.tree_width = treewidth;
	my_clique.density = density;
	my_clique.prune_number_edge = sum_num_edge / 2;
	my_clique.prune_number_node = sum_num_node;
	my_clique.prune_tree_width = max_neighbor;*/
	//return( my_clique);
}
void initilized_global_max_clique()
{
	max_clique.clear();
	max_clique_struct new_clique;
	new_clique.score = 0;
	new_clique.vertices.push_back(-1);
	max_clique.push_back(new_clique);
}
void clique(vector<int> u, double cur_weight, vector<int> cur_vertice,  int num_of_sheet)
{
	if (u.size() == 0)
	{
		if ((cur_weight >= max_clique[max_clique.size()-1].score) &&(cur_vertice.size()==num_of_sheet))//??????? need =
		{
				max_clique_struct new_clique;
				new_clique.score = cur_weight;
				new_clique.vertices = cur_vertice;
				max_clique.push_back(new_clique);			
			
		}
	}
	else
	{
		
		int fin = 0;
		while ((u.size()> 0) && (fin == 0))
		{
			double sum = 0;
			for (int j = 0; j < u.size(); j++)
			{
				sum += vertices_score[u[j]];
			}
			if (((cur_weight + sum) >= max_clique[max_clique.size() - 1].score))
			{
				vector<int> new_u;
				double new_cur_weight;
				vector<int> new_cur_vertice;
				new_cur_vertice = cur_vertice;
				new_cur_vertice.push_back(u[0]);
				new_cur_weight = cur_weight;
				new_cur_weight += vertices_score[u[0]];
				for (int i = 0; (i < neighbors_list[u[0]].size() && neighbors_list[u[0]][i] != -1); i++)
				{
					int temp = neighbors_list[u[0]][i];
					int done = 1;
					if (new_cur_weight + vertices_score[temp] + maximum_possibility[temp] >= max_clique[max_clique.size() - 1].score)
					{
						for (int j = 1; (j < u.size() && (done == 1)); j++)
						{
							if (temp == u[j])
							{
								new_u.push_back(temp);
								done = 0;
							}
						}
					}
				}
				clique(new_u, new_cur_weight, new_cur_vertice,num_of_sheet);
				u.erase(u.begin());
			}
			else
			{
				fin = 1;
			}
		}
	}//else
}
void max_weight_clique( int num_of_sheet)
{
	double cur_weight;
	vector<int> cur_vertice;
	int cur_neighbor;
	vector<int> u;
	//max_clique.clear();
	for (int i = 0; i < (vertices_score.size()); i++)
	{
		cur_vertice.clear();
	//	cur_vertice.push_back(indexes[i]);//sort
		cur_vertice.push_back(i);//no sort
		cur_weight = vertices_score[cur_vertice[0]];
		if (cur_weight != (infinite*-1))
		{
			if (cur_weight + maximum_possibility[cur_vertice[0]] >= max_clique[max_clique.size() - 1].score)
			{
				u.clear();
				for (int j = 0; (j < neighbors_list[cur_vertice[0]].size() && neighbors_list[cur_vertice[0]][j] != (-1)); j++)
				{
						cur_neighbor = neighbors_list[cur_vertice[0]][j];
						if(cur_neighbor>i)
						{
							if (cur_weight + (vertices_score[cur_neighbor] + maximum_possibility[cur_neighbor]) >= max_clique[max_clique.size() - 1].score)
							{
								u.push_back(cur_neighbor);
							}
						}
					
				}
				clique(u, cur_weight, cur_vertice,num_of_sheet);
			}
		}
	}
}


clique_info clique_algorithm(vector<vector<vector<double>>> SCORE_matrix, vector<vector<vector<vector<int>>>> PATH_matrix,  int num_row,int num_of_sheet, char * clique_fold_path_str)
{
	int k_parameter = 4;// num_row / 2;
	max_clique_struct new_clique;
	clock_t begin, end;
	double time_spent = 0;
	vector<max_clique_struct> global_max_clique;
	clique_info my_clique;	
	//my_clique=make_sheets_graph(k_parameter,SCORE_matrix,PATH_matrix,num_row);
	make_sheets_graph(k_parameter, SCORE_matrix, PATH_matrix, num_row);
	begin = clock();
	for (int i = 0; i < num_of_sheet; i++)
	{
		//cout << "clique" << num_of_sheet << "\n";
		initilized_global_max_clique();
		max_weight_clique(i + 1);
		for (int j = 0; j < max_clique.size(); j++)
			if(max_clique[j].score!=0)
			   global_max_clique.push_back(max_clique[j]);
	}
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	my_clique.run_time = time_spent;
	//cout << " Time:" << time_spent<<"\n";
	sort(global_max_clique.begin(), global_max_clique.end(), my_function2);
	reverse(global_max_clique.begin(), global_max_clique.end());
	FILE* ofp_conformation = fopen(clique_fold_path_str, "w");
	if (ofp_conformation == NULL)
	{
		//cout << "\n error: output file of clique \n";
	}
	else
	{
		fprintf(ofp_conformation, "Predicted beta-topology:\n");
			for (int w = 0; (w < 1); w++)
			{
				//fprintf(ofp_conformation, "Independent Set ID:\n%d\n\n", w+1);
				new_clique = global_max_clique[w];
				for (int z = 0; (z < new_clique.vertices.size()); z++)
				{
					int row = (new_clique.vertices[z] / (num_row*(num_row + 1)));
					int col = (new_clique.vertices[z] - row*(num_row*(num_row + 1))) / (num_row + 1);
					int step = (new_clique.vertices[z] - row*(num_row*(num_row + 1)) - col*(num_row + 1));
					func_print_sheet_clique(PATH_matrix[row][col][step], row, col, step, num_row, ofp_conformation);
					if(z!= new_clique.vertices.size()-1)
						fprintf(ofp_conformation, "\n");
				}
				fprintf(ofp_conformation, "*\n");
			}
		fclose(ofp_conformation);
	}
	vertices_score.clear();
	for (int i = 0; i < bool_PATH_matrix.size(); i++)
		bool_PATH_matrix[i].clear();
	bool_PATH_matrix.clear();
	for (int i = 0; i < neighbors_list.size(); i++)
		neighbors_list[i].clear();
	neighbors_list.clear();
	indexes.clear();
	return(my_clique);
}


