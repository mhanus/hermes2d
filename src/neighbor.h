#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "common.h"
#include "mesh.h"
#include "quad.h"
#include "solution.h"


class HERMES2D_API Neighbor
{
public:

	// Common constructor. The least what user has to provide is active element and mesh.
	// If he wants also function values, he must provide a solution. The space is for improvement of the algorithm
	// for choosing the order on the edge. Both, solution and space, has to be defined over the given mesh.

	Neighbor(Element* e, Mesh* mesh, Solution* sln = NULL, Space* space = NULL);

	~Neighbor();

	//set active edge and compute all needed informations from neighbors, edge is local number of the edge
	void set_active_edge(int edge);

	// number of neighbor elements on edge
	int number_of_neighbs();

	// return array of transformations of neighbor or central el.
	int* get_transformations(int part_edge);

	// return function values of neighbor at integration points
	scalar* get_fn_values_neighbor(int part_edge);

	// return function values of active at integration points
	scalar* get_fn_values_central(int part_edge);

	// return number of integration points;
	int get_n_integ_points(int part_edge);

	// return pointer to the vector of neighbors id.
	std::vector<int>* get_neighbors();

private:
	const static int max_n_trans = 20;    //number of allowed transformations, see "push_transform" in transform.h

	int n_neighbors; // number of neighbors
	Quad2D* quad;
	Mesh* mesh;
	Element* central_el; // active element
	Element* neighb_el;  // actual neighbor we are working with
	Solution* sol;
	Space* space;
	int transformations[max_n_trans][max_n_trans];	// table of transformations for all neighbors
	int n_trans[max_n_trans];  // number of transformations for every neighbor;
	int active_edge;			     // edge where we are searching for neighbors
	int neighbor_edge;		   	// edge of the neighbor respective to active_edge
	scalar* fn_values[max_n_trans]; //function values for active element
	int np[max_n_trans];						// number of integration points for every neighbor
	scalar* fn_values_neighbor[max_n_trans]; //function values for active element
	int central_order; //order of active element
	int neighbor_order; //order of active element
	// vector containing id's of all neighbors
	std::vector<int> neighbors_id;

	// way up for finding neighbor element, from smaller to larger
	void finding_act_elem( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices);

	// way down for finding neighbor elements, from larger to smaller
	void finding_act_elem( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert);

	// setting the sequence of function values of neighbor in same direction as on central element
	void set_correct_direction(int parent1, int parent2, int part_of_edge);

	// set order on the edge. Depends on if the space is given.
	int get_max_order();

	int get_edge_order(Element* e, int edge);

	int get_edge_order_internal(Node* en);

	// cleaning before usage of given edge
	void clean_all();

	enum Trans_flag{
		H2D_NO_TRANSF = 0,  // don't use any transformation, the edge has on both sides active element
		H2D_WAY_DOWN = 1, 	// transformation of active element, against the edge neighbor has some sons
		H2D_WAY_UP = 2			// transformation of neighbor element, central element is son
	};

	// fill function values of central and neighbor element
	void set_fn_values(Trans_flag flag);

	int solution_flag;  // if 1 then function values are computed.

};


#endif /* NEIGHBOR_H_ */
