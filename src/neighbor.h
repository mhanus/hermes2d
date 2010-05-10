#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "common.h"
#include "mesh.h"
#include "quad.h"
#include "solution.h"
#include "forms.h"


class HERMES2D_API NeighborSearch
{
public:

	// Common constructor. The least what user has to provide is active element and mesh.
	// If he wants also function values, he must provide a solution. The space is for improvement of the algorithm
	// for choosing the order on the edge. Both, solution and space, has to be defined over the given mesh.

	NeighborSearch(Element* e, Mesh* mesh, MeshFunction* sln = NULL, Space* space = NULL);

	~NeighborSearch();

	//set active edge and compute all needed informations from neighbors, edge is local number of the edge
	void set_active_edge(int edge);

	// number of neighbor elements on edge
	int get_number_of_neighbs();

	// return array of transformations of neighbor or central el.
	int* get_transformations(int part_edge);

	// return function values of neighbor at integration points
	scalar* get_fn_values_neighbor(int part_edge);

	// return function values of active at integration points
	scalar* get_fn_values_central(int part_edge);

	// return pointer to function which contains all information from central
	Func<scalar>* get_values_central(int part_edge);

	// return pointer to function which contains all information from neighbor
	Func<scalar>* get_values_neighbor(int part_edge);


	// return number of integration points;
	int get_n_integ_points(int part_edge);

	// return pointer to the vector of neighbors id.
	std::vector<int>* get_neighbors();

	// return local number of neighbor edge
	int get_number_neighb_edge(int part_edge);

	// return orientation of neighbor edge
	int get_orientation_neighb_edge(int part_edge);

	// each member of the vector contains maximum of orders of central and neighbor element.
	std::vector<int>* get_orders();

	// set the order of the integration 
	void set_order_of_integration(int order);

	// set meshfunction (solution), this method is used in case you want to get values of different solutions over the same edge.
	// Don't have to make another instance of the class.
	void set_solution(MeshFunction* solution);


private:
	const static int max_n_trans = 20;    //number of allowed transformations, see "push_transform" in transform.h

	int n_neighbors; // number of neighbors
	Quad2D* quad;
	Mesh* mesh;
	Element* central_el; // active element
	Element* neighb_el;  // actual neighbor we are working with
	MeshFunction* sol;
	Space* space;
	int transformations[max_n_trans][max_n_trans];	// table of transformations for all neighbors
	int n_trans[max_n_trans];  // number of transformations for every neighbor;
	int active_edge;			     // edge where we are searching for neighbors
	int neighbor_edge;		   	// edge of the neighbor respective to active_edge
	scalar* fn_values[max_n_trans]; //function values for active element
	int np[max_n_trans];						// number of integration points for every neighbor
	scalar* fn_values_neighbor[max_n_trans]; // function values for active element
	int central_order;  // order of active element
	int neighbor_order; // order of active element
	int max_of_orders;  // initial set equal to -1, else set by method set_order_of_integration().
	int way_flag; // this flag holds which way was used on the active edge.


	// vector containing id's of all neighbors
	std::vector<int> neighbors_id;
	std::vector<Element*> neighbors;
	std::vector<int> orders;

	// vectors of all values (function, derivatives, etc.) for central and neighbors
	std::vector<Func<scalar>*> values_central;
	std::vector<Func<scalar>*> values_neighbor;

	// way up for finding neighbor element, from smaller to larger
	void finding_act_elem_up( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices);

	// way down for finding neighbor elements, from larger to smaller
	void finding_act_elem_down( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert);

	// setting the sequence of function values of neighbor in same direction as on central element.
	void set_correct_direction(int index);

	// set order on the edge. Depends on if the space is given.
	int get_max_order();

	int get_edge_order(Element* e, int edge);

	int get_edge_order_internal(Node* en);

	// Just reverse values in vector
	void reverse_vector(scalar* vector, int n);

	// Structure containing all needed info about neighbor's edge
	struct NeighborEdgeInfo{
		NeighborEdgeInfo(){
			local_num_of_edge = -1;
			orientation = -1;
		}
		int local_num_of_edge; // number of the edge at the neighbor element
		int orientation; // if the orientation is same as on active edge is 0 otherwise 1
	};

	// find the orientation of the neighbor edge in relation with central edge.
	void direction_neighbor_edge(int parent1, int parent2, int part_of_edge, NeighborEdgeInfo* edge_info);

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


	std::vector<NeighborEdgeInfo> neighbor_edges; // vector containing all neighbor edges related to active edge

	void compute_fn_values();

};


#endif /* NEIGHBOR_H_ */
