#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "mesh.h"
#include "quad.h"
#include "solution.h"


class PUBLIC_API Neighbor
{
public:
	Neighbor(Element* e, Solution* sln){
		central_el = e;
		sol = sln;
		quad = sln->get_quad_2d();
		mesh = sln->get_mesh();
		sol->set_active_element(central_el);
		active_order = sol->get_fn_order();

		for(int i = 0; i < max_n_trans; i++){
			fn_values[i] = NULL;
			fn_values_neighbor[i] = NULL;			
			}
	}

	//set active edge and compute all needed informations from neighbors, edge is local number of the edge
	void set_active_edge(int edge);
	// number of neighbor elements on edge
	int number_of_neighbs();
	// return pseudo-order for use in getting integration points of order
//	int get_edge_points(int i);
	// return array of transformations of neighbor or central el.
	int* get_transformations(int i);
	// return function values of neighbor at integration points
	double* get_fn_values_neighbor(int i);
	// return function values of active at integration points
	double* get_fn_values_active(int i);

private:
	const static int max_n_trans = 20;    //number of allowed transformations, see push_transform() in transform.h

	int n_neighbors; // number of neighbors
	Quad2D* quad;
	Mesh* mesh;
	Element* central_el; // active element
	Element* neighb_el;  // actual neighbor we are working with
	Solution* sol;
	int transformations[max_n_trans][max_n_trans];	// table of transformations for all neighbors
	int n_trans[max_n_trans];  // number of transformations for every neighbor;
	int active_edge;			     // edge where we are searching for neighbors
	int neighbor_edge;		   	// edge of the neighbor respective to active_edge
	scalar* fn_values[max_n_trans]; //function values for active element
	int np[max_n_trans];						// number of integration points for every neighbor
	scalar* fn_values_neighbor[max_n_trans]; //function values for active element
	int active_order; //order of active element
	int neighbor_order; //order of active element

	// way up for finding neighbor element, from smaller to larger
	void finding_act_elem( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices);

	// way down for finding neighbor elements, from larger to smaller
	void finding_act_elem( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert);

	// setting the sequence of function values of neighbor in same direction as on central element
	void set_correct_direction(int parent1, int parent2, int part_of_edge);

	// cleaning before usage of given edge
	void clean_all();

	enum Trans_flag{
		H2D_NO_TRANSF = 0,  // don't use any transformation, the edge has on both sides active element
		H2D_WAY_DOWN = 1, 	// transformation of active element, against the edge neighbor has some sons
		H2D_WAY_UP = 2			// transformation of neighbor element, central element is son
	};

	// fill function values of central and neighbor element
	void set_fn_values(Trans_flag flag);


};




#endif /* NEIGHBOR_H_ */
