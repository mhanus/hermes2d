#include "neighbor.h"

NeighborSearch::NeighborSearch(Element* e, Mesh* given_mesh, MeshFunction* given_sln, Space* given_space)
{
	central_el = e;
	mesh = given_mesh;
	sol = given_sln;
	space = given_space;
	if(sol != NULL){
		quad = sol->get_quad_2d();
		sol->set_active_element(central_el);
		central_order = sol->get_fn_order();
	}
	else{
		quad = NULL;
		central_order = -1;
	}
	for(int i = 0;  i < max_n_trans; i++)
	{
		fn_values[i] = NULL;
		fn_values_neighbor[i] = NULL;
	}

	max_of_orders = -1;

	n_neighbors = 0;
	neighbors_id.reserve(20 * e->nvert);
	neighbors.reserve(20);
	values_central.reserve(20);
	values_neighbor.reserve(20);

	neighbor_edges.reserve(20);
};


NeighborSearch::~NeighborSearch()
{
	for(int i = 0; i < max_n_trans; i++){
		if(fn_values[i] != NULL)
			delete[] fn_values[i];
		if(fn_values_neighbor[i] != NULL)
			delete[] fn_values_neighbor[i];
	}

	values_central.clear();
	values_neighbor.clear();
	neighbor_edges.clear();
	orders.clear();
	neighbors.clear();
};


void NeighborSearch::set_active_edge(int edge)
{
	// Erase all data from previous edge or element.
	clean_all();

	active_edge = edge;

	debug_log("central element: %d", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		// Finding active element.
		neighb_el = central_el->get_neighbor(active_edge);

		//First case : The neighboring element is of the same size as the central one.
		if (neighb_el != NULL)
		{
			debug_log("active neighbor el: %d", neighb_el->id);
			// Finding the correct edge.
			for (int j = 0; j < neighb_el->nvert; j++)
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
					neighbor_edge = j;
					break;
				}

			// Get orientation of neighbor edge.
			int p1 = central_el->vn[active_edge]->id;
			int p2 = central_el->vn[(active_edge + 1) % central_el->nvert]->id;

			NeighborEdgeInfo local_edge_info;
			local_edge_info.local_num_of_edge = neighbor_edge;
			local_edge_info.orientation = direction_neighbor_edge(p1, p2, 0);

			// Add the local_edge_info into the vector.
			neighbor_edges.push_back(local_edge_info);

			// There is only one neighbor in this case.
			n_neighbors = 1;

			// Add neighbor id to neighbors_id.
			neighbors_id.push_back(neighb_el->id);
			neighbors.push_back(neighb_el);
			
			// No need for transformation, since the neighboring element is of the same size.
			way_flag = H2D_NO_TRANSF;
			
			if(sol != NULL)
				compute_fn_values();
		} 
		else
		{
			Node* vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,	central_el->en[active_edge]->p2);
			int orig_vertex_id[2];
			orig_vertex_id[0] = central_el->vn[active_edge]->id;
			orig_vertex_id[1]	= central_el->vn[(active_edge + 1) % central_el->nvert]->id;
			if (vertex == NULL)
			{
				// way up
				Element* parent = central_el->parent;

				// array containing vertices we went through
				Node** road_vertices = new Node*[max_n_trans]; 
				// number of used vertices
				int n_road_vertices = 0; 

				for (int j = 0; j < max_n_trans; j++)
					road_vertices[j] = NULL;

				finding_act_elem_up(parent, orig_vertex_id, road_vertices, n_road_vertices);

				delete[] road_vertices;

				way_flag = H2D_WAY_UP;
				
				if(sol != NULL)
						compute_fn_values();
			} 
			else
			{
				// way down
				int road[max_n_trans]; // array for temporal transformation
				int n_road = 0; // number of used transformations

				finding_act_elem_down( vertex, orig_vertex_id, road, n_road);

				way_flag = H2D_WAY_DOWN;
				
				if(sol != NULL)
						compute_fn_values();
				
				debug_log("number of neighbors on the way down: %d ", n_neighbors);
			}
		}
	}
	else
		error("The given edge isn't inner");
};



//way up
/*! \brief Function for finding "bigger" neighbor.
 *
 * We use recurrence in this way.
 * If the neighbor is "bigger" then this means central element is descendant of some inactive elements. We go threw this parents and
 * stop when against an edge, which has same local number as the original edge, we have active element.
 * Important is that all sons have same orientation as parent, so local number of the edge is same.
 */
void NeighborSearch::finding_act_elem_up( Element* elem, int* orig_vertex_id, Node** road_vertices, int n_road_vertices)
{
	Node* edge = NULL;
	Node* vertex = NULL;
	int p1, p2; // id of parents of the edge

	// Order parents in direction of parent element (needed for transformation of solution).
	p1 = elem->vn[active_edge]->id;
	p2 = elem->vn[(active_edge + 1) % elem->nvert]->id;

	int id_of_par_orient_1 = p1;
	int id_of_par_orient_2 = p2;

	// Find if between parents p1 and p2 is used edge (is used by neighbor element).
	edge = mesh->peek_edge_node(p1, p2);

	// When we are on parent, we take middle vertex on the edge and add it to road_vertices. This is for consequent transformation of solution
	// on neighbor element.
	vertex = mesh->peek_vertex_node(p1, p2);
	if(vertex != NULL)
	{
		if (n_road_vertices == 0)
			road_vertices[n_road_vertices++] = vertex;
		else
			if (n_road_vertices == max_n_trans - 1)
				error("n_road_vertices == max_n_trans - 1 in NeighborSearch::finding_act_elem_up");
			else
				if(road_vertices[n_road_vertices - 1]->id != vertex->id)
					road_vertices[n_road_vertices++] = vertex;
	}
	// If we have not found the parent neighboring an active element.
	if ((edge == NULL) || (central_el->en[active_edge]->id == edge->id))
		finding_act_elem_up(elem->parent, orig_vertex_id, road_vertices, n_road_vertices);
	else
		for (int i = 0; i < 2; i++)
		{
			// this condition test if on one of sides is some element and if the element is active, because it may happen that
			// something is found even thought it's not an active element.
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1))
			{
				// Getting to the correct edge.
				debug_log("way up neighbor: %d", edge->elem[i]->id);
				neighb_el = edge->elem[i];

				// Finding the correct edge.
				neighbor_edge = -1;
				for(int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
					{
						neighbor_edge = j;
						break;
					}
				if(neighbor_edge == -1) error("edge wasn't found");
			
				Node* n = NULL;
				// Filling the array containing number of transformation needed for this neighbor.
				// There is always only one neighbor (since the central element is the smaller of the two neighboring ones).
				assert(n_neighbors == 0);
				n_trans[n_neighbors] = n_road_vertices;
	
				//debugging
				for(int k = 0 ; k < n_road_vertices; k++)
					debug_log("vertices on the way: %d", road_vertices[k]->id);
				debug_log("\n");
				
				// Go between the inactive parent neighboring an active element and the central 
				// element and save correct transformation along the way down.
				for(int j = n_road_vertices - 1; j > 0; j-- )
				{
						n = mesh->peek_vertex_node(road_vertices[j]->id, p1);
						if(n == NULL)
						{
							transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
							p1 = road_vertices[j]->id;
						}
						else
						{
								if(n->id == road_vertices[j-1]->id)
								{
									transformations[n_neighbors][n_road_vertices - j - 1] = (neighbor_edge + 1) % neighb_el->nvert;
									p2 = road_vertices[j]->id;
								}
								else
								{
									transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
									p1 = road_vertices[j]->id;
								}
						}
				}

				// Final transformation on central element.
				int test = 0;
				if (orig_vertex_id[0] == road_vertices[0]->id)
					test = 1;

				if(test == 1){
					transformations[n_neighbors][n_road_vertices - 1] = neighbor_edge;
				}
				else{
					transformations[n_neighbors][n_road_vertices - 1] = (neighbor_edge + 1) % neighb_el->nvert;
				}

				NeighborEdgeInfo local_edge_info;
				local_edge_info.local_num_of_edge = neighbor_edge;
				local_edge_info.orientation = direction_neighbor_edge(id_of_par_orient_1, id_of_par_orient_2, 0);

				// Add the local_edge_info into the vector.
				neighbor_edges.push_back(local_edge_info);

				// There is only one neighbor.
				n_neighbors = 1;

				// Add neighbors to their respective vectors.
				neighbors_id.push_back(neighb_el->id);
				neighbors.push_back(neighb_el);
			}
		}
};

/*! \brief On active edge we have more neighbors. Gives us information from all neighbors.
 *
 *	Again we use recurrence in this way. In every step we take middle vertex of the edge (starting with active edge). This vertex split the edge
 *	on two parts. On every part (an edge) we test if the new edge is active. If not, the middle vertex is found and the method is called
 *	again with this new vertex on this part.
 */
//way down
void NeighborSearch::finding_act_elem_down( Node* vertex, int* par_vertex_id, int* road, int n_road)
{
	int son = vertex->id;
	int parents[2];
	parents[0] = par_vertex_id[0];
	parents[1] = par_vertex_id[1];

	for (int i = 0; i < 2; i++)
	{
		road[n_road] = (active_edge + i) % central_el->nvert;

		Node* edge = mesh->peek_edge_node(son, parents[i]);

		// Test if edge is inactive. Means there is no active element on each side.
		if (edge == NULL)
		{
			Node * n = mesh->peek_vertex_node(son, parents[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else
			{
				if(i == 0) 
					par_vertex_id[1] = son;
				else 
					par_vertex_id[0] = son;
				
				finding_act_elem_down( n, par_vertex_id, road, ++n_road);
			}
		} 
		else
		{
			// Test if on one of sides is active element.
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1)
					{
							neighb_el = mesh->get_element(edge->elem[j]->id);

						  debug_log("way down neighbor: %d", edge->elem[j]->id);

							// Getting to correct edge.
							neighbor_edge = -1;
							for(int k = 0; k < neighb_el->nvert; k++)
								if(neighb_el->en[k] == edge)
								{
									neighbor_edge = k;
									break;
								}
							if(neighbor_edge == -1) error("edge wasn't found");

							// Filling transformation.
							for(int k = 0; k <= n_road; k++) transformations[n_neighbors][k] = road[k];

							// + 1 is because how to n_road is computed it's one less then number of transformations.
							n_trans[n_neighbors] = n_road + 1;

							NeighborEdgeInfo local_edge_info;
							local_edge_info.local_num_of_edge = neighbor_edge;
							local_edge_info.orientation = direction_neighbor_edge(parents[0], parents[1], i);

							// Add the local_edge_info into the vector.
							neighbor_edges.push_back(local_edge_info);

							// Raise number of neighbors.
							n_neighbors = n_neighbors++;

							// Add neighbors to vectors.
							neighbors_id.push_back(neighb_el->id);
							neighbors.push_back(neighb_el);
					}
			}
		}
	}
};

/*! \brief Fill function values / derivatives for the central element and one of its neighbors.
 *
 */

void NeighborSearch::set_fn_values(Trans_flag flag, int neigh)
{

	int number_integ_points = 0;

	switch(flag)
	{
		case H2D_NO_TRANSF:
			{
				sol->set_active_element(neighb_el);

				int max_order;
				if (max_of_orders == -1)
				{
					neighbor_order = sol->get_fn_order();
					max_order = get_max_order();
				}
				else
					max_order = max_of_orders;

				orders.push_back(max_order);

				int eo = quad->get_edge_points(neighbor_edge, max_order);
				number_integ_points = quad->get_num_points(eo);

				scalar* local_fn_values_n = new scalar[number_integ_points];
				scalar* local_fn_values_c = new scalar[number_integ_points];

				// Fill function values of neighbor.
				sol->set_quad_order(eo);

				RefMap* rm = sol->get_refmap();
				Func<scalar>* func;
				func = init_fn(sol, rm, eo);
				values_neighbor.push_back(func);

				for(int i = 0; i < number_integ_points ; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
				fn_values_neighbor[neigh] = local_fn_values_n;


				// Fill the central.

				sol->set_active_element(central_el);
				eo = quad->get_edge_points(active_edge, max_order);
				sol->set_quad_order(eo);
				for(int i = 0; i < number_integ_points ; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
				fn_values[neigh] = local_fn_values_c;

				Func<scalar>* func1;
				func1 = init_fn(sol, rm, eo);
				values_central.push_back(func1);


				break;
			}
		case H2D_WAY_DOWN:
			{
			sol->set_active_element(neighb_el);

			int max_order;

			if (max_of_orders == -1)
			{
				neighbor_order = sol->get_fn_order();
				max_order = get_max_order();
			}
			else
				max_order = max_of_orders;

			orders.push_back(max_order);

			int eo = quad->get_edge_points(neighbor_edge, max_order);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// Fill function values of neighbor element.

			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
			fn_values_neighbor[neigh] = local_fn_values_n;

			RefMap* rm = sol->get_refmap();
			Func<scalar>* func;
			func = init_fn(sol, rm, eo);
			values_neighbor.push_back(func);

			// Fill the central element.

			sol->set_active_element(central_el);

			// Transform central element on appropriate part.
			for(int i = 0; i < n_trans[neigh]; i++){
				sol->push_transform(transformations[neigh][i]);
			}
			eo = quad->get_edge_points(active_edge, max_order);
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
			fn_values[neigh] = local_fn_values_c;

			Func<scalar>* func1;
			func1 = init_fn(sol, rm, eo);
			values_central.push_back(func1);

			break;
			}
		case H2D_WAY_UP:
			{
			sol->set_active_element(neighb_el);

			// Transform neighbor element on appropriate part.
			for(int i = 0; i < n_trans[neigh]; i++){
				sol->push_transform(transformations[neigh][i]);
			}

			int max_order;

			if (max_of_orders == -1)
			{
				neighbor_order = sol->get_fn_order();
				max_order = get_max_order();
			}
			else
				max_order = max_of_orders;

			orders.push_back(max_order);

			int eo = quad->get_edge_points(neighbor_edge, max_order);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// Fill function values of neighbor element.

			sol->set_quad_order(eo);

			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
			fn_values_neighbor[neigh] = local_fn_values_n;

			RefMap* rm = sol->get_refmap();
			Func<scalar>* func;
			func = init_fn(sol, rm, eo);
			values_neighbor.push_back(func);

			// Fill the central element.

			sol->set_active_element(central_el);
			eo = quad->get_edge_points(active_edge, max_order);
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
			fn_values[neigh] = local_fn_values_c;

			Func<scalar>* func1;
			func1 = init_fn(sol, rm, eo);
			values_central.push_back(func1);

			break;
			}
		default:
			error("wasn't find scheme for getting correctly transformed function values");
	}


	// Test if number of function values was assigned.
	if(number_integ_points == 0)
		error("number of integration points is 0");

	np[neigh] = number_integ_points;

	// Reset transformations.
	sol->reset_transform();

};

// Here part_of_edge can be only 0 or 1. Different than 0 is meaningful only for way down.
// Fills the info about orientation into edge_info.

int NeighborSearch::direction_neighbor_edge(int parent1, int parent2, int part_of_edge)
{
	if (part_of_edge == 0)
	{
		if (neighb_el->vn[neighbor_edge]->id != parent1)
			return 1; // means the orientation is conversely
	}
	else
		if (neighb_el->vn[neighbor_edge]->id == parent2)
			return 1; // means the orientation is conversely
	return 0;
}

void NeighborSearch::reverse_vector(scalar* vector, int n)
{
	scalar local_value;

	for(int i = 0; i < n / 2; i++)
	{
		local_value = vector[i];
		vector[i] = vector[n - i - 1];
		vector[n - i - 1] = local_value;
	}
}


void NeighborSearch::set_correct_direction(int index)
{
		scalar local_value = 0;
		Func<scalar>* local_fun = values_neighbor[index];

		reverse_vector(fn_values_neighbor[index], np[index]);

		// change of the orientation for all values at local_fun.
			if(local_fun->nc == 1){
				reverse_vector(local_fun->val, np[index]);
				reverse_vector(local_fun->dx, np[index]);
				reverse_vector(local_fun->dy, np[index]);
			}
			else if(local_fun->nc == 2){
				reverse_vector(local_fun->val0, np[index]);
				reverse_vector(local_fun->val1, np[index]);
				reverse_vector(local_fun->curl, np[index]);
			}
};

int NeighborSearch::get_max_order()
{
	if(space != NULL){
		central_order = get_edge_order(central_el, active_edge);
		neighbor_order = get_edge_order(neighb_el, neighbor_edge);
	}
		return std::max(central_order, neighbor_order);
};



void NeighborSearch::clean_all()
{
	active_edge = -1;

	n_neighbors = 0;
	neighb_el = NULL;
	neighbor_edge = -1;
	neighbor_order = -1;
	way_flag = -1;

	for(int i = 0; i < max_n_trans; i++)
	{
		n_trans[i] = 0;
		np[i] = 0;
		if(fn_values[i] != NULL)
			delete[] fn_values[i];
		if(fn_values_neighbor[i] != NULL)
			delete[] fn_values_neighbor[i];

		fn_values[i] = NULL;
		fn_values_neighbor[i]	= NULL;

		for(int j = 0; j < max_n_trans; j++)
			transformations[i][j] = -1;
	}
	values_central.clear();
	values_neighbor.clear();
	neighbor_edges.clear();
	orders.clear();
	neighbors.clear();
};



int NeighborSearch::get_number_of_neighbs()
{
	if(n_neighbors == 0) error("called before setting common edge");
	else return n_neighbors;
};

int* NeighborSearch::get_transformations(int part_edge)
{
	return transformations[part_edge];
};

scalar* NeighborSearch::get_fn_values_central(int part_edge)
{
	return fn_values[part_edge];
};
scalar* NeighborSearch::get_fn_values_neighbor(int part_edge)
{
	return fn_values_neighbor[part_edge];
};

Func<scalar>* NeighborSearch::get_values_central(int part_edge)
{
	return values_central[part_edge];
};

Func<scalar>* NeighborSearch::get_values_neighbor(int part_edge)
{
	return values_neighbor[part_edge];
};

int NeighborSearch::get_n_integ_points(int part_edge)
{
	return np[part_edge];
};


std::vector<int>* NeighborSearch::get_neighbors()
{
	return &neighbors_id;
};

std::vector<int>* NeighborSearch::get_orders()
{
	return &orders;
};


int NeighborSearch::get_number_neighb_edge(int part_edge)
{
	if(part_edge >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[part_edge].local_num_of_edge;
};

int NeighborSearch::get_orientation_neighb_edge(int part_edge)
{
	if(part_edge >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[part_edge].orientation;
};


// methods for getting order on the edge from space. Originally taken from class Space.

int NeighborSearch::get_edge_order(Element* e, int edge)
{
  Node* en = e->en[edge];
  if (en->id >= space->nsize || edge >= (int)e->nvert) return 0;

  if (space->ndata[en->id].n == -1)
    return get_edge_order_internal(space->ndata[en->id].base); // constrained node
  else
    return get_edge_order_internal(en);
}


int NeighborSearch::get_edge_order_internal(Node* en)
{
  assert(en->type == H2D_TYPE_EDGE);
  Element** e = en->elem;
  int o1 = 0, o2 = 0;
  assert(e[0] != NULL || e[1] != NULL);

  if (e[0] != NULL)
  {
    if (e[0]->is_triangle() || en == e[0]->en[0] || en == e[0]->en[2])
      o1 = H2D_GET_H_ORDER(space->edata[e[0]->id].order);
    else
      o1 = H2D_GET_V_ORDER(space->edata[e[0]->id].order);
  }

  if (e[1] != NULL)
  {
    if (e[1]->is_triangle() || en == e[1]->en[0] || en == e[1]->en[2])
      o2 = H2D_GET_H_ORDER(space->edata[e[1]->id].order);
    else
      o2 = H2D_GET_V_ORDER(space->edata[e[1]->id].order);
  }

  if (o1 == 0) return o2 == 0 ? 0 : o2;
  if (o2 == 0) return o1;
  return std::max(o1, o2);
}

void NeighborSearch::compute_fn_values()
{
	for(int i = 0; i < n_neighbors; i++)
	{
		neighb_el = neighbors[i];
		neighbor_edge = neighbor_edges[i].local_num_of_edge;
		set_fn_values((Trans_flag)way_flag, i);
		if(neighbor_edges[i].orientation == 1)
			set_correct_direction(i);
	}
}

void NeighborSearch::set_order_of_integration(int order)
{
	if(order < 0)
		error("given order is negative.");
	max_of_orders = order;
}

void NeighborSearch::set_solution(MeshFunction* solution)
{
	sol = solution;
	quad = sol->get_quad_2d();
	sol->set_active_element(central_el);
	central_order = sol->get_fn_order();

	if(active_edge == -1)
		error("the common(active) edge wasn't set.");
	orders.clear();
	values_central.clear();
	values_neighbor.clear();
	compute_fn_values();
}