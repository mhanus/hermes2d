#include "neighbor.h"

void Neighbor::set_active_edge(int edge)
{
	//erase all data from previous edge or element
	clean_all();

	active_edge = edge;

	printf("element: %d \n\n\n", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		neighb_el = central_el->get_neighbor(active_edge);
		if (neighb_el != NULL)
		{
			for (int j = 0; j < neighb_el->nvert; j++)
			{
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
					neighbor_edge = j;
					set_fn_values(H2D_NO_TRANSF);
					//raise the number of neighbors
					n_neighbors = n_neighbors++;
					printf("aktivni elem: %d \n", neighb_el->id);
				}
			}
		} else
		{
			Node* vertex = NULL;
			vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,	central_el->en[active_edge]->p2);
			int orig_vertex_id[2];
			orig_vertex_id[0] = central_el->vn[active_edge]->id;
			orig_vertex_id[1]	= central_el->vn[(active_edge + 1) % central_el->nvert]->id;
			if (vertex == NULL)
			{
				//way up

				Element* parent = NULL;
				parent = central_el->parent;

				Node** road_vertices;
				road_vertices = new Node*[max_n_trans];
				int n_road_vertices = 0; //number of used vertices

				for (int j = 0; j < max_n_trans; j++)
					road_vertices[j] = NULL;

				finding_act_elem(parent, active_edge, orig_vertex_id, road_vertices, n_road_vertices);

				delete[] road_vertices;
			} else
			{
				//way down

				int road[max_n_trans]; //array for temporal transformation
				int n_road = 0; //number of used transformations

				finding_act_elem( vertex, orig_vertex_id, road, n_road,	active_edge, central_el->nvert);
				printf("number of neighbors: %d", n_neighbors);
			}
		}
	}
};



//way up

void Neighbor::finding_act_elem( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices)
{

	Node* edge = NULL;
	Node* vertex = NULL;
	int p1, p2; //id of parents of the edge

	//order parents in direction of parent element
	p1 = elem->vn[edge_num]->id;
	p2 = elem->vn[(edge_num + 1) % elem->nvert]->id;

	edge = mesh->peek_edge_node(p1, p2);
	vertex = mesh->peek_vertex_node(p1, p2);
	road_vertices[n_road_vertices] = vertex;
	n_road_vertices = n_road_vertices++;

	if (edge == NULL)
	{
		finding_act_elem(elem->parent, edge_num, orig_vertex_id, road_vertices, n_road_vertices);
	}
	else
		for (int i = 0; i < 2; i++)
		{
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1)){

				//getting to correct edge
				printf("way up neighb: %d", edge->elem[i]->id);
				neighb_el = edge->elem[i];
				neighbor_edge = -1;
				for(int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
						neighbor_edge = j;
				if(neighbor_edge == -1) error("edge wasn't found");

				Node* n = NULL;

				//set active the neighbour
//				sol->set_active_element(neighb_el);
				// set number of tranformations
				//not sure about it !!!!!!!!!!!!!!!
				n_trans[n_neighbors] = n_road_vertices;

				// go threw between elements and set correct transformation
				for(int j = n_road_vertices; j > 0; j-- ){
					if(road_vertices[j] == NULL){
						if(j > 0)
							continue;
					}
					else{
						n = mesh->peek_vertex_node(road_vertices[j]->id, p1);
						if(n == NULL){
							n = mesh->peek_vertex_node(road_vertices[j]->id, p2);
							transformations[n_neighbors][n_road_vertices - j] = neighbor_edge;
//							sol->push_transform(neighbor_edge);
							p1 = n->id;
						}
						else{
								if(n->id == road_vertices[j-1]->id){
									transformations[n_neighbors][n_road_vertices - j] = (neighbor_edge + 1) % neighb_el->nvert;
//									sol->push_transform((neighbor_edge + 1) % neighb_el->nvert);
									p2 = n->id;
								}
								else{
									n = mesh->peek_vertex_node(road_vertices[j]->id, p2);
									transformations[n_neighbors][n_road_vertices - j] = neighbor_edge;
//									sol->push_transform(neighbor_edge);
									p1 = n->id;
								}
						}
					}
				}

				// final transformation on active element
				int test = 0;
				if (orig_vertex_id[0] == road_vertices[0]->id)
					test = 1;

				if(test == 1){
					transformations[n_neighbors][n_road_vertices] = neighbor_edge;
//					sln->push_transform(neighbor_edge);
				}
				else{
					transformations[n_neighbors][n_road_vertices] = (neighbor_edge + 1) % neighb_el->nvert;
//					sln->push_transform((neighbor_edge + 1) % neighb_el->nvert);
				}

				set_fn_values(H2D_WAY_UP);

				set_correct_direction(p1, p2, i);

				// raise the number of neighbors

				n_neighbors = n_neighbors++;
			}
		}
};


//way down
void Neighbor::finding_act_elem( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert)
{
	int son;
	int parents[2];

	Node* edge = NULL;
	Node* n = NULL;
	neighb_el = NULL;
	son = vertex->id;

	parents[0] = par_vertex_id[0];
	parents[1] = par_vertex_id[1];

	for (int i = 0; i < 2; i++)
	{
		road[n_road] = (use_edge + i) % n_vert;

		edge = mesh->peek_edge_node(son, parents[i]);
		//test if edge is active, means on one of sides can!! be active element
		if (edge == NULL)
		{
			n = mesh->peek_vertex_node(son, parents[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else{
				if(i == 0) par_vertex_id[1] = son;
				else par_vertex_id[0] = son;

				n_road = n_road++;
				finding_act_elem( n, par_vertex_id, road, n_road, use_edge, n_vert);
			}
		} else
			//test if on one of sides is active element
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1){
						//do something
							printf("way down neighb: %d", edge->elem[j]->id);
							neighb_el = mesh->get_element(edge->elem[j]->id);

							neighbor_edge = -1;
							for(int k = 0; k < neighb_el->nvert; k++)
								if(neighb_el->en[k] == edge)
									neighbor_edge = k;
							if(neighbor_edge == -1) error("edge wasn't found");

							//filling transformation
							for(int k = 0; k <= n_road; k++) transformations[n_neighbors][k] = road[k];

							// + 1 is because how to n_road is computed it's one less then number of transformations
							n_trans[n_neighbors] = n_road + 1;

							set_fn_values(H2D_WAY_DOWN);


/*							for(int k = 0; k <= n_road; k++) neighbs[*n_neighbs].transformations[k] = road[k];
							int neighb_order = sln->get_fn_order();

							neighbs[*n_neighbs].max_order = std::max(active_order, neighb_order);
							// in future need to set correct np and fna_values acc. max_order
*/
							// correct direction, means the orientation of direction of function values on neighbor has to be same
							// as on central element

							set_correct_direction(parents[0], parents[1], i);

/*							int test = 0;
							int neighb_first_vertex = neighb->vn[neighb_edge]->id;
							if(i == 0){
								if(neighb_first_vertex != parents[0])
									test = 1; // means the orientation is conversely
							}
							else{
								if(neighb_first_vertex == parents[1])
									test = 1; // means the orientation is conversely
							}
							neighbs[*n_neighbs].n_fn_values = np;
							if(test == 1)
								for(int k = 0; k <= np; k++) neighbs[*n_neighbs].fn_values[k] = fn[np-k];
							else
								for(int k = 0; k <= np; k++) neighbs[*n_neighbs].fn_values[k] = fn[k];
							*n_neighbs = (*n_neighbs)++;
*/

				// raise the number of neighbors

				n_neighbors = n_neighbors++;

					}
			}
	}
};

void Neighbor::set_fn_values(Trans_flag flag){

	int number_integ_points = 0;

	switch(flag){
		case H2D_NO_TRANSF:
			break;

		case H2D_WAY_DOWN:
			{
			sol->set_active_element(neighb_el);
			neighbor_order = sol->get_fn_order();
			int max_order = std::max(active_order, neighbor_order);

			// now it takes the "original" max order
			int eo = quad->get_edge_points(neighbor_edge);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// fill function values of neighbor
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
			fn_values_neighbor[n_neighbors] = local_fn_values_n;

			sol->set_active_element(central_el);

			//transform central on appropriate part
			for(int i = 0; i < n_trans[n_neighbors]; i++){
				sol->push_transform(transformations[n_neighbors][i]);
			}
			//fill the central
			eo = quad->get_edge_points(active_edge);
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
			fn_values[n_neighbors] = local_fn_values_c;

			break;
			}
		case H2D_WAY_UP:
			{
			sol->set_active_element(neighb_el);

			//transform neighbor on appropriate part
			for(int i = 0; i < n_trans[n_neighbors]; i++){
				sol->push_transform(transformations[n_neighbors][i]);
			}

			neighbor_order = sol->get_fn_order();
			int max_order = std::max(active_order, neighbor_order);

			// now it takes the "original" max order
			int eo = quad->get_edge_points(neighbor_edge);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// fill function values of neighbor
			sol->set_quad_order(eo);

			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
			fn_values_neighbor[n_neighbors] = local_fn_values_n;

			//fill the central
			sol->set_active_element(central_el);
			eo = quad->get_edge_points(active_edge);
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
			fn_values[n_neighbors] = local_fn_values_c;

			break;
			}
		default:
			error("wasn't find scheme for getting correctly transformed function values");
	}


	//number of function values
	if(number_integ_points == 0)
		error("number of integration points is 0");

	np[n_neighbors] = number_integ_points;

	//reset transformations
	sol->reset_transform();

};

// correct direction, means the orientation of direction of function values on neighbor has to be same
// as on central element

void Neighbor::set_correct_direction(int parent1, int parent2, int part_of_edge)
{

	int test = 0;
	int neighb_first_vertex = neighb_el->vn[neighbor_edge]->id;
	if (part_of_edge == 0)
	{
		if (neighb_first_vertex != parent1)
			test = 1; // means the orientation is conversely
	} else
	{
		if (neighb_first_vertex == parent2)
			test = 1; // means the orientation is conversely
	}

	if (test == 1)
	{
		scalar local_fn_value = 0;

//		for(int i = 0; i < np[n_neighbors]; i++) local_fn_values[i] = fn_values_neighbor[n_neighbors][i];
		for (int i = 0; i < (np[n_neighbors] - 1) / 2; i++){
			local_fn_value = fn_values_neighbor[n_neighbors][i];
			fn_values_neighbor[n_neighbors][i] = fn_values_neighbor[n_neighbors][np[n_neighbors] - i];
			fn_values_neighbor[n_neighbors][np[n_neighbors] - i] = local_fn_value;
		}
	}
};

void Neighbor::clean_all()
{
	active_edge = -1;

	n_neighbors = 0;
	neighb_el = NULL;
	neighbor_edge =-1;
	neighbor_order = -1;

	for(int i = 0; i < max_n_trans; i++)
	{
		n_trans[i] = 0;
		np[i] = 0;

		if(fn_values[i] != NULL)
			delete[] fn_values[i];
		if(fn_values_neighbor[i] != NULL)
			delete[] fn_values_neighbor[i];

		for(int j = 0; j < max_n_trans; j++)
			transformations[i][j] = -1;
	}
};
