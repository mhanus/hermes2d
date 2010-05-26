#include "H1PeriodicSpace.h"
/*
void H1PeriodicSpace::post_assign()
{
	H1Space::post_assign();

	/// reassign dofs belonging to image edges to those of source edges	

	mnmap src_nodes_map, img_nodes_map; // mapping of periodic nodes to their common markers
  map<int> markers;										// map of unique periodic node markers
  
	Node* bndnd;
	for_all_bnd_nodes(bndnd, mesh) 
	{
		int mrk = -1;
		if (bndnd->marker >= SRC_MARKERS_START && bndnd->marker < IMG_MARKERS_START) {
			printf("SRC %d\n", ndata[bndnd->id].dof);
			mrk = bndnd->marker;
			src_nodes_map.insert(pair<int,Node*>(mrk,bndnd));
		}
		else if (bndnd->marker >= IMG_MARKERS_START) {
			printf("IMG %d\n", ndata[bndnd->id].dof);
			mrk = bndnd->marker-MARKERS_DIFF;
			img_nodes_map.insert(pair<int,Node*>(mrk,bndnd));
		}
		
		if (mrk > 0)
			markers.insert(mrk);
	}
	
	printf("\n%d, %d, %d\n\n", markers.size(), src_nodes_map.size(), img_nodes_map.size());
	
	for (map<int>::iterator markers_it = markers.begin(); markers_it != markers.end(); ++markers_it) 
	{
		// for each periodic node marker, pick the corresponding source and image nodes
		pair<mnmap::iterator, mnmap::iterator> src_nodes = src_nodes_map.equal_range(*markers_it);
		pair<mnmap::iterator, mnmap::iterator> img_nodes = img_nodes_map.equal_range(*markers_it);
				
		for (	mnmap::iterator src_node_it = src_nodes.first, img_node_it = img_nodes.first;
				 	src_node_it != src_nodes.second;
				 	++src_node_it, ++img_node_it	) 
		{
			if (img_node_it == img_nodes.second) 
				error("Invalid");
				
			Node* src_node = src_node_it->second;
			Node* img_node = img_node_it->second;
			printf("%d -> %d\n", src_node->id, img_node->id);
			
			ndata[img_node->id].dof = ndata[src_node->id].dof;
		
			// also reassign dofs belonging to edge bounding vertices
			img_node = mesh->peek_vertex_node(img_node->p1, img_node->p2);
			src_node = mesh->peek_vertex_node(src_node->p1, src_node->p2);
			assert(	(img_node == NULL && src_node == NULL) || 
							(img_node != NULL && src_node != NULL && img_node->bnd && src_node->bnd)	);
		
			if (src_node != NULL) {
				printf("%d -> %d\n", src_node->id, img_node->id);
				ndata[img_node->id].dof = ndata[src_node->id].dof;
			}
		}
	}
}
*/

void H1PeriodicSpace::post_assign()
{
	H1Space::post_assign();

	/// reassign dofs belonging to image edges to those of source edges	

	mnmap src_nodes_map, img_nodes_map; // mapping of periodic nodes to their common markers
  set<int> markers;										// set of unique periodic node markers
  
  
  Element* e;
  for_all_active_elements(e, mesh)
  {
  	int order = get_element_order(e->id);
    if (order > 0)
    {
      for (unsigned int i = 0; i < e->nvert; i++)
      {
      	int mrk = -1;
      	
      	Node* vn1 = e->vn[i];
      	Node* vn2 = e->vn[e->next_vert(i)];
        Node* en = mesh->peek_edge_node(vn1->id, vn2->id);
      	
      	if (is_nat_bnd(en))
        {
 	        if (en->marker >= SRC_MARKERS_START && en->marker < IMG_MARKERS_START) {					
						mrk = en->marker;
						src_nodes_map.insert(pair<int,Node*>(mrk,en));
						
						if (is_nat_bnd(vn1))
							src_nodes_map.insert(pair<int,Node*>(mrk,vn1));
						if (is_nat_bnd(vn2))
							src_nodes_map.insert(pair<int,Node*>(mrk,vn2));
					}
					else if (en->marker >= IMG_MARKERS_START) {													
						mrk = en->marker - MARKERS_DIFF;
						img_nodes_map.insert(pair<int,Node*>(mrk,en));
						
						if (is_nat_bnd(vn1))
							img_nodes_map.insert(pair<int,Node*>(mrk,vn1));
						if (is_nat_bnd(vn2))
							img_nodes_map.insert(pair<int,Node*>(mrk,vn2));
					}
        }    		
        
		    if (mrk > 0)
					markers.insert(mrk);
		  }
    }
  }
	
	printf("\n%d, %d, %d\n\n", markers.size(), src_nodes_map.size(), img_nodes_map.size());
	
	for (set<int>::iterator markers_it = markers.begin(); markers_it != markers.end(); ++markers_it) 
	{
		// for each periodic node marker, pick the corresponding source and image nodes
		pair<mnmap::iterator, mnmap::iterator> src_nodes = src_nodes_map.equal_range(*markers_it);
		pair<mnmap::iterator, mnmap::iterator> img_nodes = img_nodes_map.equal_range(*markers_it);
				
		for (	mnmap::iterator src_node_it = src_nodes.first, img_node_it = img_nodes.first;
				 	src_node_it != src_nodes.second;
				 	++src_node_it, ++img_node_it	) 
		{
			if (img_node_it == img_nodes.second) 
				error("Invalid");
				
			Node* src_node = src_node_it->second;
			Node* img_node = img_node_it->second;
			printf("%d <-> %d\n", ndata[img_node->id].dof, ndata[src_node->id].dof);
			
			ndata[img_node->id].dof = ndata[src_node->id].dof;
		}
	}
	
	// remove unused DOF's
	
}
