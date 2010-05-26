#ifndef __H1_PERIODIC_SPACE
#define __H1_PERIODIC_SPACE

#include <map>
#include <set>
#include "hermes2d.h"

using namespace std;

/// Edges assigned to periodic BCs form pairs with markers >= 1001 (lesser values
/// are used by regular BCs). One edge from each pair serves as a source edge 
/// (marker 1xxx), the other represents its image (marker 2xxx).
/// For each function in the periodic space, its values over the image edge must match
/// those over the corresponding source edge. This periodicity is enforced 
/// (as an essential condition) by assigning the same DOFs to nodes belonging
/// to the source and image edges. Equality of normal derivatives of the periodic	
/// basis functions is satisfied naturally by the weak form.

const int SRC_MARKERS_START = 1001;
const int MARKERS_DIFF = 1000;
const int IMG_MARKERS_START = 2001;

#define for_all_bnd_nodes(n, mesh) \
				for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
				  if (((n) = (mesh)->get_node(_id))->used) \
				    if ((n)->bnd)
				    
typedef multimap<int,Node*> mnmap;

class H1PeriodicSpace : public H1Space
{
	public:
		H1PeriodicSpace(Mesh* mesh = NULL, Shapeset* shapeset = NULL) : H1Space(mesh, shapeset) {};
		
		// source nodes DOFs must be known during the assigning of DOFs for the periodic basis,
		// so it should be done after the regular assignment
		void post_assign();
		
		inline bool is_nat_bnd(Node* nd) { return (nd->bnd && ndata[nd->id].dof >= 0); }
};

#endif
