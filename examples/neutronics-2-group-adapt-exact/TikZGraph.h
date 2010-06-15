#ifndef __TIKZ_GRAPH
#define __TIKZ_GRAPH

#include "hermes2d.h"

class TikZGraph : public Graph
{
	public:
		TikZGraph(const char* title = NULL, const char* x_axis_name = NULL, const char* y_axis_name = NULL, 
							double plot_width = 9.5, double plot_height = 6.5, const double plot_lines_width = 2)
		     : Graph(title, x_axis_name, y_axis_name) { plw = plot_lines_width; w = plot_width; h = plot_height; legend_pos = ""; }

		virtual void save(const char* filename);
		void set_legend_pos(const char* posspec) { if (!legend) legend = true; legend_pos = posspec; }
		
	private:
		double w, h;	// plot size
		double plw;		// line widths for graphs
		std::string legend_pos;	// legend position (see documentation to gnuplot)
};

#endif
