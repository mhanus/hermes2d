#include "PNGGraph.h"

static void get_style_types(std::string line, std::string mark, std::string col, int& lt, int& pt, int& ct)
{
  if      (line == "-")  lt = 1; // solid
  else if (line == ":")  lt = 4; // dotted
  else if (line == "-.") lt = 5; // dash dot
  else if (line == "--") lt = 2; // dashed
  else lt = 1;

  if      (mark == ".") pt = 7;  // full circle
  else if (mark == "o") pt = 6;  // empty circle
  else if (mark == "O") pt = 7;  // full circle
  else if (mark == "x") pt = 2;  // cross
  else if (mark == "+") pt = 1;  // cross
  else if (mark == "*") pt = 3;  // star
  else if (mark == "s") pt = 4;  // empty square
  else if (mark == "S") pt = 5;  // full square
  else if (mark == "d") pt = 10; // empty diamond
  else if (mark == "D") pt = 11; // full diamond
  else if (mark == "v") pt = 12; // empty triangle down
  else if (mark == "V") pt = 13; // full triangle down
  else if (mark == "^") pt = 9;  // full triangle up
  else if (mark == "<") pt = 12; // empty triangle down
  else if (mark == ">") pt = 8;  // empty triangle up
  else if (mark == "p") pt = 14; // empty pentagon
  else if (mark == "P") pt = 15; // full pentagon
  else pt = 0;

  if      (col == "k") ct = -1;  // black
  else if (col == "b") ct = 3;   // blue
  else if (col == "g") ct = 2;   // green
  else if (col == "c") ct = 5;   // cyan
  else if (col == "m") ct = 4;   // magenta
  else if (col == "y") ct = 6;   // yellow
  else if (col == "r") ct = 1;   // red
  else ct = -1;
}

void PNGGraph::save(const char* filename)
{
  int j;

  if (!rows.size()) error("No data rows defined.");

  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Error writing to %s", filename);

	if (trans)
	  fprintf(f, " set terminal png font arial 14 transparent truecolor size %g,%g crop\n", this->w, this->h);
	else
		fprintf(f, " set terminal png font arial 14 size %g,%g crop\n", this->w, this->h);
		
  int len = strlen(filename);
  AUTOLA_OR(char, outname, len + 10);
  strcpy(outname, filename);
  char* slash = strrchr(outname, '/');
  if (slash != NULL) strcpy(outname, ++slash);
  char* dot = strrchr(outname, '.');
  if (dot != NULL && dot > outname) *dot = 0;
  strcat(outname, ".png");

  fprintf(f, "set output '%s'\n", (char*)outname);

  if (logx && !logy)
    fprintf(f, "set logscale x\n");
  else if (!logx && logy)
    fprintf(f, "set logscale y\n");
  else if (logx && logy)
  {
    fprintf(f, "set logscale x\n");
    fprintf(f, "set logscale y\n");
  }

  if (grid) fprintf(f, "set grid\n");

  if (title.length()) fprintf(f, "set title '%s'\n", title.c_str());
  if (xname.length()) fprintf(f, "set xlabel '%s'\n", xname.c_str());
  if (yname.length()) fprintf(f, "set ylabel '%s'\n", yname.c_str());
  if (legend && legend_pos.length()) fprintf(f, "set key %s\n", legend_pos.c_str());

  fprintf(f, "plot");
  for (unsigned int i = 0; i < rows.size(); i++)
  {
    int ct, lt, pt;
    get_style_types(rows[i].line, rows[i].marker, rows[i].color, lt, pt, ct);

    if (lt == 0)
      fprintf(f, " '-' w p pointtype %d", pt);
    else if (ct < 0)
      fprintf(f, " '-' w lp linewidth %g linetype %d pointtype %d", this->plw, lt, pt);
    else
      fprintf(f, " '-' w lp linewidth %g linecolor %d linetype %d pointtype %d", this->plw, ct, lt, pt);

		if (legend)
			fprintf(f, " title '%s' ", rows[i].name.c_str());
		else
			fprintf(f, " notitle ", rows[i].name.c_str());
			
    if (i < rows.size() - 1) fprintf(f, ", ");
  }
  fprintf(f,"\n");

  for (unsigned int i = 0; i < rows.size(); i++)
  {
    int rsize = rows[i].data.size();
    for (j = 0; j < rsize; j++)
      fprintf(f, "%.14g  %.14g\n", rows[i].data[j].x, rows[i].data[j].y);
    fprintf(f, "e\n");
  }

  fprintf(f, "set terminal x11\n");
  fclose(f);

  verbose("Graph saved. Process '%s' via gnuplot.", filename);
}
