// Numerical codes for Water-in-oil compound drop impact on a solid surface (Basilisk should be installed before running these codes)
// Zheng Zheng, June 2021
// Post-processing code to input the resulting dump files and to output all the interface points for post-processing

#include "view.h"  
#include "0constants-cs.h"

int main()
{
	run();
	return 1;
}

event init(i = 0)
{
	char filename[100];
	double tbegin = 5.550;
	double tend = 5.558;
	double t;
	double tstep = 0.001;
	char figurename[100];
	int i = 0;
	for(t = tbegin; t <= tend; t += tstep)
	{ 
		sprintf(filename,"data_t%.3f", t);
		restore(file = filename);
		sprintf(figurename,"interface-%.3f.txt", t);
		FILE* fp = fopen (figurename, "w");
  		output_facets (f1, fp); //f1 for water-oil; f2 for oil-air interface
		i++;
	}
}

event end(i=0)
{
}

/*
command to compile for 1 core 
qcc -O2 -Wall output_position.c -o output_position -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
./output_position
*/
