// Numerical codes for Water-in-oil compound drop impact on a solid surface (Basilisk should be installed before running these codes)
// Zheng Zheng, June 2021
// Post-processing code to input the resulting dump files and to generate images at each time step in order to make a video

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
	double tbegin = 0.000;
	double tend = 5.555;
	double t;
	double tstep = 0.005;
	char figurename[100];
	int i = 0;
	for(t = tbegin; t <= tend; t += tstep)
	{ 
		sprintf(filename,"data_t%.3f", t);
		restore(file = filename);
		sprintf(figurename,"vof_%03d.png", i); 
		view(width = 2000, height = 2000, fov = 19, ty = -0.5,quat = {0, 0, -0.707, 0.707});
		draw_vof ("f1", edges = true);
		draw_vof ("f2", edges = true);
		squares ("area", spread = -1);
    		mirror ({0,1}) {
		    	draw_vof ("f1", edges = true);
		    	draw_vof ("f2", edges = true);
				squares ("area", spread = -1);
	  		}
		save(figurename);
		i++;
	}
}
event end(i = 0)
{
}

/*
command line to compile for 1 core run
qcc -O2 -Wall images.c -o images -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
./images

//from images to video
ffmpeg -r 8 -i vof_%03d.png name.mp4
*/
