// Numerical codes for Water-in-oil compound drop impact on a solid surface (Basilisk should be installed before running these codes)
// Zheng Zheng, June 2021
// Post-processing code to calculate air bubble volume and oil volume 

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
	double tbegin = 8.450;
	double tend = 8.470;
	double t;
	double tstep = 0.01;
	char figurename[100];
	int i = 0; 
	double volume_dim;
	for(t = tbegin; t <= tend; t += tstep)
	{ 
		double volume = 0.;
		sprintf(filename,"data_t%.3f", t);
		restore(file = filename);
		sprintf(figurename, "volume-%.3f.txt", t);
		FILE* fp = fopen (figurename, "w");
		foreach ()
		{
  		if (x > 0.02451 && y > 0.336 && y < 0.5884) //zone of air bubble or oil 
		   {
			volume += 2*R_PI*dv()*(1 - f2[]); //air bublble
			//volume += 2*R_PI*dv()*(f2[]*(1 - f1[])); //oil volume
		   }
		}
		volume_dim = volume*(2.59*2.59*2.59);//dimensional in mm^3
		fprintf(fp, "%.10f \n", volume_dim);
		fclose(fp);
		i++;
	}
}

event end(i = 0)
{
}

/*
command to compile for 1 core 
qcc -O2 -Wall volume_bubble.c -o volume_bubble -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa 
./volume_bubble
*/
