// Numerical codes for Water-in-oil compound drop impact on a solid surface (Basilisk should be installed before running these codes)
// Zheng Zheng, June 2021
// 0impactCompoundSolid-water.c file, associated with 0constants-cs.h and 0two-phaseDOD.h

//#include "view.h"  //not compatible with parallel run for movies
//#include "draw.h"  //not compatible with parallel run for images
#include "0constants-cs.h"

#if dimension == 3
#include "lambda2.h"
#endif

static void remove_droplets1(scalar c, bool droplet, double dropsize);
struct CFDValues cfdbv;

int main(int argc, char **argv)
{
	numericalmainvalues(&cfdbv);
	size(cfdbv.domainsize); 
#if AXI
#else
	origin(0, -cfdbv.domainsize / 2., -cfdbv.domainsize / 2.);
#endif
	int initialgrid = pow(2, LEVELmin);   
	init_grid(initialgrid);    		
	//density
	rho1 = cfdbv.rhoW; //water
	rho2 = cfdbv.rhoO; //oil
	rho3 = cfdbv.rhoA; //air 
	//viscosity
	mu1 = cfdbv.muW;   //water 
	mu2 = cfdbv.muO;   //oil
	mu3 = cfdbv.muA;   //air
	// 2 tracers
	f1.sigma = cfdbv.Sigma_WO;  //water-oil interface
	f2.sigma = cfdbv.Sigma_OA;  //oil-air interface

	//Boundary conditions imposed 
	f1[left] = 0.;    //BC for water, 0 for 180 degree, 1 for 0 degree
	f2[left] = 1.;    //BC for oil, 0 for 180 degree, 1 for 0 degree
	
	//No-slip boundary condition: both the velocity normal to the boundary and the velocity parallel to the boundary are set equal to zero.
	u.t[left] = dirichlet(0);  //In Basilisk, left boundary is the solid wall 
	u.n[left] = dirichlet (0); 
	// Right and top boundaries 
	u.n[right] = neumann(0);
	p[right] = dirichlet(0);
	pf[right] = dirichlet(0);
	u.n[top] = neumann(0);
	p[top] = dirichlet(0);
	pf[top] = dirichlet(0);

	TOLERANCE = 1e-6;
	run();
	return 1;
}

event init(i = 0)
{
	if (restore (file = FILENAME_LASTFILE))
	{
#if AXI
		boundary((scalar *){fm}); // boundary bottom symmetric condition  
#endif
	}
	else
	{
		double x0 = cfdbv.initialdis + cfdbv.radius; //center for compound drop
		double x1 = cfdbv.initialdis + cfdbv.radius + cfdbv.gap; //center for water core	
		refine(sq(x - x0) + sq(y) < sq(cfdbv.radius + cfdbv.refinegap) && sq(x - x0) + sq(y) > sq(cfdbv.radius - cfdbv.refinegap) && level < LEVELmax);
		fraction(f2, sq(cfdbv.radius) - sq(x - x0) - sq(y)); //tracer for oil-air
		refine(sq(x - x1) + sq(y) < sq(cfdbv.radiusi + cfdbv.refinegap) && sq(x - x1) + sq(y) > sq(cfdbv.radiusi - cfdbv.refinegap) && level < LEVELmax);
		fraction(f1, sq(cfdbv.radiusi) - sq(x - x1) - sq(y)); //tracer for water-oil
		foreach()
		{
		  if (sq(x - x0) + sq(y) < sq(cfdbv.radius))
		     {
		     u.x[] = -cfdbv.vel; //Initial velocity, imposed 
		     }
		}
	}
}

//Removal of small bubbles 
#if REMOVE_BUBBLE_YESNO == 'y'
event drop_remove(t = 0; t += 0.001; t <= JET_TIME - 0.4) //only during a period, before jetting 
{
	if (t > cfdbv.contacttime)
	{
		remove_droplets1(f2, 0, REMOVE_BUBBLE_SIZE); //removal of air bubble in the oil shell
	}
}
#endif

//Adaptation of the mesh according to the error on the volume fraction and velocity field
event adapt(i++)
{
	double refine[3];
	refine[0] = pow(10.0, REFINE_VALUE_0);
	refine[1] = pow(10.0, REFINE_VALUE_1);
	refine[2] = pow(10.0, REFINE_VALUE_2);
	refine[3] = pow(10.0, REFINE_VALUE_3);
	if (t <= JET_TIME - 0.4) 
	{
	adapt_wavelet(REFINE_VAR, (double[]) {refine[0],refine[1],refine[2],refine[3]}, maxlevel = LEVELmax, minlevel = LEVELmin); 
	}
	if (t > JET_TIME - 0.4 && t <= JET_TIME + 1.0) //increased mesh level for the jet ejection (just before and after)
	{
	adapt_wavelet(REFINE_VAR, (double[]) {refine[0],refine[1],refine[2],refine[3]}, maxlevel = LEVELgeo, minlevel = LEVELmin); 
	}
	if (t > JET_TIME + 1.0 && t <= 10.0)
	{
	adapt_wavelet(REFINE_VAR, (double[]) {refine[0],refine[1],refine[2],refine[3]}, maxlevel = LEVELmax, minlevel = LEVELmin); 
	}
}

event showiteration(i += SHOWITERATION)
{
	switch (pid()) 
	{
	case 0:
	{
		char name[500], tmp[500];
		sprintf(name, "i = %05d, dt = %.3e, t = %.3f, P =  %03d", i, dt, t, (int)(100.0 * t / MAX_TIME));
		sprintf(tmp,", Re = %.2f, We = %.2f", (double)cfdbv.Reynolds, (double)cfdbv.Weber);
		strcat(name, tmp);
		runinfoname(name, t);
		printf("%s\r\n", name);
	}
	}
}

event end(t = MAX_TIME)
{
;
}

event outputfiles(t += SAVE_FILE_EVERY)
{
	char name[500];
	// write data dump file
	foreach(){
		pressure[] = p[]; //pressure field
		water[] = f1[]*f2[]; //to see water phase in Bview
		oil[] = f2[]*(1 - f1[]); //to see oil phase in Bview
		air[] = 1 - f2[];  //to see air phase in Bview
		density[] = water[]*rho1 + oil[]*rho2 + air[]*rho3; //to see density distribution in Bview
		visco[] = water[]*mu1 + oil[]*mu2 + air[]*mu3; //to see viscosity distribution in Bview
		area[] = water[]*-1 + oil[]*1 + air[]*0;  //to change the color of water-oil-air phase in Bview
	}
	sprintf(name, "data");
	runinfoname(name, t);
	dump(file = name);
	dump(file = FILENAME_LASTFILE);
}

//Function to remove small bubbles or droplets 
static void remove_droplets1(scalar c, bool droplet, double dropsize)
{
	int j, n;
	scalar m[];
	const double THR = R_VOFLIMIT; //THRESHOLD
	const double delta = cfdbv.domainsize / pow(2.0, LEVELmax);
	double realD2;
	foreach ()
		m[] = (droplet ? (c[] > THR) : (c[] < (1. - THR)));
	n = tag(m);
	double v[n];
	for (j = 0; j < n; j++)
		v[j] = 0.0;
	foreach_leaf()
	{
		if (m[] > 0)
		{
			j = m[] - 1;
#if dimension == 3
			v[j] += Delta * Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#else
			v[j] += Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#endif
		}
	}
#if _MPI
	MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	foreach ()
	{
		if (m[] > 0)
		{
			j = m[] - 1;
			realD2 = 4.0 * v[j] / R_PI;
			if (realD2 < dropsize * dropsize * delta * delta)
				c[] = droplet ? 0. : 1.;
		}
	}
	boundary({c});
}

//We generate a movie of the interface using Basilisk View, only for non parallel run
/*event movie (t += 0.01)
{
  //view (fov = 13.1165, quat = {0,0,0,1}, tx = 0, ty = 0, bg = {1,1,1}, width = 1000, height = 1000, samples = 1);
  view (width = 500, height = 500, fov = 20, ty = -0.5,quat = {0, 0, -0.707, 0.707});
  clear();
  draw_vof ("f1", lw=2);
  draw_vof ("f2", lw=2);
  box(notics = true);
  mirror ({0,1}) {
    draw_vof ("f1", lw = 2);
    draw_vof ("f2", lw = 2);
    box (notics = true);
  }
  save ("final-v2.5.mp4");
  //We can optionally visualise the results while we run.
  //static FILE * fp = popen ("bppm","w");
  //save (fp = fp);
}*/


/* Command lines to compile and to run with multiple cores (specify core number)
qcc -source -D_MPI=1 0impactCompoundSolid-water.c 
mpicc -O2 -Wall -std=c99 -D_MPI=1 _0impactCompoundSolid-water.c -o cs -lm
mpirun -np 27 ./cs

//command to compile for 1 core 
qcc -O2 -Wall 0impactCompoundSolid-water.c -o cs -lm -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
./cs

*/
 
