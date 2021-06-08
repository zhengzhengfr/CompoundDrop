// Numerical codes for Water-in-oil compound drop impact on a solid surface (Basilisk should be installed before running these codes)
// Zheng Zheng, June 2021
// 0constants-cs.h file, associated with 0impactCompoundSolid-water.c and 0two-phaseDOD.h

#include "axi.h"
#include "navier-stokes/centered.h"
#include "0two-phaseDOD.h"
#include "tension.h"
#include "curvature.h"
#include "tag.h" 

// Water, oil and air properties
#define VELOCITY                        2.4 //impact velocity in m/s
#define DROP_RADIUS             	1.295e-3 //m, outer drop radius
#define DROP_DIAMETER                   (2.*DROP_RADIUS)//m, outer drop diameter
#define DROP_RADIUS_INNER             	0.8695e-3 //m, inner drop radius
#define DROP_DIAMETER_I                 (2.*DROP_RADIUS_INNER)//m, inner drop diameter
#define GAP				0 //m, negative or positive

// Fluid porperties
#define RHO_W				998.0 // water in kg/m3
#define RHO_O				913.0 // oil in kg/m3
#define RHO_A                           1.21  // air in kg/m3

#define MU_W				0.89e-3 	//Pa.s
#define MU_O				4.57e-3    	//Pa.s
#define MU_A                            0.01837e-3 	//Pa.s

#define SIGMA_OA                        2.0e-2 //oil-air  N/m
#define SIGMA_WO                        4.2e-2 //water-oil N/m
#define RATIO_SIGMA			SIGMA_OA/SIGMA_WO //=0.4762            

// non-dimensional numbers
#define REYNOLDS			(RHO_W*VELOCITY*DROP_DIAMETER/MU_W) //based on the fluid properties of water and diameter of oil
#define WEBER				(RHO_W*VELOCITY*VELOCITY*DROP_DIAMETER/SIGMA_WO) //based on rho water, diameter of oil and sigma w-o

// simulation set-up
#define INITAL_GRID_LEVEL		7   //7	
#define MAX_GRID_LEVEL		 	12  //maximum grid level
#define MAX_GEO_LEVEL		 	14  //maximum grid level just before and after the jetting 
#define JET_TIME			6.0 //dimensionless time when jet occurs (estimated and it's quite accurate)
#define MAX_TIME			10.0 //simulation time (non-dimensional)
#define SAVE_FILE_EVERY			0.005 //dump file output every 0.005

// geometry 
#define DOMAIN_WIDTH			10e-3 	 //m, 3.8 times the diameter
#define INITIAL_DISTANCE		0.4e-3   //m, 0.154 times the diameter
#define REFINE_GAP			0.048e-3 //m, 0.0185 times the diameter
#define SHOWITERATION			1

// refine fields VARIABLES
#define REFINE_VAR			{f1, u.x, u.y, f2} // f1, f2 two tracers and velocity fields
#define REFINE_VAR_TEXT		  	"f1, u.x, u.y, f2" // not necessary 
#define REFINE_VALUE_0			-6
#define REFINE_VALUE_1			-3
#define REFINE_VALUE_2			-3
#define REFINE_VALUE_3			-6

#define FILENAME_LASTFILE		"lastfile"

//removal of small bubbles/droplets
#define REMOVE_BUBBLE_YESNO		'y'
#define REMOVE_BUBBLE_SIZE		320.0 // equivalent diameter based on the maximum refinement, 5 is ok for grid (7,9) for the big bubble in the bottom
#define REMOVE_BUBBLE_PERIOD		10

// Other constant values
#define R_VOFLIMIT			1.0e-6 
#define R_PI				3.1415926535897932384626433832795

int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL, LEVELgeo = MAX_GEO_LEVEL;
double Reynolds = REYNOLDS, Weber = WEBER;
double Velocity = VELOCITY;
scalar pressure[], water[], oil[], air[], density[], visco[], area[];
double maxruntime = HUGE; // not necessary 

// These are the parameters/values that will be non-dimensionalized and used in the simulation! 
struct CFDValues {
	double rhoW, rhoO, rhoA, rhoavg, muW, muO, muA, Sigma_WO, Sigma_OA, vel;
	double Reynolds, Weber;
	double radius, diameter, domainsize, refinegap, initialdis, gap;
	double radiusi, diameteri;
	double tau;
	double contacttime;
};

int numericalmainvalues(struct CFDValues *bvalues)
{
	bvalues->rhoW = 1.0; 
	bvalues->rhoO = (RHO_O/RHO_W);
	bvalues->rhoA = (RHO_A/RHO_W);

	bvalues->vel = 1.0;
	bvalues->tau = DROP_DIAMETER/VELOCITY;

	bvalues->radius = 0.5; // = DROP_RADIUS/DROP_DIAMETER;
	bvalues->diameter = 1.0; // = 2.0 * bvalues->radius;
	bvalues->radiusi = DROP_RADIUS_INNER/DROP_DIAMETER;
	bvalues->diameteri = 2.0 * bvalues->radiusi;
	
	bvalues->Reynolds = REYNOLDS;
	bvalues->Weber = WEBER;

	bvalues->muW = (1.0 / bvalues->Reynolds);
	bvalues->muO = (MU_O / MU_W) * bvalues->muW; 
	bvalues->muA = (MU_A / MU_W) * bvalues->muW; 

	bvalues->Sigma_WO = (1.0 / bvalues->Weber);
	bvalues->Sigma_OA = (RATIO_SIGMA) * bvalues->Sigma_WO;
	
	// All length based on the diameter of oil (outer)
	bvalues->domainsize = DOMAIN_WIDTH / DROP_DIAMETER;
	bvalues->initialdis = INITIAL_DISTANCE / DROP_DIAMETER;
	bvalues->refinegap  = REFINE_GAP / DROP_DIAMETER; 
	bvalues->gap = GAP / DROP_DIAMETER; 

	bvalues->contacttime = bvalues->initialdis / bvalues->vel;
	return 1;
}

int runinfoname(char *name, double t)
{
	char tmp[500];
	if (t > -100.0)
	{
		sprintf(tmp, "_t%.3f", t);
		strcat(name, tmp);
	}
	return true;
}

