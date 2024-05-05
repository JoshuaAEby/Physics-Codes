//===========================================================================
//------------------------------Constants------------------------------------
//===========================================================================

//----------------------------- Misc. Constants -----------------------------
#define eq			sqrt(4*M_PI/137)       	// electromagnetic coupling
#define mn			0.9395656		  		// neutron mass in GeV
#define jchi		0.5					  	// DM spin
#define dp			2.8					  	// proton magneton
#define dn			-1.91				  	// neutron magneton

//-------------------------- Escape velocity, can be +-50 -------------------
#define vEsc         550.0               	// Central value of v_esc       [km / sec] 

//--------------------------- NEVER CHANGE, used throughout -----------------
#define cLight 		(3.0E5)             	// Speed of Light			[km/sec]
#define cSqu		pow(cLight,2.0)		  	// c^2				[km^2 / sec^2]
#define hbarc       1.97e-19              	// hbar*c in units km*GeV

#define R_Earth 	6371.0              	// Radius of Earth 			[km]
#define R_mantle    (R_Earth - 30.0)      	// Crust - Mantle boundary
#define R_core      3483.0                  // Mantle - Core boundary

#define v0 		    220.0
#define f_norm      (M_PI * v0 * v0 * (sqrt(M_PI)*v0*erf(vEsc/v0) - 2*vEsc*exp(-pow(vEsc/v0,2))))

//Earth details, updated 5/4/18
#define FE_crust 	1.97e4			  // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]
#define FE_mantle 	3.07e4 		  	  // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]
#define FE_core 	1.01e6 		  	  // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]

#define Pb_crust 	.843	          // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]
#define Pb_mantle 	.0243 		  	  // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]
#define Pb_core 	.127 		  	  // NEEDS TO BE MULTIPLIED BY dens_mult	[km^-3]

#define Si_crust     1.68e5           // NEEDS TO BE MULTIPLIED BY dens_mult    [km^-3]
#define Si_mantle    2.06e5           // NEEDS TO BE MULTIPLIED BY dens_mult    [km^-3]
#define Si_core      1.41e5           // NEEDS TO BE MULTIPLIED BY dens_mult    [km^-3]

//------------------ ROTATIONS ----------------------------------------------
#define a_GP		degrees_to_radians(192.85948)
#define d_GP		degrees_to_radians(27.12825)
#define l_CP		degrees_to_radians(122.932)

//------------------ DETECTOR STUFF -----------------------------------------
#define r_Bor 		.0100881
#define pos_Bor 	(R_Earth - 1.4)
#define theta_Bor 	47.3333		//DEGREES, latitute 42.6

#define r_CJPL 		.0100881
#define pos_CJPL 	(R_Earth - 2.4)
#define theta_CJPL 	61.8468		//DEGREES

#define r_INO 		.0100881
#define pos_INO 	(R_Earth - 1.3)
#define theta_INO 	79.989		//DEGREES

#define r_ANDES 	.0100881
#define pos_ANDES 	(R_Earth - 1.75)
#define theta_ANDES 114.611		//DEGREES

#define r_SUPL 		.0100881
#define pos_SUPL 	(R_Earth - 1.025)
#define theta_SUPL 	127.05		//DEGREES, latitute -37.05

#define r_JUNO       .0177000
#define pos_JUNO     (R_Earth - .48125)
#define theta_JUNO   67.892  	//DEGREES, latitute 22.1

#define r_Boulby 		.0100881
#define pos_Boulby 	(R_Earth - 1.1)
#define theta_Boulby 35.67		//DEGREES, latitute 54.33

//------------------ LARGE MULTIPLIERS, put in at the end -------------------
#define sigma_0 	(1.0E-50)		  	  // DM-Nucleon Cross Section		[km^2]
#define DM_massdensity  (3.0E14 / cSqu)
#define dens_mult 	(1.0E32)
