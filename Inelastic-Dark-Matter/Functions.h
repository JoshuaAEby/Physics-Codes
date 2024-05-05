
// //===========================================================================================
//---------------- Define Vectors -----------------------------------------------------------
//===========================================================================================

struct param_type 
{ 
    double year, month, day, hour;
    char detector, element[5];
    long double mDM, delta;
	double jN, lifemult;
    int Znum, Anum; 
	long double Dcore, Dmantle, Dcrust;
	double Wcoeffs[6][2][2][10];
	double logprec, gM;
};

//==========================================================================================
//----------------------Functions on Vectors------------------------------------------------
//=========================================================================================

// Convert degrees to radians
double degrees_to_radians(double angle)
{
    return angle * M_PI / 180.;
}

// Magnitude of a vector in cartesian coordinates
double vector_magnitude(double vec[3])
{
    return sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));
}

// Binary operations on two vectors:
    // op == 'a' means addition
    // op == 's' means subtraction
    // op == 'm' means multiplication
    // op == 'd' means division
double * vector_operation(double vec_1[3], double vec_2[3],char op)
{
    static double vec_out[3];
    if(op == 'a' || op == 'A')
    {   vec_out[0] = vec_1[0] + vec_2[0];
        vec_out[1] = vec_1[1] + vec_2[1];
        vec_out[2] = vec_1[2] + vec_2[2];
        return vec_out;
    }
    else if(op == 's' || op == 'S')
    {   vec_out[0] = vec_1[0] - vec_2[0];
        vec_out[1] = vec_1[1] - vec_2[1];
        vec_out[2] = vec_1[2] - vec_2[2];
        return vec_out;
    }
    else if(op == 'm' || op == 'M')
    {   vec_out[0] = vec_1[0] * vec_2[0];
        vec_out[1] = vec_1[1] * vec_2[1];
        vec_out[2] = vec_1[2] * vec_2[2];
        return vec_out;
    }
    else
    {
        printf("Not a known operation!\n");
        return 0;
    }
}

/*
// Add two vectors in cartesian coordinates
double * vector_addition(double vec_1[3], double vec_2[3])
{
    static double vec_out[3];
    vec_out[0] = vec_1[0] + vec_2[0];
    vec_out[1] = vec_1[1] + vec_2[1];
    vec_out[2] = vec_1[2] + vec_2[2];
    return vec_out;
}

// Subtract two vectors in cartesian coordinates
double * vector_subtraction(double vec_1[3], double vec_2[3])
{
    static double vec_out[3];
    vec_out[0] = vec_1[0] - vec_2[0];
    vec_out[1] = vec_1[1] - vec_2[1];
    vec_out[2] = vec_1[2] - vec_2[2];
    return vec_out;
}
*/

double * rotate_cartesian(double vec_in[3], double rot_angle, char axis)
{
    static double vec_out[3];
    
    double cs = cos(rot_angle);
    double sn = sin(rot_angle);
    
    if(axis == 'x')
    {
        vec_out[0] = vec_in[0];
        vec_out[1] = cs * vec_in[1] - sn * vec_in[2];
        vec_out[2] = sn * vec_in[1] + cs * vec_in[2];
    }
    else if(axis == 'y')
    {
        vec_out[0] = cs * vec_in[0] + sn * vec_in[2];
        vec_out[1] = vec_in[1];
        vec_out[2] = -sn * vec_in[0] + cs * vec_in[2];
    }
    else if(axis == 'z')
    {
        vec_out[0] = cs * vec_in[0] - sn * vec_in[1];
        vec_out[1] = sn * vec_in[0] + cs * vec_in[1];
        vec_out[2] = vec_in[2];
    }
    else
        printf("What axis is that??");
    
    return vec_out;
}

// Rotate a cartesian vector around its x axis
double * rotate_cartesian_x(double vec_in[3], double rot_angle)
{
    static double vec_out[3];

    double cs = cos(rot_angle);
    double sn = sin(rot_angle);

    vec_out[0] = vec_in[0];
    vec_out[1] = cs * vec_in[1] - sn * vec_in[2];
    vec_out[2] = sn * vec_in[1] + cs * vec_in[2];

    return vec_out;
}

// Rotate a cartesian vector around its y axis
double * rotate_cartesian_y(double vec_in[3], double rot_angle)
{
    static double vec_out[3];

    double cs = cos(rot_angle);
    double sn = sin(rot_angle);

    vec_out[0] = cs * vec_in[0] + sn * vec_in[2];
    vec_out[1] = vec_in[1];
    vec_out[2] = -sn * vec_in[0] + cs * vec_in[2];

    return vec_out;
}

// Rotate a cartesian vector around its z axis
double * rotate_cartesian_z(double vec_in[3], double rot_angle)
{
    static double vec_out[3];

    double cs = cos(rot_angle);
    double sn = sin(rot_angle);

    vec_out[0] = cs * vec_in[0] - sn * vec_in[1];
    vec_out[1] = sn * vec_in[0] + cs * vec_in[1];
    vec_out[2] = vec_in[2];

    return vec_out;
}

double * cart_to_sph(double vec_Cart[3]) // (x,y,z) to (r, cos_theta, phi)
{
    static double vec_Sph[3];
    double mag = vector_magnitude(vec_Cart);

    vec_Sph[0] = mag;
    vec_Sph[1] = vec_Cart[2] / mag;
    vec_Sph[2] = atan2(vec_Cart[1] , vec_Cart[0]);

    return vec_Sph;
} 

double * sph_to_cart(double vec_Sph[3]) // (r, cos_theta, phi) to (x,y,z)
{
    static double vec_Cart[3];

  //CHECK THIS, TEST: Why don't these different calculations agree? They do.
  //vec_Cart[0] = vec_Sph[0] * sin(acos(vec_Sph[1])) * cos(vec_Sph[2]);
  //vec_Cart[1] = vec_Sph[0] * sin(acos(vec_Sph[1])) * sin(vec_Sph[2]);  
    vec_Cart[0] = vec_Sph[0] * sqrt(1.0 - pow(vec_Sph[1],2)) * cos(vec_Sph[2]);
    vec_Cart[1] = vec_Sph[0] * sqrt(1.0 - pow(vec_Sph[1],2)) * sin(vec_Sph[2]);
    vec_Cart[2] = vec_Sph[0] * vec_Sph[1];

//  printf("sin(acos(theta)) = %f, sqrt(1 - theta^2) = %f\n", sin(acos(vec_Sph[1])), sqrt(1.0 - pow(vec_Sph[1],2)));

    return vec_Cart;
} 

double * cart_to_sph2(double vec_Cart[3]) 			// (x,y,z) to (r, cos_theta, phi)
{
    static double vec_Sph[3];
    double mag = vector_magnitude(vec_Cart);

    vec_Sph[0] = mag;
    vec_Sph[1] = acos(vec_Cart[2] / mag);
    vec_Sph[2] = atan2(vec_Cart[1] , vec_Cart[0]);

    return vec_Sph;
} 

double * sph_to_cart2(double vec_Sph[3]) 			// (r, theta, phi) to (x,y,z)
{
    static double vec_Cart[3];

    vec_Cart[0] = vec_Sph[0] * sin(vec_Sph[1]) * cos(vec_Sph[2]);
    vec_Cart[1] = vec_Sph[0] * sin(vec_Sph[1]) * sin(vec_Sph[2]);
    vec_Cart[2] = vec_Sph[0] * cos(vec_Sph[1]);

    return vec_Cart;
} 

// Cosine of the angle between two vectors in Cartesian coordinates
double cos_angle_between(double v1[3], double v2[3])
{
    double cos_angle;
    double mag1 = vector_magnitude(v1);
    double mag2 = vector_magnitude(v2);

    cos_angle = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    cos_angle = cos_angle / mag1 / mag2;

    return cos_angle;
}

//=========================================================================
//-------------Display Functions-------------------------------------------
//=========================================================================

void print ( double x ) 
{
    printf("%f\n", x);
}
void display_results (char *title, double result, double error)
{
    printf ("          ~~~~~~~~~~~~~~~~~~~~~ %s ~~~~~~~~~~~~~~~~~~~~~~\n", title);
    printf ("          Result           = % .3e\n", result);
    printf ("          Method Error     = % .3e\n", error);
    printf ("          Fractional Error = % .3f percent\n", error/result*100);
    printf ("          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

void detector_banner(char det_letter) 
{
	printf("\n\n");
  printf("==========================================================================\n");
    if( det_letter == 'B') {
      printf("---------------------------- Detector: Borexino --------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees ---------------\n",-pos_Bor+R_Earth,90-theta_Bor);
    }
    else if ( det_letter == 'C') {
      printf("----------------------------- Detector: CJPL ----------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees ---------------\n",-pos_CJPL+R_Earth,90-theta_CJPL);
    }
    else if ( det_letter == 'I') {
      printf("------------------------------ Detector: INO ----------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees ---------------\n",-pos_INO+R_Earth,90-theta_INO);
    }
    else if ( det_letter == 'S') {
      printf("------------------------------ Detector: SUPL ---------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees --------------\n",-pos_SUPL+R_Earth,90-theta_SUPL);
    }
    else if ( det_letter == 'A') {
      printf("----------------------------- Detector: ANDES ---------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees --------------\n",-pos_ANDES+R_Earth,90-theta_ANDES);
    }
    else if ( det_letter == 'J') {
      printf("----------------------------- Detector: JUNO ----------------------------\n");
      printf("---------------------- Depth %.1f km, Latitude %.1f degrees ---------------\n",-pos_JUNO+R_Earth,90-theta_JUNO);
    }
  printf("==========================================================================\n");
	printf("\n\n");
}

//==========================================================================================
//--------------------- Set Detector Params ------------------------------------------------
//==========================================================================================

void detector_choice (char det_letter, double det_vars_out[3]) 
{
    if (det_letter == 'B') {			    			 // Borexino, should probably call it Gran Sasso
        det_vars_out[0]   	= r_Bor;        		     // Radius of "spherical" detector [km]
        det_vars_out[1]     = pos_Bor;		 	     	 // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_Bor;  			     // Latitude 42.666667 N
    }

    else if (det_letter == 'C') {		  	     		 // CJPL, "JinPing"
        det_vars_out[0]   	= r_CJPL;       		     // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_CJPL;			     	 // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_CJPL;  		     	 // Latitude 28.153232 N
    }

    else if (det_letter == 'I') {			     		 // INO, "India-based Neutrino Observatory"
        det_vars_out[0]   	= r_INO;         		     // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_INO;		 	     	 // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_INO;			     // Latitude 10.010986 N
    }
    else if (det_letter == 'A') {			     		 // ANDES, "ANDES Underground Lab"
        det_vars_out[0]   	= r_ANDES;         		     // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_ANDES;			     // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_ANDES;			     // Latitude -24.6114 S
    }
    else if (det_letter == 'S') {		   	     		 // SUPL, "Stawell Underground Physics Lab"
        det_vars_out[0]   	= r_SUPL;         		     // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_SUPL;			     	 // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_SUPL;       		 // Latitude -37.05 S
    }
    else if (det_letter == 'J') {		   	     		 // JUNO, "Jiangmen Underground Neutrino Observatory"
        det_vars_out[0]   	= r_JUNO;         		     // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_JUNO;			     	 // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_JUNO;       		 // Latitude 22.11 N
    }
    else if (det_letter == 'E') {		   	     		 // England, "Boulby Mine"
        det_vars_out[0]   	= r_Boulby;         		 // Radius of "spherical" detector 	[km]
        det_vars_out[1]     = pos_Boulby;			     // Detector distance from Earth center [km]
        det_vars_out[2]    	= theta_Boulby;       		 // Latitude 54.33 N
    }
    return;
}

//==========================================================================================
//---------------------Dynamics of DM Speed-------------------------------------------------
//==========================================================================================

double date_to_J2000(double year, double month, double day, double hour)
{
    double yy, mm, dd;

    if ( month == 1 || month == 2){
		yy = year - 1;
		mm = month + 12;
	}
    else{
		yy = year;
		mm = month;
	}
    dd = day + hour / 24.0;
 
    return floor(365.25 * yy) + floor(30.61 * (mm + 1)) + dd - 730563.5;
}

//------------------------------- Coordinate rotations ----------------------------------------

double * coordinates_gal_to_equat(double vec_galactic[3])
{
    static double vec_out[3];

 //This matrix is from McCabe, arXiv: 1312.1355. It is the inverse of M in (B3) of the Appendix.

    vec_out[0] = 		   vec_galactic[0] * ( - cos(a_GP) * cos(l_CP) * sin(d_GP) - sin(a_GP) * sin(l_CP) );
    vec_out[0] = vec_out[0] + vec_galactic[1] * ( - cos(a_GP) * sin(l_CP) * sin(d_GP) + sin(a_GP) * cos(l_CP));
    vec_out[0] = vec_out[0] + vec_galactic[2] * cos(a_GP) * cos(d_GP);

    vec_out[1] = 		   vec_galactic[0] * ( - sin(a_GP) * cos(l_CP) * sin(d_GP) + cos(a_GP) * sin(l_CP));
    vec_out[1] = vec_out[1] + vec_galactic[1] * ( - sin(a_GP) * sin(l_CP) * sin(d_GP) - cos(a_GP) * cos(l_CP));
    vec_out[1] = vec_out[1] + vec_galactic[2] * sin(a_GP) * cos(d_GP);

    vec_out[2] = 		   vec_galactic[0] * cos(d_GP) * cos(l_CP);
    vec_out[2] = vec_out[2] + vec_galactic[1] * cos(d_GP) * sin(l_CP);
    vec_out[2] = vec_out[2] + vec_galactic[2] * sin(d_GP);

    return vec_out;
}

double * coordinates_equat_to_gal(double vec_equ[3])
{
    static double vec_out[3];

 //This matrix is from McCabe, arXiv: 1312.1355. It is M in (B3) of the Appendix.

    vec_out[0] = 		   vec_equ[0] * ( - cos(a_GP) * cos(l_CP) * sin(d_GP) - sin(a_GP) * sin(l_CP));
    vec_out[0] = vec_out[0] + vec_equ[1] * ( - sin(a_GP) * cos(l_CP) * sin(d_GP) + cos(a_GP) * sin(l_CP));
    vec_out[0] = vec_out[0] + vec_equ[2] * cos(d_GP) * cos(l_CP);

    vec_out[1] = 		   vec_equ[0] * ( - cos(a_GP) * sin(l_CP) * sin(d_GP) + sin(a_GP) * cos(l_CP));
    vec_out[1] = vec_out[1] + vec_equ[1] * ( - sin(a_GP) * sin(l_CP) * sin(d_GP) - cos(a_GP) * cos(l_CP));
    vec_out[1] = vec_out[1] + vec_equ[2] * cos(d_GP) * sin(l_CP);

    vec_out[2] = 		   vec_equ[0] * cos(a_GP) * cos(d_GP);
    vec_out[2] = vec_out[2] + vec_equ[1] * sin(a_GP) * cos(d_GP);
    vec_out[2] = vec_out[2] + vec_equ[2] * sin(d_GP);

    return vec_out;
}

//------------------------------------------- Velocity in Earth Frame ---------------------------------------------------

// Compute Earth velocity at time n_date days from J2000.0 (noon GMT on 1/1/2000)
// Using McCabe, arXiv: 1312.1355
double * compute_Earth_velocity(double year, double month, double day, double hour)
{
    static double v_out[3];

    double n_date	= date_to_J2000(year, month, day, hour);

    double T_date 	= n_date / 36525.0;

    double eps_x[3] = {.054876 - .024232 * T_date, -.494109 -.002689 * T_date, .867666 + .000001546 * T_date};
    double eps_y[3] = {.993821 + .001316 * T_date, .110992 - .011851 * T_date, .000352 + .021267    * T_date};

    double speed_E  = 29.79;												// average Earth speed

    double e	    = degrees_to_radians(.9574);							// orbital eccentricity
    //double g 	 = degrees_to_radians(357.528 + .9856003 * n_date);			// mean anomaly
    double omega	= degrees_to_radians(282.932 + .0000471 * n_date);  	// perihelion longitude
    double L        = degrees_to_radians(280.460 + .9856474 * n_date); 		// mean longitude
 
    //double nu 	 = g + 2 * e * sin(g) + (5/4) * pow(e,2) * sin(2 * g);  // true anomaly
    //double ell	 = omega + nu;											// ecliptic longitude

 // Velocity IN GALACTIC COORDINATES
    v_out[0]        = - speed_E * ( sin(L) + e * sin(2 * L - omega) ) * eps_x[0];
    v_out[0]        = v_out[0] + speed_E * ( cos(L) + e * cos(2 * L - omega) ) * eps_y[0];
    v_out[1]        = - speed_E * ( sin(L) + e * sin(2 * L - omega) ) * eps_x[1];
    v_out[1]        = v_out[1] + speed_E * ( cos(L) + e * cos(2 * L - omega) ) * eps_y[1];
    v_out[2]        = - speed_E * ( sin(L) + e * sin(2 * L - omega) ) * eps_x[2];
    v_out[2]        = v_out[2] + speed_E * ( cos(L) + e * cos(2 * L - omega) ) * eps_y[2];

    return v_out;
}

double * compute_Earth_velocity_TOT(double year, double month, double day, double hour) 
{
    double v_LSR[3] 	= {0.0,  v0,   0.0};								// Local Standard of Rest velocity
    double v_pec[3] 	= {11.1, 12.2, 7.3};								// Sun peculiar velocity
    double *p	= compute_Earth_velocity(year,month,day,hour);				// Earth peculiar velocity
    
    static double u_E[3];
    for (int i = 0; i < 3; i++ ) u_E[i] = *(p+i);
    p			= vector_operation(u_E, v_LSR,'a');//vector_addition(u_E, v_LSR);
    for (int i = 0; i < 3; i++ ) u_E[i] = *(p+i);
    p			= vector_operation(u_E, v_pec,'a');//vector_addition(u_E, v_pec);
    for (int i = 0; i < 3; i++ ) u_E[i] = *(p+i);

    return u_E;
}


//========================================================================================================================
//----------------------------------------------- Old Form Factor Functions ----------------------------------------------
//========================================================================================================================

long double lifetime_h(long double mDM, long double delta)
{
    long double lifetime;
    lifetime = (.0678) * pow(350.0/delta,3) * pow(mDM/1000.0,2);  			//[sec]
    return lifetime;
}
/*
double FF_Helm(double momentum, double atom_A)
{
    double s = 0.9 / .197326; 									// [GeV^-1]
    double r_n = 1.14*pow(atom_A,1.0/3.0) / .197326; 			// [GeV^-1]
    double x = momentum * r_n; 									// dimensionless if momentum in GeV

    return fabs((3.0 / x) * (sin(x)/x/x - cos(x)/x) * exp(-pow(momentum*s,2)/2.0));
}


double choose_element(char element[])
{
    double at_num = 0.0;

    if (strcmp(element,"Fe") == 0)
	{
        at_num = 55.85*0.93149; 				//Isotopic mass, average over isotopes in Earth
    }
    else if (strcmp(element,"Pb") == 0)
	{
        at_num = 207.2*0.93149; 				//Isotopic mass, average over isotopes in Earth
    }
    else if (strcmp(element,"Si") == 0)
	{
        at_num = 28.08*0.93149; 				//Isotopic mass, average over isotopes in Earth
    }
    else{
        printf("I don't know that element!\n");
    }

    return at_num;
}

long double threelayer_density(char element[2], double det_pos, double dist, double cos_theta)
{
    long double element_density,el_dens_crust,el_dens_mantle,el_dens_core;

    if (strcmp(element,"Fe") == 0){
        el_dens_crust 	= FE_crust;
        el_dens_mantle 	= FE_mantle;
        el_dens_core	= FE_core;
    }
    else if (strcmp(element,"Pb") == 0){
        el_dens_crust 	= Pb_crust;
        el_dens_mantle 	= Pb_mantle;
        el_dens_core	= Pb_core;
    }
    else if (strcmp(element,"Si") == 0){
        el_dens_crust     = Si_crust;
        el_dens_mantle     = Si_mantle;
        el_dens_core    = Si_core;
    }
    else{
        el_dens_crust	= 0;
        el_dens_mantle 	= 0;
        el_dens_core	= 0;
    }

  // Figure out if it hits Earth in the crust, mantle, core... or at all!
  // Densities are modulo dens_mult

    long double cos_thetastar 	 = (pow(det_pos,2) + pow(dist,2) - pow(R_Earth,2))  / (2 * dist * det_pos);
    long double cos_thetaMantle = (pow(det_pos,2) + pow(dist,2) - pow(R_mantle,2)) / (2 * dist * det_pos);
    long double cos_thetaCore   = (pow(det_pos,2) + pow(dist,2) - pow(R_core,2))   / (2 * dist * det_pos);

    if (cos_theta < cos_thetastar) 			// Hits outside Earth sphere
    	return 0;  
    else if (dist < det_pos - R_mantle || dist > det_pos + R_mantle){
        element_density = el_dens_crust;
    }
    else if (dist < det_pos - R_core || dist > det_pos + R_core){
        if (cos_theta < cos_thetaMantle)
            element_density = el_dens_crust;
        else
            element_density = el_dens_mantle;
    }
    else{
        if (cos_theta < cos_thetaMantle)
            element_density = el_dens_crust;
        else if (cos_theta < cos_thetaCore)
            element_density = el_dens_mantle;
        else
            element_density = el_dens_core;
    }
    return element_density;
}
*/

//=========================================================================
//---------------------- Lifetime (new) -----------------------------------
//=========================================================================

long double lifetime_MiDM(long double mDM, long double delta, double gM)
{
    long double lifetime;
    lifetime = pow(5.65/delta,3) * pow(mDM/1000.0,2) / gM / gM;  //[sec]
    return lifetime;
}

//========================================================================================================================
//----------------------------------------------- Spin-dependent Form Factors --------------------------------------------
//========================================================================================================================

double FF_Wick(double momIN, char element[], int atom_A)
{
    char filename[40];
    char line[100]; 					// temp buffer
    float momentum = 0, FF_res = 0;
    float mom_last = 0;
    
    //sprintf(filename, "isodata.txt");
    sprintf(filename,"%0d%s/WM00.txt", atom_A, element);
    //printf("file = %s\n",filename);
    
    FILE *fptr = fopen(filename, "r"); // take file pointer
    
    if(fptr == NULL) {
        perror("Error opening file");
        return(-1);
    }
    
    //if( feof(fptr) ) return -2 ;

    //printf("before: %f %f\n",momentum,FF_res);

    while((fgets(line, 512, fptr)) != NULL){
        //printf("hey\n");
        sscanf(line, "%f %f", &momentum, &FF_res);
        if(momIN < momentum && momIN > mom_last) {
            //printf("%f < %f < %f\n",mom_last,momIN,momentum);
            //printf("Data from the file:\n%s", line);
            //printf("Closing the file %s\n",filename);
            fclose(fptr);
            return FF_res;
        }
        //printf("%f %f %f\n",mom_last,momIN,momentum);

        //printf("after: %f %f\n",momentum,FF_res);
        
        mom_last = momentum;
        
    }
    //printf("%f\n",FF_res);
    //printf("Closing the file %s\n",filename) ;
    //free(fptr);
    fclose(fptr);
    return FF_res;
}

//============================================= Read details of isotope ==================================================

double * read_isodata(char eltIN[], int Ain)
{
    static double eltdata[8];
    int eltIN_len = strlen(eltIN);

    //printf("%d\n",strlen(element));
    int Znum, Anum, Znum_last = 0;
    float spin, moment, abundance;
    long double thalf, Dcore, Dmantle, Dcrust;
    char elt[2];
    //int test;
    
    char line[1000]; // temp buffer
    
    FILE *fptr = fopen("isodata.txt", "r"); // take file pointer
    
    int find_result = 0;
    
    if( fptr == NULL ) 
	{
        perror( "isodata.txt" );
        return eltdata;
        //exit( EXIT_FAILURE );
    }
    
    // read the file and find the match
    
    while((fgets(line, 512, fptr) != NULL))
	{     // && (find_result == 0)) {
        sscanf(line, "%d %s %d %Lf %f %f %f %Lf %Lf %Lf", &Znum, elt, &Anum, &thalf, &spin, &moment, &abundance, &Dcore, &Dmantle, &Dcrust);
        
        //printf("%s \n",line);
        if(Znum > Znum_last && find_result > 0)
		{
            fclose(fptr);
            //printf("%f\n",eltdata[6]);
            return eltdata;
        }
        
        if(elt[0] == eltIN[0] && Ain == Anum)
		{
            if(eltIN_len==1 || elt[1] == eltIN[1])
			{

            //printf("%d %s %d %f %f %f %f %lf %lf %lf\n", Znum, eltIN, Anum, thalf, spin, moment, abundance, Dcore, Dmantle, Dcrust);
            
                eltdata[0] = Znum;
                if(find_result==0)
				{
                    eltdata[1] = Anum;
                    eltdata[2] = spin;
                    eltdata[3] = moment;
                    eltdata[4] = abundance;
                    eltdata[5] = Dcore/dens_mult;
                    eltdata[6] = Dmantle/dens_mult;
                    eltdata[7] = Dcrust/dens_mult;
                    //return eltdata;
                    //printf("%f\n",eltdata[7]);
                }
				/*
                else if(find_result==1)
				{
					eltdata[8] = Znum;
                    eltdata[9]  = Anum;
                    eltdata[10]  = spin;
                    eltdata[11] = moment;
                    eltdata[12] = abundance;
                    eltdata[13] = Dcore/dens_mult;
                    eltdata[14] = Dmantle/dens_mult;
                    eltdata[15] = Dcrust/dens_mult;
                }
                */
                Znum_last = Znum;
                find_result++;
            }
        }
        
        if(Znum == 82 && find_result == 0){
            printf("Error! Not a valid element. \n");
            //free(fptr);
            fclose(fptr);
            //printf("%f\n",eltdata[7]);
            return eltdata;
        }
    }
    
    printf("Something's wrong.\n");
    //free(fptr);
    fclose(fptr);
    return eltdata;
}

double whichlayer(double det_pos, double dist, double cos_theta)
{
    // Figure out if it hits Earth in the crust, mantle, core... or at all!
    // Miss --> 0; Crust --> 1; Mantle --> 2; Core --> 3
    
    double cos_thetastar   = (pow(det_pos,2) + pow(dist,2) - pow(R_Earth,2))  / (2 * dist * det_pos);
    double cos_thetaMantle = (pow(det_pos,2) + pow(dist,2) - pow(R_mantle,2)) / (2 * dist * det_pos);
    double cos_thetaCore   = (pow(det_pos,2) + pow(dist,2) - pow(R_core,2))   / (2 * dist * det_pos);
    
    if (cos_theta < cos_thetastar) return 0;                                    // Hits outside Earth sphere
    else if (dist < det_pos - R_mantle || dist > det_pos + R_mantle) return 1;  // Hits Crust
    else if (dist < det_pos - R_core || dist > det_pos + R_core)
	{
        if (cos_theta < cos_thetaMantle) return 1;                              // Hits Crust
        else return 2;                                                          // Hits Mantle
    }
    else
	{
        if (cos_theta < cos_thetaMantle) return 1;                              // Hits Crust
        else if (cos_theta < cos_thetaCore) return 2;                           // Hits Mantle
        else return 3;                                                          // Hits Core
    }
}

//========================================================================================================================
//----------------------------------------------- Response Functions and MeltSqu -----------------------------------------
//========================================================================================================================

//---------------------------------- Reads output file from Bigstick + Mathematica ---------------------------------------

double * getWcoefficients(int Zin, int Ain, 
							int Oselector, int iso1, int iso2) 
{
	static double W[10];
	
    char line[1000]; // temp buffer
	
	int Znum, Anum, Op, tau, taup, find_result=0;
    
    FILE *fptr = fopen("NucResFuncCoefficientsISO.txt", "r"); // take file pointer
    
    //int find_result = 0;
    
    if( fptr == NULL ) 
	{
        perror( "NucResFuncCoefficients.txt" );
        return 0;
        //exit( EXIT_FAILURE );
	}
	
    while((fgets(line, 512, fptr) != NULL) && (find_result == 0)) 
	{
        sscanf(line, "%d,%d,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", 
					&Znum, &Anum, &Op, &tau, &taup, &W[0], &W[1], &W[2], &W[3], &W[4], &W[5], &W[6], &W[7], &W[8], &W[9]);
		
		if ( Znum == Zin && Anum == Ain && Op == Oselector && tau == iso1 && taup == iso2) 
		{
			fclose(fptr);
			return W;
		}
	}
    
	printf("We did not compute response functions for that element (yet)!");
	fclose(fptr);
	return 0;
}

//-------------------------- Wilson Coefficients (right now, not in use) -----------------------

long double cWilson(int counter, int isospin, long double mDM, long double momentum)
{
	long double cp, cn;
	//long double Dp = (eq/2/mn)*dp;
	//long double Dn = (eq/2/mn)*dn;
	//printf("q = %Lf\n",momentum);
	
	if(counter == 1)
	{
		cp = pow(eq/2/mDM/cSqu,2);
		cn = 0;
	}
	else if(counter == 4)
	{
		cp = dp*eq*eq/mn/mDM/cSqu;
		cn = dn*eq*eq/mn/mDM/cSqu;
	}
	else if(counter == 5)
	{
		cp = -mn/mDM/cSqu*pow(eq/momentum,2);
		cn = 0;
	}
	else if (counter == 6)
	{
		cp = -dp*mn/mDM/cSqu*pow(eq/momentum,2);
		cn = -dn*mn/mDM/cSqu*pow(eq/momentum,2);
	}
	else 
	{
		//printf("This operator has zero DM response for MiDM!");
		cp = 0;
		cn = 0;		
	}
	
	if(isospin == 0) return (cp + cn)/2.0;
	else if(isospin == 1) return (cp - cn)/2.0;
	else {
		printf("Invalid Isospin Value!");
		return 0;
	}
}

long double cWsqu(int iIN, int iso1, int iso2, long double mDM, long double momentum)
{
	return cWilson(iIN,iso1,mDM,momentum)*cWilson(iIN,iso2,mDM,momentum);
}


//-------------------------- Computes the polynomial function that determines the nuclear response -----------------------

long double NuclearResponse(int Zin, int Ain, 
							double Win[6][2][2][10], int Oselector, int iso1, int iso2, 
							double y)
{
	double W[10];
	int nmax_y=9;
	
	for(int i=0; i<=nmax_y; i++) W[i] = Win[Oselector][iso1][iso2][i];
	
	long double sum = 0;
	for(int i=0; i<=nmax_y; i++) sum += W[i]*pow(y,i);
	
	//return 1; // TESTING
	return exp(-2*y)*sum;
		//(W[0]+W[1]*y+W[2]*y^2+W[3]*y^3+W[4]*y^4+W[5]*y^5+W[6]*y^6)
}

//============================================= DM Response functions ====================================================
//---------- Oselector chooses from 5 available operators: M, Sigma'', Sigma', Delta, DeltaSigma' = 1, 2, 3, 4, 5 --------
//------------------------- dimensionful inputs mDM, delta, mN, ER are in units of GeV/c^2 -------------------------------

//double DMResponse(double mDM, double delta, double mN, double vin, double ER, int Oselector, int iso1, int iso2)
long double DMResponse(long double mDM, long double delta, long double mN, 
						long double ER, double vin, 
						int Oselector, int iso1, int iso2)
{
	double q 		= sqrt(2*(mN*cSqu)*ER);
	double muT 		= mDM*mN/(mDM+mN)*cSqu;
	double vT 		= vin/cLight;//sqrt(pow(vin/cLight,2) + pow(q/2/muT,2));
	double vminT 	= q/(2*muT) + delta*cSqu/q;
	long double temp;
	
	// TESTING
	// printf("mDM, muT, ER, q   = %Lf, %f, %Lf, %f\n",mDM*cSqu,muT,ER,q);
	// printf("term1, term2 = %f, %Lf\n",q/2/muT,delta*cSqu/q);
	// printf("vf [km/sec], q [GeV]: vT [c] , vTmin [c], diffvT [c]     = %4.3f, %4.5f: %4.7f, %4.7f, %4.5fE-6\n\n",vin, q, vT,vminT,(pow(vT,2.0) - pow(vminT,2.0))*pow(10,6));
	
	if (Oselector == 1)									// Operator M
	{ 
		temp = jchi*(jchi+1)/3 * cWilson(5, iso1, mDM, q)*cWilson(5, iso2, mDM, q);
		temp *= pow(q/mn,2) * (pow(vT,2.0) - pow(vminT,2.0));
		temp += cWilson(1, iso1, mDM, q)*cWilson(1, iso2, mDM, q);
		return temp;
	}
	else if (Oselector == 2)							// Operator Sigma'' (Note: response is 0 for MiDM)
	{ 
		return 0;
	}
	else if (Oselector == 3)							// Operator Sigma'
	{ 
		temp = jchi*(jchi+1)/12 * cWilson(4, iso1, mDM, q)*cWilson(4, iso2, mDM, q);
		return temp;
	}
	else if (Oselector == 4)							// Operator Delta
	{ 
		temp = jchi*(jchi+1)/3 * cWilson(5, iso1, mDM, q)*cWilson(5, iso2, mDM, q) * pow(q/mn,2);
		return temp;
	}
	else if (Oselector == 5)							// Operator DeltaSigma' (Note: asymmetric in \tau, \tau')
	{ 	
		temp = jchi*(jchi+1)/3 * cWilson(5, iso1, mDM, q)*cWilson(4, iso2, mDM, q);		
		return temp;
	}
	else
	{
		printf("Invalid operator!\n");
		return 0;
	}
}

long double ResponseProduct(double Zin, double Ain, 
								long double mDM, long double delta, long double mN, 
								long double ER, double vin, 
								double Win[6][2][2][10], int Oselector)
{
	long double temp = 0;
	
	long double q = sqrt(2*(mN*cSqu)*ER);
	long double b = 5.0677*sqrt(41.467/(-25*pow(Ain,-2.0/3.0) + 45*pow(Ain,-1.0/3.0)));
	long double y = pow(b*q/2,2);
		
	for (int tau=0; tau<=1; tau++)
	{
		for (int taup=0; taup<=1; taup++)
		{
			temp += DMResponse(mDM, delta, mN, ER, vin, Oselector,tau,taup)
								 * NuclearResponse(Zin, Ain, Win, Oselector,tau,taup, y);
			
			// TESTING
			// printf("Nuc res = %Lf \n", NuclearResponse(Zin, Ain, Win, Oselector,tau,taup, y));
			// printf("DM res = %Lf \n", DMResponse(mDM, delta, mN, ER, vin, Oselector,tau,taup));
		}
	}
	
	//printf("b, y1 = %Lf, %Lf\n\n",b,y);
	
	return temp;
}

//============================================== Matrix Element Squared ==================================================


//--------------------------------- From product of DM * Nuc response functions ------------------------------------------

long double MatEltSqu(double Zin, double Ain, 
						long double mDM, long double delta, long double mN, 
						long double ER, double vin,
						double Win[6][2][2][10])
{
	long double q = sqrt(2*(mN*cSqu)*ER);
	long double temp=0;
	
	temp += ResponseProduct(Zin,Ain,mDM,delta,mN,ER,vin,Win,1); 					// Operator M
	//printf("product1 = %Lf\n",temp);
	temp += ResponseProduct(Zin,Ain,mDM,delta,mN,ER,vin,Win,2); 					// Operator Sigma''
	//printf("product2 = %Lf\n",temp);
	temp += ResponseProduct(Zin,Ain,mDM,delta,mN,ER,vin,Win,3); 					// Operator Sigma'
	//printf("product3 = %Lf\n",temp);
	temp += pow(q/mn,2)*ResponseProduct(Zin,Ain,mDM,delta,mN,ER,vin,Win,4); 		// Operator Delta
	//printf("product4 = %Lf\n",temp);
	temp += pow(q/mn,2)*ResponseProduct(Zin,Ain,mDM,delta,mN,ER,vin,Win,5); 		// Operator DeltaSigma'
	//printf("product5 = %Lf\n",temp);
	//printf("DM=%Lf, Nuc=%Lf\n",DMResponse(mDM, delta, mN, ER, vin, 1,0,0),NuclearResponse(Zin, Ain, 1,0,0, y));
	return temp;
}

// ----------------------------------------- Using "copy-paste" method ---------------------------------------------------

double matrixelementSQU(char element[], int Anum, 
						long double mDM, long double delta, long double mN, 
						double speedtot, long double ER)
{
    // Runs only for 27Al
	// Need to check output
    
    double E = exp(1);
 
    if (strcmp(element,"Al") == 0){ // Checked: Matches perfectly FormFactors_JE.nb for same inputs
        return pow(E,-87.13130844161698*ER*mN)*pow(ER,-2)*pow(mDM,-2)*pow(mN,-2)*
                (delta*ER*mDM*(mDM + mN)*(-1.9376975290370302 + ER*mN*(138.5305596071742 +
                                                           ER*mN*(-3373.7871900307237 + ER*mN*(31889.09782540305 - 105728.80991532124*ER*mN))))
                 + (-0.9688487645185151 + ER*mN*(69.2652798035871 + ER*mN*(-1686.8935950153618
                            + ER*mN*(15944.548912701524 - 52864.40495766062*ER*mN))))
                 * pow(delta,2)*pow(mDM,2) + ER*mN*(ER*mN*(0.9688487645185152 +
                                             ER*mN*(-69.26527980358712 + ER*mN*(1686.893595015362 + ER*mN*(-15944.548912701528 + 52864.40495766063*ER*mN))) +
                                             (0.07501681445451591 + ER*mN*(-7.30799565567299 + ER*mN*(278.9355010962528 + ER*mN*(-4217.987459077475 + 41758.05338325841*ER*mN))))*
                                             pow(mDM,2)) + (2.1559792642980095e-11 + ER*mN*(-1.5413603491206154e-9 +
                                                                                            ER*mN*(3.753844506100664e-8 + ER*mN*(-3.5481406482934024e-7 + 1.1763916627876455e-6*ER*mN))))*pow(mDM,2)*pow(speedtot,2)));
    }
    else
        return 0;
}

//=========================================================================
//--------------- Compute The Integrand -----------------------------------
//=========================================================================

double sixDim_Integrand (double *k, size_t dim, void *my_params)
{
    struct param_type params = *(struct param_type *)my_params;
	
	double gM = params.gM; // Can update this coupling value later

  	//---------------------------- DM sector constants -------------------------------------------
    long double mDM 		= params.mDM / cSqu; 									// [GeV/c^2]
    long double delta 		= params.delta * pow(10.0,-6.0) / cSqu; 						// [GeV/c^2]
    long double lifetime    = lifetime_MiDM(params.mDM,params.delta,gM) * params.lifemult;   				// [sec]
    long double M_Exc 	    = (mDM + delta);	   	 		 		     // DM Excited State Mass 	// [GeV/c^2]
	 
	//printf("lifetime = %Lf \n",lifetime);
    
	//----------------------------- Target constants-------------------------------------------
    int Anum = params.Anum;
	int Znum = params.Znum;
    long double atom_A     = Anum*0.93149;
    long double mN         = atom_A / cSqu;            					// Nuclear Mass         	  // [GeV/c^2]
    long double m_Red      = mDM * mN / (mDM + mN);    					// DM-Nucleus Reduced Mass	 // [GeV/c^2]
	//printf("mRed = %Lf\n",m_Red*cSqu);

  	//----------------------------- Choose Detector ------------------------------------------
    double det_vars[3];
    detector_choice(params.detector,det_vars);

    double det_radius   = det_vars[0];							// Radius of detector		// [meters]
    double det_pos	    = det_vars[1];						// Depth of detector		// [meters]
    double det_theta	= degrees_to_radians(det_vars[2]);				// Angular position of detector (latitude -> polar angle)

  	//---------------------------- Integration Variables -------------------------------------
  	//       Integration over scatter sites: "Detector Coordinates" (detector at origin, z points towards core)
    double dist_DET	     = k[0];	   					// Scatter site distance to detector // [km]
    double cos_theta_DET = k[1];	   					// Scatter site polar angle to detector
    double phi_r_DET     = k[2];	  					// Scatter site azimuthal angle to detector

    double r_DET_Sph[3] = {dist_DET, cos_theta_DET, phi_r_DET};

  	//        Integration over DM velocities: Galactic Coordinates
    double speed_in_GAL  = k[3];    					// Incoming DM speed		[km/sec]
    double cos_gamma_GAL = k[4];	   					// Cosine of incoming DM polar angle
    double phi_v_GAL     = k[5];	   					// Incoming DM aximuthal angle
	
    double v_in_GAL_Sph[3] = {speed_in_GAL, cos_gamma_GAL, phi_v_GAL};
	
  	// 				Boltzmann distribution
    double fMB 		= exp(-pow(speed_in_GAL / v0 , 2)) / f_norm;
	
	// TESTING
  	//dist_DET		= 5369.600652;  // TEST
  	//cos_theta_DET 	= cos(-0.000003);   // TEST
  	//phi_r_DET		= -1.311157;	// TEST
  	//speed_in_GAL 	= 300;  // TEST
  	//cos_gamma_GAL 	= .1;   // TEST
  	//phi_v_GAL 		= M_PI; // TEST
	    
  	// Rotate v_in_GAL to cartesian coordinates
    double *p 		= sph_to_cart(v_in_GAL_Sph);
    double v_in_GAL_Cart[3];
    for (int i = 0; i < 3; i++) v_in_GAL_Cart[i] = *(p+i);

  	// Calculate Earth's velocity u_E at given time
    p    = compute_Earth_velocity_TOT(params.year, params.month, params.day, params.hour);
    double u_E[3];
    for (int i = 0; i < 3; i++ ) u_E[i] = *(p+i);
    
    //printf("v_E [km/sec] = {%f,%f,%f}\n",u_E[0],u_E[1],u_E[2]);
  
  	// Add Earth's velocity u_E to v_in_GAL to find v_tot_GAL
    p = vector_operation(v_in_GAL_Cart,u_E,'s');
    double v_tot_GAL_Cart[3];
    for (int i = 0; i < 3; i++ ) v_tot_GAL_Cart[i] = *(p+i);
    
    double speed_tot	= vector_magnitude(v_tot_GAL_Cart);
  
  // ===================================================================
  // ||||||||||||||||| BUNCH OF ROTATIONS ||||||||||||||||||||||||||||||
  // ===================================================================

  	// Earth's rotation period is sidereal day
    double nTropical = date_to_J2000(params.year, params.month, params.day, params.hour); // number of days since J2000
    double nSidereal = nTropical * 24.0 / 23.9345;								 // number of sidereal days since J2000

    // Convert v_tot from Galactic to Equatorial coordinates (Cartesian)
    p             = coordinates_gal_to_equat(v_tot_GAL_Cart);
    double v_tot_EQU_Cart[3];
    for (int i = 0; i < 3; i++ ) v_tot_EQU_Cart[i] = *(p+i);

    // Convert v_tot from Equatorial to Detector coordinates (Cartesian) (it's a velocity, doesn't care about origin!)
    p            = rotate_cartesian_z(v_tot_EQU_Cart, -nSidereal * 2.0 * M_PI);
	//p            = rotate_cartesian_z(v_tot_EQU_Cart, -2.0*M_PI*params.hour/24.0); // for tropical days (test)
    double v_tot_DET_Cart[3];
    for (int i = 0; i < 3; i++ ) v_tot_DET_Cart[i] = *(p+i);

    p            = rotate_cartesian_y(v_tot_DET_Cart, M_PI+det_theta);
    for (int i = 0; i < 3; i++ ) v_tot_DET_Cart[i] = *(p+i);

    // Convert r_DET from Spherical to Cartesian basis
    p            = sph_to_cart(r_DET_Sph);
    double r_DET_Cart[3];
    for (int i = 0; i < 3; i++ ) r_DET_Cart[i] = *(p+i);

    double cos_alpha 	= - cos_angle_between(v_tot_DET_Cart,r_DET_Cart);
	
	//double r_minus[3] 	= {-r_DET_Cart[0],-r_DET_Cart[1],-r_DET_Cart[2]};
	//double cos_alpha 	= cos_angle_between(v_tot_DET_Cart,r_minus);
	
    //double alpha = acos(cos_alpha);
    //double sin_alpha = sin(alpha);
    
	//------------------- Elemental density, normalized by dens_mult (re-multiplied in mult factor below) -------------

    long double elt_density;	

    int layer = whichlayer(det_pos, dist_DET, cos_theta_DET);
    if(layer == 0)      elt_density = 0;
    else if(layer == 1) elt_density = params.Dcrust;
    else if(layer == 2) elt_density = params.Dmantle;
    else if(layer == 3) elt_density = params.Dcore;
	//if(layer != 0) printf("%d %Lf\n",layer,elt_density);
    
    // ------------------------------------ Lower bound on velocity! -----------------------------------
    
    double vmin         	= cLight * sqrt(2*delta / (mDM * pow(1 - m_Red / mN,2) + pow(m_Red,2) / mN));
    double vLower       	= vmin - vector_magnitude(u_E);
    if(vLower > vEsc) return 0;
	//printf("vmin = %f\n",vLower);
	//printf("v = %f\n",speed_tot);

    //---------------------------------------- Other Kinematics  ---------------------------------------
	
    if(m_Red * pow(speed_tot/cLight,2) <= 2 * delta) return 0;					// enforces real vf_CM
	
	//printf("vmin = %f\n",vLower);
	//printf("v = %f\n",speed_tot);
	
	
	//printf("speed_tot, m_N/mDM, 2mN delta/mDM^2/vtot^2 = %f,%Lf,%Lf\n",speed_tot,mN/mDM,2*delta*mN*cSqu/mDM/mDM/speed_tot/speed_tot);
    	
    double cos_alphastar    = 1 - mN/mDM + 2*mN*delta*pow(cLight / mDM / speed_tot,2);
    cos_alphastar        	= sqrt((1 + mN/M_Exc)*cos_alphastar); 				// Maximum Outgoing Scattering Angle
    double alpha_star       = acos(cos_alphastar);
	double alpha			= acos(cos_alpha);
	
	//printf("%f\n",cos_alphastar);
	
	// -------------__ TEST TEST TEST __---------------------------
    //dist_DET = 100.;
	//speed_tot = 730;
	//cos_alpha = 0.9998;
	//cos_alphastar = -1.0;//-0.99999999;
	
	
    // Require outgoing cos_alpha is kinematically allowed 
    if (cos_alpha < cos_alphastar) return 0.0;
        	
	// Difference between (cos of) actual scattering angle and maximal one    
    long double temp     	= sqrt(fabs(pow(cos_alpha,2) - pow(cos_alphastar,2)));

    //Two solutions for the outgoing DM speed, fast and slow one, + and - [km/sec]
    double vf[2]     		= {mDM * speed_tot / (M_Exc + mN) * (cos_alpha + temp),
        					   mDM * speed_tot / (M_Exc + mN) * (cos_alpha - temp)};				
		

	long double jacobian[2];
	
	double vf_CM 	= cLight * sqrt( mN / M_Exc / (M_Exc + mN)) * sqrt(m_Red * pow(speed_tot/cLight,2) - 2*delta);

    temp 				= (mDM * speed_tot) / ((mDM + mN) * vf_CM); // this is u/w
	long double temp2	= 1.0 - temp*temp + pow(temp*cos(alpha),2);

	if (temp2 > 0.0) // equiv to w > sinalpha*u ?
	{
	  temp2		= ( 1.0 - temp*temp + 2 * pow(temp*cos(alpha),2) ) / sqrt(temp2);
    	  jacobian[0] = fabsl(2.0 * temp * cos(alpha) + temp2) * fabs(sin(alpha));
		  jacobian[1] = fabsl(2.0 * temp * cos(alpha) - temp2) * fabs(sin(alpha));
		  //jacobianPlus 	= fabs(2.0 * temp * cos(alpha) + temp2) * fabs(sin(alpha));
    	  //jacobianMinus	= fabs(2.0 * temp * cos(alpha) - temp2) * fabs(sin(alpha));
	}			
        		
    long double E_recoil[2] = 
			{(cSqu/2.0) * mDM * pow(speed_tot/cLight,2) - (cSqu/2.0) * M_Exc * pow(vf[0]/cLight,2) - delta*cSqu,
             (cSqu/2.0) * mDM * pow(speed_tot/cLight,2) - (cSqu/2.0) * M_Exc * pow(vf[1]/cLight,2) - delta*cSqu};

		// Test "copy/paste" method against new method of reading input files for response functions
			
	// ---------------- Differential cross section dsigma/dE = C * P_tot ----------------------------
			 
	long double probE[2];
		   		        
        // This is the equivalent of P_{tot} in Catena+Schwabe, for each of the two final recoil energies
    probE[0] = 2 * exp(-dist_DET / (lifetime * vf[0]));
	probE[0] *= sinh(det_radius / (2*lifetime * vf[0])) * jacobian[0];
	//printf("probE0 second = %Lf\n",probE[0]);
    probE[0] *= 4*M_PI/(2*params.jN+1);
	//printf("probE0 third = %Lf\n",probE[0]);
	probE[0] *= MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[0], speed_tot, params.Wcoeffs);
	//printf("probE0 final = %Lf\n",probE[0]);

    probE[1] = 2 * exp(-dist_DET / (lifetime * vf[1]));
	//printf("probE1 = %Lf\n",probE[1]);
	probE[1] *= sinh(det_radius / (2*lifetime * vf[1])) * jacobian[1];
	//printf("probE1 second = %Lf\n",probE[1]);
    probE[1] *= 4*M_PI/(2*params.jN+1);
	probE[1] *= MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[1], speed_tot, params.Wcoeffs);
		
	//printf("m_Red = %f\n",m_Red*cSqu);

//------------------ Final Expressions ----------------------------------
    
    long double suppression 			= pow(det_radius,2) / pow(alpha_star * dist_DET ,2);

	// Units depend on choice of cWilson normalization; think I have this right
	// hbar*c factor below implies that the final rate will be in 1/sec as desired
	long double dsigdcosthCM_prefactor	= pow(gM*m_Red*cSqu,2)/(2*M_PI);
	dsigdcosthCM_prefactor				*= sqrt(1-2*delta*cSqu/(m_Red*speed_tot*speed_tot));

    long double mult				= dsigdcosthCM_prefactor * suppression * dens_mult * DM_massdensity / mDM * hbarc * hbarc;
		//sigma_prefactor * suppression * dens_mult * DM_massdensity / mDM; //sigma_0	
	
	long double q1 = sqrt(2*(mN*cSqu)*E_recoil[0]);
	long double q2 = sqrt(2*(mN*cSqu)*E_recoil[1]);
	long double b = 5.0677*sqrt(41.467/(-25*pow(Anum,-2.0/3.0) + 45*pow(Anum,-1.0/3.0)));
	long double y1 = pow(b*q1/2,2);
	long double y2 = pow(b*q2/2,2);
	
	printf("(mDM,delta,gM(true))    (vinlab, voutlab1, voutlab2, voutCM)   (costhetalab)   10^10*(dsigdcosthetaCM(voutlab1),dsigdcosthetaCM(voutlab2))  (y1, y2) = (%4.0Lf GeV, %3.0Lf keV, %3.3f)     (%f, %f, %f, %f)    (%f)    (%Lf, %Lf)  (%Lf, %Lf)\n",
					mDM*cSqu, delta*cSqu*pow(10,6),sqrt(2)*gM, speed_tot, vf[0], vf[1], vf_CM, cos_alpha, 
					2*M_PI/(2*params.jN+1)*dsigdcosthCM_prefactor*MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[0], speed_tot, params.Wcoeffs)*pow(10,10), 
					2*M_PI/(2*params.jN+1)*dsigdcosthCM_prefactor*MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[1], speed_tot, params.Wcoeffs)*pow(10,10), y1, y2);
	
	//printf("sigpre = %Lf\n",dsigdcosthCM_prefactor);
	
    //double integrand = mult * pow(dist_DET*speed_in_GAL,2) * fMB * suppression * (dsigma_dcosCM * speed_tot) * elt_density * (prob[0] + prob[1]);
	long double integrand 			= mult * pow(dist_DET*speed_in_GAL,2) * fMB * speed_tot * 
												elt_density * (probE[0] + probE[1]); //sigma_prefactor *  suppression 
	
	/* TESTING
	printf("(m,delta,tau) = %4.0Lf GeV, %3.0Lf keV, %3.3Lf sec\n",mDM*cSqu, delta*cSqu*pow(10,6),lifetime);
	printf("(year,month,hour) = %4.0f, %2.0f, %2.0f\n",params.year, params.month, params.hour);
	printf("(rho,v_esc,v_0, ????) = %3.3f GeV/cm^3, %3.0f km/sec, %3.0f km/sec, ????\n",DM_massdensity*pow(10,-15)*cSqu,vEsc,v0);
	printf("atomic number of target: %d\n",Anum);
	printf("nuclear spin of target:  %f\n",params.jN);
	printf("elemental densities (core,mantle,crust): %3.3LfE32,%3.3LfE32,%3.3LfE32\n\n", params.Dcore, params.Dmantle, params.Dcrust);
	*/	
	/*
	printf("integrand = %.10Lf\n\n", integrand);
	
	printf("ResponseProduct [10^{-10} GeV^{-4}] = %Lf, %Lf, %Lf, %Lf, %Lf \n", 
	ResponseProduct(params.Znum,params.Anum,params.mDM/cSqu, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, E_recoil[0], speed_tot,params.Wcoeffs,1)*pow(10,10),
	ResponseProduct(params.Znum,params.Anum,params.mDM/cSqu, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, E_recoil[0], speed_tot,params.Wcoeffs,2)*pow(10,10),
	ResponseProduct(params.Znum,params.Anum,params.mDM/cSqu, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, E_recoil[0], speed_tot,params.Wcoeffs,3)*pow(10,10),
	ResponseProduct(params.Znum,params.Anum,params.mDM/cSqu, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, E_recoil[0], speed_tot,params.Wcoeffs,4)*pow(10,10),
	ResponseProduct(params.Znum,params.Anum,params.mDM/cSqu, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, E_recoil[0], speed_tot,params.Wcoeffs,5)*pow(10,10));
	*/
							
    return integrand;
	

	// if TESTING THE OUTPUTS THOROUGHLY, use the code below:
	/*
	double q[2]		= {sqrt(2*(mN*cSqu)*E_recoil[0]), sqrt(2*(mN*cSqu)*E_recoil[1])};
	//double muT 		= mDM*mN/(mDM+mN)*cSqu;
	//double vT[2] 	= {sqrt(pow(speed_tot/cLight,2) + pow(q[0]/2/muT,2)), 
	//				sqrt(pow(speed_tot/cLight,2) + pow(q[1]/2/muT,2))};
	//double vminT[2] = {q[0]/(2*muT) + delta*cSqu/q[0],
	//						q[1]/(2*muT) + delta*cSqu/q[1]};	
							
	long double b = 5.0677*sqrt(41.467/(-25*pow(params.Anum,-2.0/3.0) + 45*pow(params.Anum,-1.0/3.0)));
	long double y[2] = {pow(b*q[0]/2,2),pow(b*q[1]/2,2)};

	temp  = MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[0], speed_tot, params.Wcoeffs);
	temp *= 4*M_PI/(2*params.jN+1);
	
	temp2 = MatEltSqu(Znum, Anum, mDM, delta, mN, E_recoil[1], speed_tot, params.Wcoeffs);
	temp2 *= 4*M_PI/(2*params.jN+1);							
	
	printf("~~~~~~~~~~~~~~~~~~ BENCHMARK POINT ~~~~~~~~~~~~~~~~~~\n\n");
	printf("speed earth frame |  costhetamax | cos_to_scattersite | voutCM \n %f        |  %f    | %f          | %f \n\n",
	speed_tot, cos_alphastar, cos_alpha, vf_CM);
	
	printf("  Z,  A,  mDM [GeV], delta [keV]\n  %d,  %d,  %.0Lf,        %.0Lf \n\n ",
					Znum, Anum, mDM*cSqu, delta*cSqu*pow(10,6));
										
	printf("------------ Solution 1 ------------ \n\n             q = %f\n" ,sqrt(2*(mN*cSqu)*E_recoil[0]));
	printf("vin [km/sec], vf [km/sec] = %4.3f, %4.3f\n", speed_tot, vf[0]);
	printf("c[0] = %LfE-9, %LfE-5, %LfE-2, %LfE-3\n",
			cWilson(1, 0, mDM, q[0])*pow(10,9),
			cWilson(4, 0, mDM, q[0])*pow(10,5),
			cWilson(5, 0, mDM, q[0])*pow(10,2),
			cWilson(6, 0, mDM, q[0])*pow(10,3));
	printf("c[1] = %LfE-9, %LfE-5, %LfE-2, %LfE-3\n",
			cWilson(1, 1, mDM, q[0])*pow(10,9),
			cWilson(4, 1, mDM, q[0])*pow(10,5),
			cWilson(5, 1, mDM, q[0])*pow(10,2),
			cWilson(6, 1, mDM, q[0])*pow(10,3));
	for(int i = 0; i<2; i++)
		for(int j = 0; j<2; j++)
			printf("(t,tp) = (%d,%d) -- \n    RM, WM                     = %.4LfE-13 | %.4Lf, \n    Rsigmap, Wsigmap           = %.4LfE-10 | %.4Lf, \n    Rdelta, Wdelta             = %.4LfE-7  | %.4Lf, \n    Rdeltasigmap, Wdeltasigmap = %.4LfE-7  | %.4Lf \n",
				i,j,
				DMResponse(mDM, delta,mN, E_recoil[0], speed_tot, 1, i,j)*pow(10,13),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 1, i,j, y[0]),
				DMResponse(mDM, delta,mN, E_recoil[0], speed_tot, 3, i,j)*pow(10,10),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 3, i,j, y[0]),
				DMResponse(mDM, delta,mN, E_recoil[0], speed_tot, 4, i,j)*pow(10,7),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 4, i,j, y[0]),
				DMResponse(mDM, delta,mN, E_recoil[0], speed_tot, 5, i,j)*pow(10,7),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 5, i,j, y[0])
				);
				//printf("y = %Lf %Lf\n\n",y[0],y[1]);
								
	printf("MatEltSqu1 [10^{-10} GeV^-4] = %Lf\n\n", temp*pow(10,10));

	printf("------------ Solution 2 ------------ \n\n             q = %f\n" ,sqrt(2*(mN*cSqu)*E_recoil[1]));
	printf("vin [km/sec], vf [km/sec] = %4.3f, %4.3f\n", speed_tot, vf[1]);
	printf("c[0] = %LfE-9, %LfE-5, %LfE-2, %LfE-3\n",
			cWilson(1, 0, mDM, q[1])*pow(10,9),
			cWilson(4, 0, mDM, q[1])*pow(10,5),
			cWilson(5, 0, mDM, q[1])*pow(10,2),
			cWilson(6, 0, mDM, q[1])*pow(10,3));
	printf("c[1] = %LfE-9, %LfE-5, %LfE-2, %LfE-3\n",
			cWilson(1, 1, mDM, q[1])*pow(10,9),
			cWilson(4, 1, mDM, q[1])*pow(10,5),
			cWilson(5, 1, mDM, q[1])*pow(10,2),
			cWilson(6, 1, mDM, q[1])*pow(10,3));
	for(int i = 0; i<2; i++)
		for(int j = 0; j<2; j++)
			printf("(t,tp) = (%d,%d) -- \n    RM, WM                     = %.4LfE-13 | %.4Lf, \n    Rsigmap, Wsigmap           = %.4LfE-10 | %.4Lf, \n    Rdelta, Wdelta             = %.4LfE-7  | %.4Lf, \n    Rdeltasigmap, Wdeltasigmap = %.4LfE-7  | %.4Lf \n",
				i,j,
				DMResponse(mDM, delta,mN, E_recoil[1], speed_tot, 1, i,j)*pow(10,13),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 1, i,j, y[1]),
				DMResponse(mDM, delta,mN, E_recoil[1], speed_tot, 3, i,j)*pow(10,10),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 3, i,j, y[1]),
				DMResponse(mDM, delta,mN, E_recoil[1], speed_tot, 4, i,j)*pow(10,7),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 4, i,j, y[1]),
				DMResponse(mDM, delta,mN, E_recoil[1], speed_tot, 5, i,j)*pow(10,7),
				NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 5, i,j, y[1])
				);
			
	printf("MatEltSqu2 [10^{-10} GeV^-4] = %Lf\n\n", temp2*pow(10,10));
	
	//printf("dsig/dcosthLAB [10^{-10} GeV^-2] = %Lf\n\n",
		 			//dsigdcosthLAB_prefactor*temp*pow(10,10));
	
	printf("------------ RESULT ------------\n\n");
	printf("dsig/dcosthCM [10^{-50} km^2]  = %Lf\n\n",
		 			dsigdcosthLAB_prefactor*(temp+temp2)* hbarc * hbarc * pow(10,50));
	printf("dsig/dcosthLAB [10^{-50} km^2] = %Lf\n\n",
					dsigdcosthLAB_prefactor*(temp*jacobian[0]+temp2*jacobian[1])* hbarc * hbarc * pow(10,50));
		*/
}

//========================================================================
//------------------ Monte Carlo -----------------------------------------
//========================================================================

double *MC_Ints(void *my_params, double results[2])
{
    struct param_type params = *(struct param_type *)my_params;

    int runcount		= 0; // Keeps track of Vegas run #

//--------------THIS IS ALL A REPEAT OF EARLIER, come back later to make nicer-------------------
  
    long double mDM 		= params.mDM / cSqu;
    long double delta 		= params.delta*pow(10.0,-6) / cSqu; // [keV/c^2]

    
long double mN         = params.Anum*0.93149/cSqu; 							// Nuclear Mass [GeV/c^2]
						//isodata[1]*0.93149/cSqu;       //choose_element(params.element)/cSqu; 
    long double m_Red 		= mDM * mN / (mDM + mN);	 // DM-Nucleon Reduced Mass	[GeV/c^2]

    double vmin 			= cLight * sqrt(2*delta / (mDM * pow(1 - m_Red / mN,2) + pow(m_Red,2) / mN));

    double * p	= compute_Earth_velocity_TOT(params.year, params.month, params.day, params.hour);
    double u_E[3];
    for (int i = 0; i < 3; i++ ) u_E[i] = *(p+i);

    double vLower		= vmin - vector_magnitude(u_E);
    if(vLower > vEsc)
	{
        vLower=vEsc-.0000001;
        printf("------> Kinematically Not Allowed! <--------\n");
        runcount = 1000;
        
        results[0] = 0.0;
        results[1] = 0.0;
    }
    else
	{
  
    //double vLower = 0.;
    
    double res, err;
  	// vars are 	   ( L, 	    costheta,  phir, 	    v,       cosgamma,  phiv)
    double xl[6] 	= { 0.0, 	    -1.0, 	   0.0, 	    vLower,  -1.0, 	    0.0};
    double xu[6] 	= { 2*R_Earth, 	 1.0, 	   2.0*M_PI,	vEsc, 	 1.0, 	    2.0*M_PI}; // TEST: 2*R_Earth - 1 = 12741
  
    const gsl_rng_type *T;
    gsl_rng *r;
      
    //gcc -L/usr/lib iDM_ThreeDets_ATLAS.o -lgsl -lm -lgslcblas;
    gsl_monte_function G 	= { &sixDim_Integrand, 6, &params};

    size_t calls 		= pow(10,params.logprec);//10000000; // "max" is 8 zeros!
    //printf(" | | | | | | | | > Precision:    Log(calls) = %.1f   < | | | | | | | |\n\n", log10(calls));
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    //gsl_rng_set(r,0);
                
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
        
    //s->verbose = 0; /* default = -1 */
        s->alpha = 0.8; /* default = 1.5 */

    gsl_monte_vegas_integrate (&G, xl, xu, 6, calls, r, s, &res, &err);
        
    int totruns = 175;
        
    long double res_save[3] = {0.0,0.0,0.0};
   
    printf ("          ~~~~~~~~~~~~~~~~~~~~~ Vegas Warmup ~~~~~~~~~~~~~~~~~~~~~\n");
    printf ("          ~~~~~~~~~~~~~~~~~~~~~ alpha = %.1f ~~~~~~~~~~~~~~~~~~~~~\n",s->alpha);
    do
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 6, calls/5, r, s, &res, &err);
        printf ("          Run %d: result = % .6f | sigma = % .6f | ""chisq/dof = %.1f\n", runcount+1, res, err, s->chisq);
        runcount += 1;
        res_save[0] = res_save[1];
        res_save[1] = res_save[2];
        res_save[2] = res;
        
        //printf("------>  res = %Lf, %Lf, %Lf\n",res_save[0],res_save[1],res_save[2]);
    }
    while (fabs(s->chisq - 1.0) > 0.5 && runcount<=totruns);
    while (res_save[0] == 0.0  && runcount<=totruns);
      
    if(runcount >= totruns)
	{
        results[0] = 0.0;
        results[1] = .01;
    }
    else
	{
        results[0] = res;
        results[1] = err;
        gsl_monte_vegas_free (s);
    }
    gsl_rng_free (r);
    }

    return 0;
}