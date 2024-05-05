// !!! IMPORTANT  !!
//       Depending on the system you are using, you might use:
//  gcc -Wall -I/sw/include/ -c iDM.c; gcc -L/sw/lib iDM.o -lgsl -lm -lgslcblas; mv a.out iDM.out
//  ./iDM.out <det> <elt> <A> <gM> <mDM> <delta> <lifemult> <Month> <logprec>

//	OR
// gcc -Wall -I/usr/include -c iDM.c; gcc -L/usr/lib iDM.o -lgsl -lm -lgslcblas; mv a.out iDM.out
//  ./iDM.out <det> <elt> <A> <gM> <mDM> <delta> <lifemult> <Month> <logprec>

//       OR
// gcc -Wall -I/opt/local/include -c iDM.c; gcc -L/opt/local/lib iDM.o -lgsl -lm -lgslcblas; mv a.out iDM.out
//  ./iDM.out <det> <elt> <A> <gM> <mDM> <delta> <lifemult> <Month> <logprec>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_monte.h> 
#include <gsl/gsl_monte_plain.h> 
#include <gsl/gsl_monte_miser.h> 
#include <gsl/gsl_monte_vegas.h>
#include <time.h>

#include "Constants.c"
#include "Functions.h"

//-----------------------------------------------------------------------------------//
//   LOOPS FOR ONE DAY, COMPUTING ONCE PER HOUR, STARTING ON INPUTTED YEAR/MONTH/DAY
//-----------------------------------------------------------------------------------//

int main (int argc, char *argv[] ) // nVars = 0, vars = {lifetime}
{
    double results[2];
    struct param_type params;
  
	//--------------------------------------- Loop variables ------------------------------------------------
	
    int i_HOURS_MIN  = 0;  					// Start hour
    int i_HOURS_MAX  = 24; 					// End hour
    int i_HOURS_TOT  = 24; 					// Number of loops
    int i_HOURS_STEP = (i_HOURS_MAX-i_HOURS_MIN)/i_HOURS_TOT;
                            				// Step size, (End - Start)/(Number of loops)
	
	// -----------------> IF LOOPING OVER mDM : 
    //int count_mDM = 3; 											// Needs to match length in next line
    //double mDM_Range[3]	= {1000.0, 200.0, 500.0};
    //for(int i_mDM = 1; i_mDM <= count_mDM; i_mDM++) {
    //printf("\n\nLooping mDM: Run Number %d of %d\n", i_mDM, count_mDM);

	// -----------------> IF LOOPING OVER delta :
    //int count_delta = 4; 											// Needs to match length in next line
    //double delta_Range[4]	= {100.0, 125.0, 150.0, 175.0};
    //for(int i_delta = 1; i_delta <= count_delta; i_delta++) {
    //printf("\n\nLooping delta: Run Number %d of %d\n", i_delta, count_delta);

	// -----------------> IF LOOPING OVER lifemult :
	//int count_tau = 3; 											// Needs to match length in next line
    //double tau_range[3]	= {0.5  , 1.0  , 2.0};
    //for(int i_tau = 1; i_tau <= count_tau; i_tau++) {
    //printf("\n\nLooping lifemult: Run Number %d of %d\n", i_delta, count_delta);
	
	//-------------------------------- Read user-input variables -------------------------------------------

    //params.element[0] = 'F';//argv[1][0];
    //params.element[1] = 'e';//argv[1][1];
    params.detector = argv[1][0];
    strncpy(params.element, argv[2],5);
    params.Anum	  	= atof(argv[3]);
	params.gM		= atof(argv[4]);   								// Typical: 1.0 or 0.00633 = 1/16pi^2
    params.mDM	 	= atof(argv[5]);
	//params.mDM	= mDM_Range[i_mDM-1];							// IF LOOPING OVER mDM
   	params.delta	= atof(argv[6]); 			
	//params.delta	= delta_Range[i_delta-1];						// IF LOOPING OVER delta
   	params.lifemult = atof(argv[7]); //
	//params.lifemult = tau_range[i_tau-1];							// IF LOOPING OVER lifemult
    params.month    = atof(argv[8]); 
	params.logprec	= atof(argv[9]);
	
	params.day		= 1;
	params.year 	= 2000;
	//printf("date = %.0f/%.0f/%.0f\n",params.day,params.month,params.year);

	//----------------------------- Read isotope-specific variables ----------------------------------------
	
    double *p 		= read_isodata(params.element, params.Anum);
    double isodata[30];
    for (int i = 0; i <= 8; i++ ) isodata[i] = *(p+i);
  	params.Znum 	= isodata[0];
  	params.jN		= fabs(isodata[2]);
    params.Dcore 	= isodata[4]*isodata[5];
    params.Dmantle 	= isodata[4]*isodata[6];
    params.Dcrust 	= isodata[4]*isodata[7];
	
	double *lp;
  	for(int iOP=1; iOP<=5; iOP++)
	{
  		for(int iso1=0; iso1<=1; iso1++)
		{
  			for(int iso2=0; iso2<=1; iso2++)
			{
  				lp = getWcoefficients(params.Znum,params.Anum, iOP, iso1, iso2);
  				for (int i = 0; i <= 9; i++ ) params.Wcoeffs[iOP][iso1][iso2][i] = *(lp+i);
  			    //printf("Wcoeff = %d %d %f\n",iOP, iso1, params.Wcoeffs[iOP][iso1][iso2][0]);			
  			}
  		}
  	}	
		

	//for (int i=0; i<10; i++)	printf( "r[%d] = %f\n", i, W[i]);
	 /* 	
	printf("nucres = %Lf, %Lf, %Lf, %Lf, %Lf\n",
	NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 1, 0, 0, 1),
	NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 2, 0, 0, 1),
	NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 3, 0, 0, 1),
	NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 4, 0, 0, 1),
	NuclearResponse(params.Znum, params.Anum, params.Wcoeffs, 5, 0, 0, 1));
	  	
	//printf("DM response = %Lf\n", DMResponse(5, 1, 0, 1000, 200, 5, 500, 50));
	printf("DM response = %Lf\n %Lf \n %Lf \n %Lf \n %Lf \n", 
	DMResponse(params.mDM, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, 1, 0, 0),
	DMResponse(params.mDM, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, 2, 0, 0),
	DMResponse(params.mDM, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, 3, 0, 0),
	DMResponse(params.mDM, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, 4, 0, 0),
	DMResponse(params.mDM, params.delta*pow(10,-6)/cSqu,params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, 5, 0, 0));
	  */
	
		  
	//printf("MelSqu = %Lf\n",MatEltSqu(params.Znum, params.Anum, params.mDM, params.delta*pow(10,-6)/cSqu, params.Anum*0.93149/cSqu, params.delta*pow(10,-6), 600, params.Wcoeffs));
	
      time_t rawtime;
      struct tm * timeinfo;

      /* generate filename, and open the file */
	  
      FILE *fp;

      char filename[80];

      sprintf(filename, "MiDM_det=%c_%d%s_mDM=%.0Lf_delta=%.0Lf_tMult=%.1f_TIME=2000-%.0f-1_prec%.1f_gM=%.3f.txt", //.txt",
	  				params.detector, params.Anum, params.element, params.mDM, params.delta, params.lifemult, params.month,params.logprec,params.gM); //);

      printf("\nWrite to file: %s\n",filename);
      fp 					= fopen(filename, "w");

      /* write to the file */
	  
      //printf("-------------------------------------\n");
      printf("Start Date:    %.0f/%.0f/%.0f \n", params.day,params.month,params.year);
      printf("Hour Range:    %d to %d in steps of %d\n", i_HOURS_MIN,i_HOURS_MAX,i_HOURS_STEP);
	  printf("Precision:     10^%.1f\n",params.logprec);
	  printf("Coupling:      gM = %.5f\n",params.gM);	
	  // BUG! Today's time not appearing, why April 19??  

      printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTEGRATION BEGIN ~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	  fprintf(fp, "Run performed:   %s    \n",ctime(&rawtime));
	  fprintf(fp, "Precision level: 10^%.1f  ; gM = %.5f  ; alpha = 0.8 \n",params.logprec, params.gM);
	  fprintf(fp, "-------------------------------------------------------------------------\n");
      fprintf(fp, "(Year, Month,   Day,  Hour)        VEGAS                VEGAS err.\n");

      params.hour	  		= i_HOURS_MIN-i_HOURS_STEP;		// + 0.5; // params.hour = -1 + 0.5
    
      float int_time 		= 0;
      int time_start 		= 0;
      int time_end 			= 0;
      float tot_time 		= 0;
    
	//------------------------ Perform calculation for range of hours -----------------------------------
	
      for(int i_HOURS = i_HOURS_MIN; i_HOURS < i_HOURS_MAX; i_HOURS=i_HOURS+i_HOURS_STEP) 
	  {
          printf("\n | | | | | | | | >                 Run: Begin            < | | | | | | | | \n");
        
          // Start time
          time ( &rawtime );
          timeinfo 			= localtime ( &rawtime );
        
          time_start 		= 24*(timeinfo->tm_hour) + 60*(timeinfo->tm_min) + timeinfo->tm_sec;
        
          //printf(" | | | | | | | | > System Time: %d/%d/%d, %d:%d:%d < | | | | | | | | \n",timeinfo->tm_mday, timeinfo->tm_mon + 1, timeinfo->tm_year + 1900, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
          //printf ("                  System Time: %s", asctime (timeinfo) );
          //printf("%d\n",timeinfo->tm_sec);

          params.hour 		= params.hour + i_HOURS_STEP;
          printf(" | | | | | | | | >        Current Hour: h = %.2f         < | | | | | | | | \n\n", params.hour);
  		  MC_Ints(&params, results);
  		  display_results ("Vegas Final", results[0], results[1]);

  		  fprintf(fp, "2000,    %.0f,      1,   %.2f,    %.12f,       %.12f\n", 
		  					params.month, params.hour, results[0],results[1]);
        
          time ( &rawtime );
          timeinfo 			= localtime ( &rawtime );
          time_end 			= 24*(timeinfo->tm_hour) + 60*(timeinfo->tm_min) + timeinfo->tm_sec;
          int_time 			= time_end-time_start;
          if(int_time > 0) 
			  tot_time 		+= int_time;
		  
          printf("\n | | | | | | | | >                 Run: Finished         < | | | | | | | | \n");
          if(int_time < 0)
			  printf(" | | | | | | | | >    Integration Time: ERROR! (t<0)     < | | | | | | | | \n\n");
		  else if (int_time < 60)
			  printf(" | | | | | | | | >    Integration Time: %.2f seconds     < | | | | | | | | \n\n",int_time);
          else if(int_time < 3600)
			  printf(" | | | | | | | | >    Integration Time: %.2f minutes     < | | | | | | | | \n\n",int_time/60);
          printf("==========================================================================\n");
      	  printf("==========================================================================\n");

     }

    /* close the file, end the run */
	 
    fclose(fp);

    printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTEGRATION END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    
    printf("\n | | | | | | | | >  Total Integration Time: %.2f minutes  < | | | | | | | | \n\n",tot_time/60);
    
    return 0; 
}
