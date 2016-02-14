// 
// A Tempo2 plug-in to map the chi-squared space of pulsar timing data.
//
// Joris P.W. Verbiest, MPIfR Bonn, 28 March 2012

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "tempo2.h"
#include <fstream>
#include <cpgplot.h>
#include "T2toolkit.h"
#include <iostream>
#include <vector>
#include <iomanip>

#define DAY2S 86400
#define PARSEC 3.08568026e16

using namespace std;

int verb = 0; // generic verbosity flag

// Generic functions:
void Feedback( const char param[], char *parval );
void Feedback( int nread, const char param[], int parval, const char label[] );
void Feedback( int nread, const char param[], double parval, 
               const char label[] );
void Feedback( int nread, const char param[], long double parval, 
               const char label[] );
void GetGalCoord( char name[], long double *Gal );
long double GetKz( long double zz );
void GiveMeABreak();
void help();
void palett();
void PlotResiduals( pulsar *psr );
void Progress( time_t start, int counter, double size );
int ReadHeader( int form, char *fname, char *names[], char *psr, char *par,
                char *tim );
void SetOmDot( pulsar *psr );
int SetDist( pulsar *psr, long double dd );
void WriteResults( long int counter, long double val, long double xx, 
                   long double yy, long double zz, long double dd, int *map, 
                   int form, FILE *fout );

// Class for easy identification of parameters
class Param{
  char longlabel[100];
  char shortlabel[100];
  int ID;
  int derivID;
  long double min;
  long double max;
  long double sz;
  int nbins;
public:
  Param(){
    ID = -1;
    derivID = -1;
    min = -1.0L;
    max = -1.0L;
    sz = -1.0L;
    nbins = 1;
  }
  void Init( char *name, pulsar *psr ){
    // Initialise parameter
    strcpy( shortlabel, name );
    if( strcmp( shortlabel, "COSI" ) == 0 ){
      ID = -2;
      derivID = -2;
      strcpy( longlabel, "cos( i )" );
    }else if( strcmp( shortlabel, "H3" ) == 0 &&
              strcmp( psr[0].binaryModel, "DDH" ) != 0 &&
              strcmp( psr[0].binaryModel, "ELL1H" ) != 0 ){
      ID = -3; 
      derivID = -3;
      strcpy( longlabel, "h3" );
    }else if( strcmp( shortlabel, "DIST" ) == 0 ){
      ID = -2;
      derivID = -2;
      strcpy( longlabel, "Distance (kpc)" );
    }else{
      for( int p1ct = 0; p1ct < MAX_PARAMS; p1ct++ ){
        for( int p2ct = 0; p2ct < psr[0].param[p1ct].aSize; p2ct++ ){
          if( strcmp( shortlabel, psr[0].param[p1ct].shortlabel[p2ct] ) == 0 ){
            ID = p1ct;
            derivID = p2ct;
            strcpy( longlabel, psr[0].param[ID].label[derivID] );
          }
        }
      }
    }

    if( ID == -1 ){
      cout << "Did not understand parameter " << name << ".\n";
      cout << "Use option \"-all\" for a list of all known parameters.\n";
    }
    if( verb > 3 ){
      cout << "===============================\n"
           << "Found parameter: \n";
      cout << "\t longlabel: " << longlabel << endl
           << "\t shortlabel: " << shortlabel << endl
           << "\t ID: " << ID << endl
           << "\t deriv: " << derivID << endl
           << "\t min,max,sz: " << min << " " << max << " " << sz << endl
           << "\t nbins: " << nbins << endl; 
      cout << "===============================\n";
    }
  }// end of initialisation function.

  void Values( char *val1, char *val2, char *val3 ){
    // Reads values into min, sz and max.
    int nread; // checks if values are read in properly

    // First minimum:
    nread = sscanf( val1, "%Lg", &min );
    Feedback( nread, "minimum", min, "-param name min step max" );
    
    // Now step size:
    nread = sscanf( val2, "%Lg", &sz );
    Feedback( nread, "stepsize", sz, "-param name min step max" );

    // Now maximum:
    nread = sscanf( val3, "%Lg", &max );
    Feedback( nread, "maximum", max, "-param name min step max" );

    // Now calculate number of bins:
    nbins = (int)( ( max - min ) / sz + 0.5 ) + 1;
    nread = 1;
    Feedback( nread, "nbins", nbins, "-param name min step max" );
    if( nbins < 1 ){
      cout << "ERROR! For parameter " << shortlabel 
           << " you specified we should use " << nbins << " bins.\n"
           << "       This doesn't make sense! \n"
           << "       (Specified range was " << min << " to " << max << endl;
      cout << "EXITING...\n";
      fflush( stdout );
      exit( 1 );
    }
      
  }
  void Values( long double val1, long double val2, long double val3 ){
    // Reads values into min, sz and max.
    int nread; // checks if values are read in properly

    // First minimum:
    min = val1;
    
    // Now step size:
    sz = val2;

    // Now maximum:
    max = val3;

    // Now calculate number of bins:
    nbins = (int)( ( max - min ) / sz ) + 1;
  }
  int Size(){ return( nbins ); }
  long double I(){ return( min ); }
  long double step(){ return( sz ); }
  long double A(){ return( max ); }
  int par(){ return( ID ); }
  int der(){ return( derivID ); }
  char *name(){ return( shortlabel ); }
  char *fullname(){ return( longlabel ); }
}; // end of Param class definition

// Class for easy passing of input parameters
// Note that standard Tempo2 variables are not placed in this class
class InputParams{
  char dev[100];
  int parDef; // Counts how many parameters have been defined.
  char outfile[100];
  int resume; // Resume flag (0: start from scratch; 
              //              1: resume an interrupted run)
  int plotonly; // 0: also calculate; 1: only plot what's in the input file.
  int AscOut; // 0: binary output; 1: ascii output
  int RealPlot; // 1: no realtime plotting
                // *2: realtime map updating
                // *3: realtime residual plots
  //int RealPlot; // 0: no realtime plotting
  //              // 1: realtime map updating
  //              // 2: realtime residual plots
  //              // 3: both realtime plots
  vector<float> contlevels; // Contour levels
  int UseOmdot; // Whether or not to hardwire Omdot to its GR value
  int nits; // number of fitting iterations
public:
  InputParams( int argc, char **argv, Param *XYZ, pulsar *psr ){
    int ct; // generic counter
    // Initialise variables
    strcpy( dev, "11/xs" );
    parDef = 0;
    strcpy( outfile, "chisq.out" );
    resume = 0;
    plotonly = 0;
    AscOut = 0;
    RealPlot = 1;
    UseOmdot = 0;
    nits = 2;

    // Read command-line arguments
    int nread = 0; // To check if arguments are read correctly
    if( argc < 4 ){ // only used "tempo2 -gr ChiSqCube"
      help();
      exit( 0 );
    }

    for( ct = 1; ct < argc; ct++ ){ // Loop over all input arguments
      //cout << "\n\n\t JORIS: Reading " << argv[ct] << endl;
      // Verbosity level
      if( strcmp( argv[ct], "-v" ) *
          strcmp( argv[ct], "-verb" ) == 0 ){
        // Check if the next value exists and is not a new argument
        if( ct < argc - 1 && argv[ct+1][0] != '-' ){
          nread = sscanf( argv[++ct], "%d", &verb );
          Feedback( nread, "verbosity", verb, "-v" );
        }
      }

      // Contour levels
      else if( strcmp( argv[ct], "-cont" ) == 0 ){
        float tempcont;
        while( argv[ct+1][0] != '-'  && ct+1 < argc ){
          nread = sscanf( argv[++ct], "%g", &tempcont );
          Feedback( nread, "contour level", tempcont, "-cont" );
          contlevels.push_back( tempcont );
        }
      }
      
      // Use GR prediction for Omdot
      else if( strcmp( argv[ct], "-UseOmdot" ) == 0 ){
        UseOmdot = 1;
      }

      // Number of fitting iterations
      else if( strcmp( argv[ct], "-nits" ) == 0 ){
        nread = sscanf( argv[++ct], "%d", &nits );
        Feedback( nread, "Number of fitting iterations", nits, "-nits" );
      }

      // Realtime map updating
      else if( strcmp( argv[ct], "-PlotRealtime" ) == 0 ){
        //RealPlot += 1;
        RealPlot *= 2.0;
      }
      
      // Realtime residual plots
      else if( strcmp( argv[ct], "-ResidualPlot" ) == 0 ){
        //RealPlot += 2;
        RealPlot *= 3.0;
      }

      // Output: binary or ascii
      else if( strcmp( argv[ct], "-ascii" ) == 0 ){
        AscOut = 1;
      }

      // Graphics device
      else if( strcmp( argv[ct], "-dev" ) *
               strcmp( argv[ct], "-grdev" ) == 0 ){
        strcpy( dev, argv[++ct] );
        Feedback( "graphics device", dev );
      }

      // Output file
      else if( strcmp( argv[ct], "-file" ) *
               strcmp( argv[ct], "-out" ) == 0 ){
        strcpy( outfile, argv[++ct] );
        Feedback( "output file", outfile );
        
        // Test whether file already exists or not
        FILE *fin = fopen( outfile, "r" );
        if( fin == NULL ){
          // File does not yet exist. 
          // The default setting (resume = 0) applies.
        }else{
          // The file does already exist: we have to resume:
          resume = 1;
          // and now close the file
          fclose( fin );
        }
      }

      // Plot only
      else if( strcmp( argv[ct], "-plotonly" ) == 0 ){
        plotonly = 1;
      }

      // Parameters
      else if( strcmp( argv[ct], "-x" ) *
               strcmp( argv[ct], "-y" ) * 
               strcmp( argv[ct], "-z" ) *
               strcmp( argv[ct], "-param" ) == 0 ){
        XYZ[parDef].Init( argv[++ct], psr );
        Feedback( "parameter name", argv[ct] );
        //cout << "\n JORIS: " << argv[ct] << endl;
        XYZ[parDef].Values( argv[++ct], argv[ct+2], argv[ct+3] );
        ct+=2;
        parDef++;
      }
      else if( strcmp( argv[ct], "-params" ) == 0 ){
        XYZ[parDef].Init( argv[++ct], psr );
        Feedback( "parameter name", argv[ct] );
        XYZ[parDef].Values( argv[++ct], argv[ct+2], argv[ct+3] );
        ct+=2;
        parDef++;
        // Second parameter:
        XYZ[parDef].Init( argv[++ct], psr );
        Feedback( "parameter name", argv[ct] );
        XYZ[parDef].Values( argv[++ct], argv[ct+2], argv[ct+3] );
        ct+=2;
        parDef++;
        // Third parameter:
        XYZ[parDef].Init( argv[++ct], psr );
        Feedback( "parameter name", argv[ct] );
        XYZ[parDef].Values( argv[++ct], argv[ct+2], argv[ct+3] );
        ct+=2;
        parDef++;
      }

      // Print all known timing model parameters
      else if( strcmp( argv[ct], "-all" ) == 0 ){
        cout << "The known parameters are: " << endl;
        for( int p1ct = 0; p1ct < MAX_PARAMS; p1ct++ ){
          for( int p2ct = 0; p2ct < psr[0].param[p1ct].aSize; p2ct++ ){
            if( strcmp( psr[0].param[p1ct].shortlabel[p2ct], "" ) != 0 )
              cout << psr[0].param[p1ct].shortlabel[p2ct] << "\t"
                   << psr[0].param[p1ct].label[p2ct] << endl;
          }
        }
      }
               
    }// End loop over arguments
    
    // Check sensibility of input arguments
    if( (int)contlevels.size() == 0 ){
      // revert to defaults
      contlevels.push_back( 0.001 ); // optimum
      contlevels.push_back( 1.0 ); // 1 sigma
      contlevels.push_back( 4.0 ); // 2 sigma
    }

  } // End initialisation function

  // Change the name of the output file
  void Rename(){
    if( verb > 1 ){
      cout << "Changing filename to avoid overwriting...\n";
      fflush( stdout );
    }
    time_t rawtime;
    struct tm * timeinfo;
    time( &rawtime );
    timeinfo = localtime( &rawtime );
    char timeval[200];
    strcpy( timeval, asctime( timeinfo ) );
    strcpy( &timeval[(int)strlen(timeval)-1], "" );
    char time[5][20];
    sscanf( timeval, "%s %s %s %s %s", time[0], time[1], time[2], time[3],
            time[4] );
    char newfile[200];
    sprintf( newfile, "%s.%s%s%s%s-%s", outfile, time[0], time[2], time[1], 
             time[4], time[3] );
    strcpy( outfile, newfile );
  }

  // Returns the graphics device
  char *grdev(){ return( dev ); }

  // The name of the output file
  char *fout(){ return( outfile ); }

  // Wheter a run needs to be resumed (1) or started from scratch (0)
  int res(){ return( resume ); }
  // Override resume flag
  void SetRes( int val ){ resume = val; }

  // Number of contour lines
  int ncont(){ return( (int)contlevels.size() ); }
  // Contour levels
  void contl( float *contlev ){ 
    int clsize = (int)contlevels.size();
    for( int ct = 0; ct < clsize; ct++ )
      contlev[ct] = contlevels[ct];
  }

  // Using GR prediction for OMDOT
  int DoOm(){ return( UseOmdot ); }

  // Output format (0: binary; 1: ascii)
  int form(){ return( AscOut ); }

  // Toggle format
  void ReForm(){ 
    if( AscOut == 0 ) AscOut = 1;
    else if( AscOut == 1 ) AscOut = 0;
  }    

  // Number of fitting iterations
  int fitit(){ return( nits ); }
  
  // Only plotting (1) or also calculating (0)?
  int plot(){ return( plotonly ); }

  // Realtime plotting options
  int Realtime(){ return( RealPlot ); }

  // Number of parameters that are defined
  int pd(){ return( parDef ); }

}; // End InputParams class definition.

// Functions that use our custom classes:
void CalcCube( Param *XYZ, long double *XX, long double *YY, long double *ZZ,
               long double *DD, long double ****Map, pulsar *psr, int *npsr, 
               FILE *fout, InputParams IP, int parmap[], Param Dist );
void CheckRange( InputParams IP, Param *XYZ, pulsar *psr );
void FillArrays( long double *XX, long double *YY, long double *ZZ,                  
                 long double *DD, Param *XYZ, Param Dist );
void FindRange( long double *ar, int NN, const char *name, long double *par,
                float level );
void MakePlots( Param *XYZ, long double ****Map, InputParams IP, Param Dist);
void Merge( InputParams *IP, int argc, char *argv[] );
long double NormaliseMap( long double ****Map, Param *XYZ, Param Dist, 
                          int *location );
void OneDProject( Param *XYZ, Param Dist, long double ****Map, long double *XX,
                  long double *YY, long double *ZZ, long double *DD);
void PlotDist( Param *XYZ, long double ****Map, InputParams IP, Param Dist );
void PlotMapXY( Param *XYZ, long double ****Map, InputParams IP, Param Dist );
void PlotMapXZ( Param *XYZ, long double ****Map, InputParams IP, Param Dist );
void PlotMapYZ( Param *XYZ, long double ****Map, InputParams IP, Param Dist );
void PutHeader( FILE *fout, Param *XYZ, InputParams IP, pulsar *psr,
                char *parFile, char *timFile );
void Resume( Param *XYZ, Param Dist, InputParams *IP, FILE *fout, pulsar *psr,
             long double *XX, long double *YY, long double *ZZ, long double *DD,
             long double ****Map, char *parFile, char *timFile,
             int parmap[] );
void Transform( InputParams *IP );

extern "C" int graphicalInterface( int argc, char *argv[], pulsar *psr,
                                   int *npsr ){
  // =========================================
  // Declaration of Tempo2-required variables:
  // =========================================
  char parFile[1][MAX_FILELEN];
  char timFile[1][MAX_FILELEN];
  double globalParameter;
  *npsr = 1;
  int ct; // generic counter

  printf( "====================================================\n" );
  printf( "Graphical Interface: ChiSqCube\n" );
  printf( "Author:              Joris P.W. Verbiest\n" );
  printf( "Version Number:      0.0 - Under Construction\n" );
  printf( "Use \"tempo2 -gr ChiSqCube -h\" for help information.\n" );
  printf( "====================================================\n" );

  // ============================
  // Read Command-line arguments:
  // ============================
  // the array holding all information on the three parameters in question:
  Param XYZ[3];
  Param Dist;
  Dist.Init( "DIST", psr );
  Dist.Values( 0.1L, 0.01L, 2.0L );
  InputParams IP( argc, argv, XYZ, psr ); // Contains most input arguments
  
  // The Tempo2-required command-line arguments:
  for( ct = 1; ct < argc; ct++ ){ // Foreach argument

    // par and tim files
    if( strcmp( argv[ct], "-f" ) == 0 ){
      strcpy( parFile[0], argv[++ct] );
      Feedback( "parfile name", parFile[0] );
      strcpy( timFile[0], argv[++ct] );
      Feedback( "timfile name", timFile[0] );
    }
    // Help
    else if( strcmp( argv[ct], "-h" ) * strcmp( argv[ct], "--help" ) == 0 ){
      help();
      exit( 0 );
    }
    
    // Distance range
    else if( strcmp( argv[ct], "-dist" ) == 0 ){
      Dist.Values( argv[ct+1], argv[ct+2], argv[ct+3] );
      if( verb ) 
        cout << "Read distance range [" << Dist.I() << "::" 
             << Dist.A() << "] with stepsize " << Dist.step() << endl;
    }

    // Transform binary <--> ascii
    else if( strcmp( argv[ct], "-transform" ) == 0 ){
      Transform( &IP );
      exit( 0 );
    }

    // Merge files
    else if( strcmp( argv[ct], "-merge" ) == 0 ){
      Merge( &IP, argc, argv );
      exit( 0 );
    }

  } // end of loop over all arguments

  // Initial execution of the basic functions
  if( verb > 10 )
    cout << "Initial fitting...\n"; fflush( stdout );

  readParfile( psr, parFile, timFile, *npsr );
  readTimfile( psr, timFile, *npsr );
  preProcess( psr, *npsr, argc, argv );
  if( IP.plot() == 0 ){
    for( ct = 0; ct < 2; ct++ ){
      formBatsAll( psr, *npsr );
      formResiduals( psr, *npsr, 1 );
      if( ct == 0 )
        doFitAll( psr, *npsr, NULL );
      else 
        calcRMS( psr, 0 );
    }
  }
  if( strcmp( psr[0].binaryModel, "DDH" ) != 0 &&
      strcmp( psr[0].binaryModel, "ELL1H" ) != 0 ){
    // GD CheckRange( IP, XYZ, psr );
    // We probably should include a check for the DDH and ELL1H models
    // as well, but I haven't gotten around to that yet.
  }
  
  // ===========================================================================
  // ACTUAL START OF THE PLUG-IN
  // ===========================================================================
  // 
  // This plug-in attempts to be the final answer to all our chi-squared mapping
  // worries. 
  // 
  // This means I will try to make things as flexible as possible.
  // ===========================================================================

  // Step 1: Determine grid for chi-squared evaluation

  // Determination of arrays to hold x,y,z values and cube data.
  // Rather than determine an actual three-dimensional cube, we will
  // work with one-dimensional arrays. The main advantage to this is
  // that such 1D arrays are required by pgplot for plotting maps,
  // ergo, we will not need to resample the cubes into arrays.
  if( verb > 10 )
    cout << "Allocating memory...\n"; fflush( stdout );

  long double *XX;
  long double *YY;
  long double *ZZ;
  long double *DD;
  long double ****Map;
  //int MapSize = XYZ[0].Size() * XYZ[1].Size() * XYZ[2].Size();
  XX = (long double *)malloc( XYZ[0].Size() * sizeof( long double ) );
  YY = (long double *)malloc( XYZ[1].Size() * sizeof( long double ) );
  ZZ = (long double *)malloc( XYZ[2].Size() * sizeof( long double ) );
  DD = (long double *)malloc( Dist.Size() * sizeof( long double ) );

  Map = (long double ****)malloc( XYZ[0].Size() * sizeof( long double ***) );
  for( ct = 0; ct < XYZ[0].Size(); ct++ ){
    Map[ct] = (long double ***)malloc( XYZ[1].Size() * sizeof( long double **) );
    for( int ct2 = 0; ct2 < XYZ[1].Size(); ct2++ ){
      Map[ct][ct2] = (long double **)malloc(  XYZ[2].Size() *
                                             sizeof( long double *) );
      for( int ct3 = 0; ct3 < XYZ[2].Size(); ct3++ ){
        Map[ct][ct2][ct3] = (long double *)malloc( Dist.Size() * 
                                                   sizeof( long double ) );
        for( int ct4 = 0; ct4 < Dist.Size(); ct4++ ){
          Map[ct][ct2][ct3][ct4] = -1.0;
        }
      }
    }
  }

  // Initialise the coordinate arrays:
  if( verb > 10 )
    cout << "Filling arrays...\n"; fflush( stdout );
  FillArrays( XX, YY, ZZ, DD, XYZ, Dist );

  // Step 2: Ensure sensible fitting choices (parameters, jumps)
  // ===========================================================
  // Turn off fitting for the parameters of choice
  psr[0].param[param_px].fitFlag[0] = 0;
  if( verb > 10 ){
    cout << "Turning off fitting for the parameters. \n";
    cout << "Specifically: " << XYZ[0].par() << " " << XYZ[0].der()
         << " and " << XYZ[1].par() << " " << XYZ[1].der() 
         << " and " << XYZ[2].par() << " " << XYZ[2].der() << endl;
    fflush( stdout );
  }
  for( ct = 0; ct < 3; ct++ ){
    if( XYZ[ct].par() != -2 && XYZ[ct].par() != -3 )
      psr[0].param[XYZ[ct].par()].fitFlag[XYZ[ct].der()] = 0;
    else if( XYZ[ct].par() == -2 ){
      // Have cos(i). Turn off fitting for SINI and KIN:
      psr[0].param[param_sini].fitFlag[0] = 0;
      psr[0].param[param_kin].fitFlag[0] = 0;
    }else if( XYZ[ct].par() == -3 ){
      // Have H3. Turn off fitting for M2:
      psr[0].param[param_m2].fitFlag[0] = 0;
    }else{
      cout << "ERROR! Don't know parameter with ID " << XYZ[ct].par()
           << " this will probably break, but we'll run with it anyway.\n";
      fflush( stdout );
    }
  }
  if( verb > 100 )
    cout << "Finished checking fit flags. Moving on.\n"; fflush( stdout );

  // Step 3a: Read input to build on.
  FILE *fout;
  // Initialise the output file
  int parmap[3]; // maps the command-line parameters to the input file params.
  parmap[0] = 0;
  parmap[1] = 1;
  parmap[2] = 2;

  if( IP.res() == 1 ){
    // Resume from file
    if( verb > 1 )
      cout << "Resuming from existing file...\n"; fflush( stdout );
    Resume( XYZ, Dist, &IP, fout, psr, XX, YY, ZZ, DD, Map, parFile[0], 
            timFile[0],parmap );
    if( IP.form() == 0 )
      fout = fopen( IP.fout(), "ab" );
    else
      fout = fopen( IP.fout(), "a" );
    if( fout != NULL && verb > 0 ){
      cout << "Opened output file (0).\n"; fflush( stdout );
    }

  }else{ 
    // Don't resume: start a new file.
    if( verb > 1 )
      cout << "Opening new output file.\n"; fflush( stdout );
    if( IP.form() == 0 )
      fout = fopen( IP.fout(), "wb" );
    else
      fout = fopen( IP.fout(), "w" );
    if( fout != NULL && verb > 0 ){
      cout << "Opened output file (1).\n"; fflush( stdout );
    }
    PutHeader( fout, XYZ, IP, psr, parFile[0], timFile[0] );
    if( verb > 1 )
      cout << "Put header.\n"; fflush( stdout );
  }
    
  // Step 3: Loop over grid to evaluate post-fit chi-squared values
  // ==============================================================
  if( IP.plot() == 0 )
    CalcCube( XYZ, XX, YY, ZZ, DD, Map, psr, npsr, fout, IP, parmap, Dist );
          

  // Step 4: Output results (interactive plotting/projecting?)
  int loc[4];
  cout << "Going to normalise...\n"; fflush( stdout );
  long double min = NormaliseMap( Map, XYZ, Dist, loc );
  if( verb ){
    cout << "Minimum Chisquared value: " << min 
         << " at " << XYZ[0].name() << " = " << XX[loc[0]]
         << "; " << XYZ[1].name() << " = " << YY[loc[1]]
         << " and " << XYZ[2].name() << " = " << ZZ[loc[2]] 
         << " Distance = " << DD[loc[3]] << endl;
    fflush( stdout );
  }

  OneDProject( XYZ, Dist, Map, XX, YY, ZZ, DD );

  MakePlots( XYZ, Map, IP, Dist );

  free( XX ), free( YY ), free( ZZ ), free( DD );
  for( ct = 0; ct < XYZ[0].Size(); ct++ ){
    for( int ct2 = 0; ct2 < XYZ[1].Size(); ct2++ ){
      for( int ct3 = 0; ct3 < XYZ[2].Size(); ct3++ )
        free(Map[ct][ct2][ct3]);
      free( Map[ct][ct2] );
    }
    free( Map[ct] );
  }
  free( Map );
    
  return( 0 );
}// end of graphicalInterface function

void OneDProject( Param *XYZ, Param Dist, long double ****Map, long double *XX,
                  long double *YY, long double *ZZ, long double *DD ){
  // First in X dimension:
  int xct, yct, zct, dct;

  long double *Xmin;
  Xmin = (long double *)malloc( XYZ[0].Size() * sizeof( long double ) );

  for( xct = 0; xct < XYZ[0].Size(); xct++ ){
    Xmin[xct] = Map[xct][0][0][0];
    for( yct = 0; yct < XYZ[1].Size(); yct++ )
      for( zct = 0; zct < XYZ[2].Size(); zct++ )
        for( dct = 0; dct < Dist.Size(); dct++ )
          if( Xmin[xct] == -1.0 || ( Map[xct][yct][zct][dct] < Xmin[xct] &&
                                     Map[xct][yct][zct][dct] != -1.0 ) )
            Xmin[xct] = Map[xct][yct][zct][dct];
  }
  FindRange( Xmin, XYZ[0].Size(), XYZ[0].name(), XX, 1.0 );
  FindRange( Xmin, XYZ[0].Size(), XYZ[0].name(), XX, 4.0 );
  free( Xmin );

  long double *Ymin;
  Ymin = (long double *)malloc( XYZ[1].Size() * sizeof( long double ) );

  for( yct = 0; yct < XYZ[1].Size(); yct++ ){
    Ymin[yct] = Map[0][yct][0][0];
    for( xct = 0; xct < XYZ[0].Size(); xct++ )
      for( zct = 0; zct < XYZ[2].Size(); zct++ )
        for( dct = 0; dct < Dist.Size(); dct++ )
          if( Ymin[yct] == -1.0 || ( Map[xct][yct][zct][dct] < Ymin[yct] &&
                                     Map[xct][yct][zct][dct] != -1.0 ) )
            Ymin[yct] = Map[xct][yct][zct][dct];
  }
  FindRange( Ymin, XYZ[1].Size(), XYZ[1].name(), YY, 1.0 );
  FindRange( Ymin, XYZ[1].Size(), XYZ[1].name(), YY, 4.0 );
  free( Ymin );

  long double *Zmin;
  Zmin = (long double *)malloc( XYZ[2].Size() * sizeof( long double ) );
  
  for( zct = 0; zct < XYZ[2].Size(); zct++ ){
    Zmin[zct] = Map[0][0][zct][0];
    for( xct = 0; xct < XYZ[0].Size(); xct++ )
      for( yct = 0; yct < XYZ[1].Size(); yct++ )
        for( dct = 0; dct < Dist.Size(); dct++ )
          if( Zmin[zct] == -1.0 || ( Map[xct][yct][zct][dct] < Zmin[zct] &&
                                     Map[xct][yct][zct][dct] != -1.0 ) )
          Zmin[zct] = Map[xct][yct][zct][dct];
  }
  FindRange( Zmin, XYZ[2].Size(), XYZ[2].name(), ZZ, 1.0 );
  FindRange( Zmin, XYZ[2].Size(), XYZ[2].name(), ZZ, 4.0 );
  free( Zmin );

  long double *Dmin;
  Dmin = (long double *)malloc( Dist.Size() * sizeof( long double ) );
  
  for( dct = 0; dct < Dist.Size(); dct++ ){
    Dmin[dct] = Map[0][0][0][dct];
    for( xct = 0; xct < XYZ[0].Size(); xct++ )
      for( yct = 0; yct < XYZ[1].Size(); yct++ )
        for( zct = 0; zct < XYZ[2].Size(); zct++ )
          if( Dmin[dct] == -1.0 || 
              ( Map[xct][yct][zct][dct] < Dmin[dct] && 
                Map[xct][yct][zct][dct] != -1.0 ) )
            Dmin[dct] = Map[xct][yct][zct][dct];
  }
  FindRange( Dmin, Dist.Size(), Dist.name(), DD, 1.0 );
  FindRange( Dmin, Dist.Size(), Dist.name(), DD, 4.0 );
  free( Dmin );
}

void FindRange( long double *ar, int NN, const char *name, long double *par,
                float level ){
  int ct;
  int id1 = -1, id2 = -1, id3 = -1, id4 = -1;
  int idmin = -1;

  for( ct = 0; ct < NN; ct++ )
    if( (float)ar[ct] == 0.0 )
      idmin = ct;
  if( idmin == -1 ){
    cout << "ERROR! No minimum found!\n";
    cout << "Exiting... \n";
    fflush( stdout );
  }

  // Find the left edge of the 1-sigma interval:
  ct = 0; 
  while( (float)ar[ct] > level && ct < NN )
    ct++;
  if( ct < NN ){
    // Have entered minimum
    id1 = ct;
  }else{
    cout << "RESULT: " << name << " min = " << par[idmin] 
         << " Entire range (" << par[0] << " :: " << par[NN-1] 
         << ") within " << sqrt( level ) << "-sigma range.\n";
    fflush( stdout );
  }


  while( (float)ar[ct] < level && ct < NN )
    ct++;
  if( ct < NN ){
    id2 = ct;
  }else{
    id2 = NN - 1;
    cout << "RESULT: " << name << " min = " << par[idmin]
         << " " << sqrt( level ) << "-sigma interval: " 
         << par[id1] << " :: " << par[NN-1] << endl;
    fflush( stdout );
  }  

  // Look for second minimum
  while( (float)ar[ct] > level && ct < NN )
    ct++;
  
  if( ct < NN ){
    // Have found second minimum!
    id3 = ct;
    
    while( (float)ar[ct] < level && ct < NN )
      ct++;
    
    if( ct < NN ){
      id4 = ct;
    }else{
      id4 = NN - 1;
    }
    cout << "RESULT: " << name <<  " min = " << par[idmin]
         << " " << sqrt( level ) << "-sigma interval: " 
         << par[id1] << "::" << par[id2] << " and " 
         << par[id3] << "::" << par[id4] << endl;
    fflush( stdout );
  }else{
    cout << "RESULT: " << name << " min = " << par[idmin]
         << " " << sqrt( level ) << "-sigma interval: " 
         << par[id1] << "::" << par[id2]
         << endl;
    fflush( stdout );
  }
  cout << endl;
}

void MakePlots( Param *XYZ, long double ****Map, InputParams IP, Param Dist ){
  
  cpgbeg( 0, IP.grdev(), 2, 2 );
  cpgpap( 0.0, 1.0 );
  cpgscir( 0, 15 );
  cpgsch( 1.0 );

  PlotMapXY( XYZ, Map, IP, Dist );
  PlotMapXZ( XYZ, Map, IP, Dist );
  PlotMapYZ( XYZ, Map, IP, Dist );
  PlotDist( XYZ, Map, IP, Dist ); 
  cpgend( );
}

void PlotDist( Param *XYZ, long double ****Map, InputParams IP, Param Dist ){
  int ct, xct, yct, zct, dct, mapct; // counters
  float *DMap;
  DMap = (float *)malloc( Dist.Size() * sizeof( float ) );
  float *Dvals;
  Dvals = (float *)malloc( Dist.Size() * sizeof( float ) );
  
  long double min;
  long double globalmin = -1.0;
  float max = -1.0;

  mapct = 0;
  for( dct = 0; dct < Dist.Size(); dct++ ){
    min = Map[0][0][0][dct];
    for( xct = 0; xct < XYZ[0].Size(); xct++ )
      for( yct = 0; yct < XYZ[1].Size(); yct++ )
        for( zct = 0; zct < XYZ[2].Size(); zct++ ){
          if( min == -1.0 ||
              (Map[xct][yct][zct][dct] < min && 
               Map[xct][yct][zct][dct] != -1.0 ) )
            min = Map[xct][yct][zct][dct];
        }
    DMap[dct] = min;
    Dvals[dct] = ( Dist.I() + (long double)dct * Dist.step() );
    if( min > max || max == -1.0 )
      max = min;
    if( globalmin == -1.0 || ( globalmin > min && min != -1.0 ) )
      globalmin = min;
  }// End loop over D.

  // Put uncalculated pixels equal to the maximum and subtract the minimum:
  //cout << endl << "JORIS: : " << DMap[0] << endl; fflush( stdout );
  for( ct = 0; ct < Dist.Size(); ct++ )
    if( DMap[ct] == -1.0 )
      DMap[ct] = max - globalmin;
    else
      DMap[ct] -= globalmin;
  max -= globalmin;
  //cout << "\tJORIS: : " << DMap[0] << " " << DMap[1] << endl; fflush( stdout );

  // Plot curve:
  char xstr[100], ystr[100], titstr[100];
  sprintf( xstr, "%s", Dist.fullname() );
  sprintf( ystr, "\\gD\\gX\\u2\\d" );
  sprintf( titstr, "Goodness of fit curve for distance" );

  cpgsci( 15 );
  cpgenv( (float)Dist.I(), (float)Dist.A(), 0.0, max, 0, 0 );
  cpglab( xstr, ystr, titstr );
  cpgline( Dist.Size(), Dvals, DMap );
  //cout << "Plotting: " << Dist.Size() << " values.\n"; fflush( stdout );
  //for( ct = 0; ct < (int)Dist.Size(); ct++ )
  //  cout << Dvals[ct] << " : " << DMap[ct] << endl; fflush( stdout );
             
}

void PlotMapXY( Param *XYZ, long double ****Map, InputParams IP, Param Dist ){
  int ct, xct, yct, zct, dct, mapct; // counters
  float *XYMap;
  XYMap = (float *)malloc( XYZ[0].Size() * XYZ[1].Size() * sizeof( float ) );

  long double min;
  long double globalmin = -1.0;
  float max = -1.0;
  
  mapct = 0;
  for( yct = 0; yct < XYZ[1].Size(); yct++ ){
    for( xct = 0; xct < XYZ[0].Size(); xct++ ){
      min = Map[xct][yct][0][0];
      for( zct = 0; zct < XYZ[2].Size(); zct++ )
        for( dct = 0; dct < Dist.Size(); dct++ ){
          if( min == -1.0 || 
              ( Map[xct][yct][zct][dct] < min && 
                Map[xct][yct][zct][dct] != -1.0 ) )
            min = Map[xct][yct][zct][dct];
        }// end loop over d
      XYMap[mapct] = min;
      if( min > max || max == -1.0 )
        max = min;
      if( globalmin == -1.0 || ( globalmin > min && min != -1.0 ) ) 
        globalmin = min;
      mapct++;
    }// end loop over x
  } // end loop over y
  
  // Put uncalculated pixels equal to the maximum:
  for( ct = 0; ct < mapct; ct++ )
    if( XYMap[ct] == -1.0 )
      XYMap[ct] = max;

  for( ct = 0; ct < mapct; ct++ )
    XYMap[ct] -= globalmin;
  max -= globalmin;

  // Plot map:
  // pgplot transformation matrix
  float tr[6];
  tr[0] = (float)(XYZ[0].I() - 0.5 * (XYZ[0].A()-XYZ[0].I())/XYZ[0].Size());
  tr[1] = (float)((XYZ[0].A() - XYZ[0].I())/XYZ[0].Size());
  tr[2] = 0.0;
  tr[3] = (float)(XYZ[1].I()-0.5*(XYZ[1].A()-XYZ[1].I())/XYZ[1].Size() );
  tr[4] = 0.0;
  tr[5] = (float)((XYZ[1].A()-XYZ[1].I())/XYZ[1].Size());

  char xstr[100], ystr[100], titstr[100];
  sprintf( xstr, "%s", XYZ[0].name() );
  sprintf( ystr, "%s", XYZ[1].name() );
  sprintf( titstr, "\\gD\\gx\\u2\\d map for %s and %s",
           XYZ[0].name(), XYZ[1].name() );

  cpgsci( 15 );
  cpgenv( (float)XYZ[0].I(), (float)XYZ[0].A(), 
          (float)XYZ[1].I(), (float)XYZ[1].A(), 0, 1 );
  cpglab( xstr, ystr, titstr );
  palett();
  cpgimag( XYMap, XYZ[0].Size(), XYZ[1].Size(), 
           1, XYZ[0].Size(), 1, XYZ[1].Size(), 
           0.0, max, tr );

  cpgscr( 16, 1.0, 0.0, 0.0 );
  cpgsci( 16 );
  cpgslw( 6 );

  float contLevels[IP.ncont()];
  IP.contl( contLevels );
  cpgcont( XYMap, XYZ[0].Size(), XYZ[1].Size(), 1, XYZ[0].Size(), 1, 
           XYZ[1].Size(), contLevels, IP.ncont(), tr );
  cpgslw( 1 );
  cpgsci( 15 );

  free( XYMap );
}

void PlotMapXZ( Param *XYZ, long double ****Map, InputParams IP, Param Dist ){
  int ct, xct, yct, zct, dct, mapct; // counters
  float *XZMap;
  XZMap = (float *)malloc( XYZ[0].Size() * XYZ[2].Size() * sizeof( float ) );

  long double min;
  long double globalmin = -1.0;
  float max = -1.0;

  mapct = 0;
  for( zct = 0; zct < XYZ[2].Size(); zct++ ){
    for( xct = 0; xct < XYZ[0].Size(); xct++ ){
      min = Map[xct][0][zct][0];
      for( yct = 0; yct < XYZ[1].Size(); yct++ )
        for( dct = 0; dct < Dist.Size(); dct++ ){
          if( min == -1.0 || 
              ( Map[xct][yct][zct][dct] < min  && 
                Map[xct][yct][zct][dct] != -1.0 ) )
            min = Map[xct][yct][zct][dct];
        }// end loop over d
      XZMap[mapct] = min;
      if( min > max || max == -1.0 )
        max = min;
      if( globalmin == -1.0 || ( globalmin > min && min != -1.0 ) )
        globalmin = min;
      mapct++;
    }// end loop over x
  }// end loop over z

  // Put uncalculated pixels equal to the maximum:
  for( ct = 0; ct < mapct; ct++ )
    if( XZMap[ct] == -1.0 )
      XZMap[ct] = max;

  for( ct = 0; ct < mapct; ct++ )
    XZMap[ct] -= globalmin;
  max -= globalmin;

  // Plot map:
  // pgplot transformation matrix
  float tr[6];
  tr[0] = (float)( XYZ[0].I() - 0.5 * 
                   ( XYZ[0].A() - XYZ[0].I() ) / XYZ[0].Size() );
  tr[1] = (float)( ( XYZ[0].A() - XYZ[0].I() ) / XYZ[0].Size() );
  tr[2] = 0.0;

  tr[3] = (float)( XYZ[2].I() - 0.5 * 
                   ( XYZ[2].A() - XYZ[2].I() ) / XYZ[2].Size() );
  tr[4] = 0.0;
  tr[5] = (float)( ( XYZ[2].A() - XYZ[2].I() ) / XYZ[2].Size() );

  cpgsci( 15 );
  cpgenv( (float)XYZ[0].I(), (float)XYZ[0].A(),
          (float)XYZ[2].I(), (float)XYZ[2].A(), 0, 1 );

  // axis labels and title:
  char xstr[100], zstr[100], titstr[100];
  sprintf( xstr, "%s", XYZ[0].name() );
  sprintf( zstr, "%s", XYZ[2].name() );
  sprintf( titstr, "\\gD\\gx\\u2\\d map for %s and %s",
           XYZ[0].name(), XYZ[2].name() );
  cpglab( xstr, zstr, titstr );
  palett();
  cpgimag( XZMap, XYZ[0].Size(), XYZ[2].Size(),
           1, XYZ[0].Size(), 1, XYZ[2].Size(), 0.0, max, tr );

  cpgscr( 16, 1.0, 0.0, 0.0 );
  cpgsci( 16 );
  cpgslw( 6 );

  float contLevels[IP.ncont()];
  IP.contl( contLevels );
  cpgcont( XZMap, XYZ[0].Size(), XYZ[2].Size(), 1, XYZ[0].Size(), 1, 
           XYZ[2].Size(), contLevels, IP.ncont(), tr );
  cpgslw( 1 );
  cpgsci( 15 );

  free( XZMap );
}

void PlotMapYZ( Param *XYZ, long double ****Map, InputParams IP, Param Dist ){
  int ct, xct, yct, zct, dct, mapct; // counters
  float *YZMap;
  YZMap = (float *)malloc( XYZ[1].Size() * XYZ[2].Size() * sizeof( float ) );
  long double min;
  long double globalmin = -1.0;
  float max = -1.0;
  
  mapct = 0;
  for( zct = 0; zct < XYZ[2].Size(); zct++ ){
    for( yct = 0; yct < XYZ[1].Size(); yct++ ){
      min = Map[0][yct][zct][0];
      for( xct = 0; xct < XYZ[0].Size(); xct++ )
        for( dct = 0; dct < Dist.Size(); dct++ ){
          if( min == -1.0 ||
              ( Map[xct][yct][zct][dct] < min && 
                Map[xct][yct][zct][dct] != -1.0 ) )
          min = Map[xct][yct][zct][dct];
      }// end loop over d
      YZMap[mapct] = min;
      if( min > max )
        max = min;
      if( globalmin == -1.0 || ( globalmin > min && min != -1.0 ) )
        globalmin = min;
      mapct++;
    }// end loop over y
  }// end loop over z

  // Put uncalculated pixels equal to the maximum:
  for( ct = 0; ct < mapct; ct++ )
    if( YZMap[ct] == -1.0 )
      YZMap[ct] = max;

  for( ct = 0; ct < mapct; ct++ )
    YZMap[ct] -= globalmin;
  max -= globalmin;

  // Plot map:
  // pgplot transformation matrix
  float tr[6];
  tr[0] = (float)( XYZ[1].I() - 0.5 *
                   ( XYZ[1].A() - XYZ[1].I() ) / XYZ[1].Size() );
  tr[1] = (float)( ( XYZ[1].A() - XYZ[1].I() ) / XYZ[1].Size() );
  tr[2] = 0.0;

  tr[3] = (float)( XYZ[2].I() - 0.5 * 
                   ( XYZ[2].A() - XYZ[2].I() ) / XYZ[2].Size() );
  tr[4] = 0.0;
  tr[5] = (float)( ( XYZ[2].A() - XYZ[2].I() ) / XYZ[2].Size() );

  cpgsci( 15 );
  cpgenv( (float)XYZ[1].I(), (float)XYZ[1].A(),
          (float)XYZ[2].I(), (float)XYZ[2].A(), 0, 1 );
  
  // axis labels and title:
  char ystr[100], zstr[100], titstr[100];
  sprintf( ystr, "%s", XYZ[1].name() );
  sprintf( zstr, "%s", XYZ[2].name() );
  sprintf( titstr, "\\gD\\gx\\u2\\d map for %s and %s", 
           XYZ[1].name(), XYZ[2].name() );
  cpglab( ystr, zstr, titstr );
  palett();
  cpgimag( YZMap, XYZ[1].Size(), XYZ[2].Size(), 
           1, XYZ[1].Size(), 1, XYZ[2].Size(), 0.0, max, tr );
  
  cpgscr( 16, 1.0, 0.0, 0.0 );
  cpgsci( 16 );
  cpgslw( 6 );

  float contLevels[IP.ncont()];
  IP.contl( contLevels );
  cpgcont( YZMap, XYZ[1].Size(), XYZ[2].Size(), 1, XYZ[1].Size(), 1, 
           XYZ[2].Size(), contLevels, IP.ncont(), tr );
  cpgslw( 1 );
  cpgsci( 15 );
  
  free( YZMap );
}

long double NormaliseMap( long double ****Map, Param *XYZ, Param Dist,
                          int *location ){
  long int xct, yct, zct, dct; // generic counters
  long double minimum = -1.0;
  long double maximum = -1.0;
  location[0] = 0;
  location[1] = 0;
  location[2] = 0;
  location[3] = 0;

  // Find minimum
  for( xct = 0; xct < XYZ[0].Size(); xct++ )
    for( yct = 0; yct < XYZ[1].Size(); yct++ )
      for( zct = 0; zct < XYZ[2].Size(); zct++ )
        for( dct = 0; dct < Dist.Size(); dct++ ){
          if( minimum == -1.0 || 
              ( Map[xct][yct][zct][dct] < minimum 
                && Map[xct][yct][zct][dct] != -1.0 ) ){
            minimum = Map[xct][yct][zct][dct];
            location[0] = xct;
            location[1] = yct;
            location[2] = zct;
            location[3] = dct;
          }
          if( Map[xct][yct][zct][dct] > maximum )
            maximum = Map[xct][yct][zct][dct];
        }
  cout << "Minimum: " << minimum << endl;
  if( verb ) cout << "Maximum: " << maximum << endl; fflush( stdout );

  // Subtract minimum
  for( xct = 0; xct < XYZ[0].Size(); xct++ )
    for( yct = 0; yct < XYZ[1].Size(); yct++ )
      for( zct = 0; zct < XYZ[2].Size(); zct++ )
        for( dct = 0; dct < Dist.Size(); dct++ )
          if( Map[xct][yct][zct][dct] != -1.0 )
            Map[xct][yct][zct][dct] -= minimum;
  
  return( minimum );
}

void SetParam( pulsar *psr, Param XYZ, long double val ){

  if( XYZ.par() != -2 && XYZ.par() != -3 ){
    psr[0].param[XYZ.par()].val[XYZ.der()] = val;
    if( strcmp(XYZ.name(),"KIN") == 0 ){
      // Working with KIN. Should also set sini!
      psr[0].param[param_sini].val[0] = sin( val/180.0*M_PI );
    }
  }else if( XYZ.par() == -2 ){
    // Working with cos(i)
    psr[0].param[param_sini].val[0] = sin( acos( val ) );
    if( psr[0].param[param_kin].paramSet[0] == 1 )
      psr[0].param[param_kin].val[0] = acos( val ) / M_PI * 180.0;
  }else if( XYZ.par() == -3 ){
    // Working with H3 without DDH or ELL1H model
    // First calculate STIG:
    long double stig;
    stig = psr[0].param[param_sini].val[0] /
      ( 1.0 + sqrt( 1.0 - pow( psr[0].param[param_sini].val[0], 2.0 ) ) );
    // Now define m2 based on the desired H3 and stigma:
    psr[0].param[param_m2].val[0] = val / pow( stig, 3.0 ) / TSUN;
    if( verb > 100 )
      cout << "Putting H3 value " << val << " as M2 = " << 
        psr[0].param[param_m2].val[0] << endl;
  }else{
    cout << "ERROR: Don't know which parameter has index " << XYZ.par() 
        << ". This will probably break the program.\n";
  }
}

void CalcCube( Param *XYZ, long double *XX, long double *YY, long double *ZZ,
               long double *DD, long double ****Map, pulsar *psr, int *npsr, 
               FILE *fout, InputParams IP, int parmap[], Param Dist ){
  if( verb )
    cout << "Calculating cube...\n"; fflush( stdout );
  int rv; // return value GD
  int ct; // generic counter
  int xct; // counter for x coordinate
  int yct; // counter for y coordinate
  int zct; // counter for z coordinate
  int dct; // counter for distance coordinate
  long int counter = 0;
  long int calculated = 0;

  time_t start;
  double size = (double)( XYZ[0].Size() * XYZ[1].Size() * XYZ[2].Size()
                          * Dist.Size() );
  time( &start );

  // Save timing model for reverting when needed:
  vector<int> parID, derID;
  vector<long double> parvals;
  for( ct = 0; ct < MAX_PARAMS; ct++ )
    for( int ct2 = 0; ct2 < psr[0].param[ct].aSize; ct2++ )
      if( psr[0].param[ct].fitFlag[ct2] != 0.0 ){
        parID.push_back( ct );
        derID.push_back( ct2 );
        parvals.push_back( psr[0].param[ct].val[ct2] );
      }
  if( verb > 5 ) 
    cout << "Saved timing model.\n"; fflush( stdout );

  if( verb > 10 ){
    cout << "Starting iteration over: " << XYZ[0].Size() << " (X); " 
         << XYZ[1].Size() << " (Y); " << XYZ[2].Size() << " (Z); "
         << Dist.Size() << " (Dist).\n"; 
    fflush( stdout );
  }
  for( xct = 0; xct < XYZ[0].Size(); xct++ ){ // loop over x parameter
    for( yct = 0; yct < XYZ[1].Size(); yct++ ){ // loop over y param.
      for( zct = 0; zct < XYZ[2].Size(); zct++ ){ // loop over z param
        for( dct = 0; dct < Dist.Size(); dct++ ){ // loop over distance
          if( Map[xct][yct][zct][dct] == -1.0 ){
            // If Map != -1.0, then this value has been read
            // from the input file and we don't need to recalculate
            
            // Set the parameters to their required values:
            SetParam( psr, XYZ[0], XX[xct] );
            SetParam( psr, XYZ[1], YY[yct] );
            SetParam( psr, XYZ[2], ZZ[zct] );
            rv = SetDist( psr, DD[dct] ); 
	    if (rv) {
              Map[xct][yct][zct][dct] = 1e6;
            } else {

            
            if( IP.DoOm() == 1 ) SetOmDot( psr );
            
            // Now refit:
            for( ct = 0; ct < IP.fitit() + 1; ct++ ){
	      updateBatsAll( psr, *npsr );
              formResiduals( psr, *npsr, 1 );
              if( ct < IP.fitit() )
                doFitAll( psr, *npsr, NULL );
              else
                calcRMS( psr, 0 );
            }
            if( fmod( IP.Realtime(), 3.0 ) == 0.0 ){
              PlotResiduals( psr );
            }
            Map[xct][yct][zct][dct] = psr[0].fitChisq;
            if( psr[0].fitChisq > 1e5 ){
              if( verb > 10 )
                cout << "Chisquared " << psr[0].fitChisq << " > 100,000.\t"
                     << " Redoing fit.\n"; fflush( stdout );
              if( dct > 0 )
                dct--;
              else if( zct > 0 ){
                zct--;
                dct = Dist.Size() - 1;
              }else if( yct > 0 ){
                yct--;
                zct = XYZ[2].Size() - 1;
                dct = Dist.Size() - 1;
              }else{
                cout << "Initial fit failed. Ignoring and moving on.\n";
                fflush( stdout );
              }
              // Reset timing model
              for( ct = 0; ct < parID.size(); ct++ )
                psr[0].param[parID[ct]].val[derID[ct]] = parvals[ct];
            }
	    } // GD
            
            counter = dct + zct * Dist.Size() + yct * XYZ[2].Size() * Dist.Size()
              + xct * XYZ[1].Size() * XYZ[2].Size() * Dist.Size();
            
            if( fmod( calculated, 100.0 ) == 0 ){
              if( verb > -1 )
                Progress( start, counter, size );
              // Plot progress
              if( fmod( IP.Realtime(), 2.0 ) == 0.0 )
                MakePlots( XYZ, Map, IP, Dist );
            }
            
            // Step 3b: Continuously write results to file for future resumption.
            if( verb > 10 ){
              cout << "Cube XX: " << XX[xct] << " YY: " << YY[yct]
                   << " ZZ: " << ZZ[zct] << " DD: " << DD[dct]
                   << " chisq: " << Map[xct][yct][zct][dct] << endl;
              fflush( stdout );
            }
            WriteResults( counter, Map[xct][yct][zct][dct], XX[xct], YY[yct], 
                          ZZ[zct], DD[dct], parmap, IP.form(), fout );
            calculated++;
          }// End of "if not read from input" clause
        } // end of loop over distance parameter
      } // end loop over Z parameter
 
      if( fmod( calculated, 20.0 ) == 0 && calculated > 0 ){
        if( verb > -1 )
          Progress( start, counter, size );
        // Plot progress
        if( fmod( IP.Realtime(), 2.0 ) == 0.0 )
          MakePlots( XYZ, Map, IP, Dist );
      }
    } // end loop over Y parameter

    if( fmod( calculated, 20.0 ) == 0.0 && calculated > 0 ){
      if( verb > -1 )
        Progress( start, counter, size );
      // Plot progress
      if( fmod( IP.Realtime(), 2.0 ) == 0.0 )
        MakePlots( XYZ, Map, IP, Dist );
    }
  }// end loop over X parameter

  fclose( fout );
  cout << endl;
} // end of CalcCube function

void FillArrays( long double *XX, long double *YY, long double *ZZ, 
                 long double *DD, Param *XYZ, Param Dist ){
  // This function fills the x, y and z values for all pixels of the
  // map to be created.

  int ct1, ct2, ct3; // generic counters

  for( ct1 = 0; ct1 < XYZ[0].Size(); ct1++ ){
    XX[ct1] = XYZ[0].I() + (long double) ct1 * XYZ[0].step();
    if( fabs( XX[ct1] ) < 1.0e-18 )
      XX[ct1] = 0.0;
  }

  for( ct2 = 0; ct2 < XYZ[1].Size(); ct2++ ){
    YY[ct2] = XYZ[1].I() + (long double) ct2 * XYZ[1].step();
    if( fabs( YY[ct2] ) < 1.0e-18 )
      YY[ct2] = 0.0;
  }

  for( ct3 = 0; ct3 < XYZ[2].Size(); ct3++ ){
    ZZ[ct3] = XYZ[2].I() + (long double) ct3 * XYZ[2].step();
    if( fabs( ZZ[ct3] ) < 1.0e-18 )
      ZZ[ct3] = 0.0;
  }

  for( ct1 = 0; ct1 < Dist.Size(); ct1++ ){
    DD[ct1] = Dist.I() + (long double) ct1 * Dist.step();
    if( fabs( DD[ct1] ) < 1.0e-18 )
      DD[ct1] = 0.0;
  }

  if( verb > 1 )
    cout << "Finished filling the coordinate arrays!\n";
}

void Resume( Param *XYZ, Param Dist, InputParams *IP, FILE *fout, pulsar *psr, 
             long double *XX, long double *YY, long double *ZZ, long double *DD,
             long double ****Map, char *parFile, char *timFile,
             int parmap[] ){
  //FILE *fout;
  // First read the header:
  char *names[3]; // Parameter names
  for( int ct = 0; ct < 3; ct++ )
    names[ct] = (char *)malloc( 20 * sizeof( char ) );
  char psrname[20]; // pulsar name
  char par[100], tim[100]; // par and tim file

  int testval = ReadHeader( IP->form(), IP->fout(), names, psrname, par, tim );
  if( verb > 100 )
    cout << "In Resume. ReadHeader result: " << testval << endl; fflush( stdout );

  if( testval != IP->form() ){
    if( verb > 10 ){
      cout << "Had >>" << IP->form() << "<< and want >>" << testval << "<<.\n";
      cout << "\t\t\t\t changing...\n"; fflush( stdout );
    }
    IP->ReForm(); // Toggle format (ascii/binary)
  }

  // Now reopen the file and get rid of the header:
  char read[1000];
  if( IP->form() == 0 ){
    // binary
    fout = fopen( IP->fout(), "rb" );
    fread( read, 100, 1, fout );
  }else{
    fout = fopen( IP->fout(), "r" );
    fgets( read, 1000, fout );
    if( verb ) 
      cout << "Purged header >>" << read << "<<..\n"; fflush( stdout );
  }
  if( verb > 100 )
    cout << "File reopened and header purged.\n"; fflush( stdout );

  int changed = 0;

  // Verify/check header information
  if( strcmp( psrname, psr[0].name ) != 0 ){
    cout << "Different pulsar! >>" << psrname << "<< from input data file "
         << " != >>" << psr[0].name << "<< from input par file.\n";
    fflush( stdout );
    changed = 1;
  }

  if( changed == 0 ){
    if( strcmp( par, parFile ) != 0 ){
      cout << "\n"
           << "WARNING!! Same pulsar but different par file.\n"
           << "          We will continue shortly, assuming everything is OK.\n"
           << "          But really you should press Ctrl-c NOW\n"
           << "          and provide another file name.\n";
      GiveMeABreak();
    }
    if( strcmp( tim, timFile ) != 0 ){
      cout << "\n"
           << "WARNING!! Same pulsar but different tim file.\n"
           << "          We will continue shortly, assuming everything is OK.\n"
           << "          But really you should press Ctrl-c NOW\n"
           << "          and provide another file name.\n";
      GiveMeABreak();
    }
  }

  // Figure out the organisation of the parameters
  if( verb > 100 )
    cout << "Figuring out the organisation of the parameters.\n"; fflush( stdout );
  if( changed == 0 ){
    if( strcmp( names[0], XYZ[0].name() ) == 0 ){
      parmap[0] = 0;
      if( strcmp( names[1], XYZ[1].name() ) == 0 ){
        parmap[1] = 1;
        if( strcmp( names[2], XYZ[2].name() ) == 0 ){
          parmap[2] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else if( strcmp( names[1], XYZ[2].name() ) == 0 ){
        parmap[2] = 1;
        if( strcmp( names[2], XYZ[1].name() ) == 0 ){
          parmap[1] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else{
        // Parameter not used in previous file.
        // Change to new output file:
        IP->Rename();
        changed = 1;
      }
    }else if( strcmp( names[0], XYZ[1].name() ) == 0 ){
      parmap[1] = 0;
      if( strcmp( names[1], XYZ[0].name() ) == 0 ){
        parmap[0] = 1;
        if( strcmp( names[2], XYZ[2].name() ) == 0 ){
          parmap[2] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else if( strcmp( names[1], XYZ[2].name() ) == 0 ){
        parmap[2] = 1;
        if( strcmp( names[2], XYZ[0].name() ) == 0 ){
          parmap[0] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else{
        // Parameter not used in previous file.
        // Change to new output file:
        IP->Rename();
        changed = 1;
      }
    }else if( strcmp( names[0], XYZ[2].name() ) == 0 ){
      parmap[2] = 0;
      if( strcmp( names[1], XYZ[0].name() ) == 0 ){
        parmap[0] = 1;
        if( strcmp( names[2], XYZ[1].name() ) == 0 ){
          parmap[1] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else if( strcmp( names[1], XYZ[1].name() ) == 0 ){
        parmap[1] = 1;
        if( strcmp( names[2], XYZ[0].name() ) == 0 ){
          parmap[0] = 2;
        }else{
          // Parameter not used in previous file.
          // Change to new output file:
          IP->Rename();
          changed = 1;
        }
      }else{
        // Parameter not used in previous file.
        // Change to new output file:
        IP->Rename();
        changed = 1;
      }
    }else{
      // Parameter not used in previous file.
      // Change to new output file:
      IP->Rename();
      changed = 1;
    }
  }

  if( verb > 100 )
    cout << "Sorted out organisation. Now going to read data.\n"; fflush( stdout );
  if( changed == 0 ){
    int xct, yct, zct, dct;
    // The same parameters are being investigated
    // Read data and fill into the array.
    long int counter;
    long double Mapval, vals[3], dval;
    while( !feof( fout ) ){

      if( IP->form() == 0 ){
        // binary data
        fread( &counter, sizeof( long int ), 1, fout );
        fread( &Mapval, sizeof( long double ), 1, fout );
        fread( &dval, sizeof( long double ), 1, fout );
        fread( &vals[0], sizeof( long double ), 1, fout );
        fread( &vals[1], sizeof( long double ), 1, fout );
        fread( &vals[2], sizeof( long double ), 1, fout );
      }else{
        // ascii data
        fscanf( fout, "%ld %Lg %Lg %Lg %Lg %Lg\n", &counter, &Mapval, &dval,
                &vals[0], &vals[1], &vals[2] );
      }
      if( fabs( dval ) < 1.0e-18 )
        dval = 0.0;
      if( (long double)fabs( vals[0] ) < (long double)1.0e-18 )
        vals[0] = 0.0;
      if( fabs( vals[1] ) < (long double)1.0e-18 )
        vals[1] = 0.0;
      if( fabs( vals[2] ) < (long double)1.0e-18 )
        vals[2] = 0.0;
      
      if( verb > 100 ){
        cout << "Resuming... Read: ct:" << counter << "\tchisq: " 
             << Mapval << "\tdist: " << dval << " X: " << vals[0] << " Y: " 
             << vals[1] << " Z: " << vals[2] << endl;
        fflush( stdout );
      }
      
      // Find read data point
      xct = (int)(0.5+ ( vals[parmap[0]] - XYZ[0].I() ) / XYZ[0].step() );
      yct = (int)(0.5+ ( vals[parmap[1]] - XYZ[1].I() ) / XYZ[1].step() );
      zct = (int)(0.5+ ( vals[parmap[2]] - XYZ[2].I() ) / XYZ[2].step() );
      dct = (int)(0.5+ ( dval - Dist.I() ) / Dist.step() );

      //cout << "JORIS: COSI: " << vals[2]; 
      //cout << "\t index: " << zct << endl; fflush( stdout );

      if( verb > 7 ){
        cout << "Found indices: x:" << xct << " y: " << yct 
             << " z: " << zct << " dist: " << dct << endl;
        fflush( stdout );
      }

      if( xct >= 0 && yct >= 0 && zct >= 0 && dct >= 0
          && xct < XYZ[0].Size() && yct < XYZ[1].Size() && zct < XYZ[2].Size() 
          && dct < Dist.Size() ){
        if( (float)XX[xct] == (float)vals[parmap[0]] && 
            (float)YY[yct] == (float)vals[parmap[1]] &&
            (float)ZZ[zct] == (float)vals[parmap[2]] &&
            (float)DD[dct] == (float)dval ){
          Map[xct][yct][zct][dct] = Mapval;
          if( verb > 10 ){
            cout << "Resuming... entered " << Mapval << endl;
            fflush( stdout );
          }
        }else{
          // No exact correspondence in the required array. 
          // Ignore input values.
          if( verb > 3 ){
            cout << "Resuming.. Entry didn't match: (X):" 
                 << XX[xct] << " != " << vals[parmap[0]]
                 << "   (Y):" << YY[yct] << " != " << vals[parmap[1]] 
               << "   (Z):" << ZZ[zct] << " != " << vals[parmap[2]] 
                 << "\t(D): " << DD[dct] << " != " << dval << endl;
            fflush( stdout );
          }// End verbosity
        }// End else-clause
      }
    } // end of "feof" loop
    
  }else if( changed == 1 ){
    // Other parameters are being investigated.
    // In this case things are very simple: 
    // We just have to start all over again.
    if( IP->form() == 0 )
      fout = fopen( IP->fout(), "ab" );
    else
      fout = fopen( IP->fout(), "a" );
    PutHeader( fout, XYZ, *IP, psr, parFile, timFile );
  }

} // End of resume function.

void PutHeader( FILE *fout, Param *XYZ, InputParams IP, pulsar *psr,
                char *parFile, char *timFile ){
  // Write header information
  char header[1000];
  if( IP.form() == 0 ){
    // binary header
    sprintf( header, "%s %s %s %s %s %s %s\n", psr[0].name, parFile, timFile, 
             XYZ[0].name(), XYZ[1].name(), XYZ[2].name(), "binary" );
    fwrite( header, 100, 1, fout );
  }else{

    // ascii header
    fprintf( fout, "%s %s %s %s %s %s %s\n", psr[0].name, parFile, timFile,
             XYZ[0].name(), XYZ[1].name(), XYZ[2].name(), "ascii" );
  }
    
  fflush( fout );
}

// =================
// Generic functions
// =================

// Feedback on the successful (or not) reading in of command-line argument
void Feedback( const char param[], char *parval ){
  //cout << "\nJORIS in Feedback with : " << verb << endl;
  if( verb > 0 )
    cout << "Read " << param << ": " << parval << endl;
}

void Feedback( int nread, const char param[], int parval, const char label[] ){
  double temp = (double)parval;
  Feedback( nread, param, temp, label );
}

void Feedback( int nread, const char param[], long double parval, 
               const char label[] ){
  double temp = (double) parval;
  Feedback( nread, param, temp, label );
}

void Feedback( int nread, const char param[], double parval, 
               const char label[] ){
  if( nread != 1 )
    cout << "===============================================\n"
         << "ERROR!\n"
         << "      Failed to read \"" << param << "\" value.\n"
         << "      Will try to continue with " << param << " = " << parval
         << ".\n"
         << "Try using [" << label << " value].\n"
         << "===============================================\n";
  if( verb > 0 )
    cout << "Read " << param << ": " << parval << endl;
}

// Help on command-line arguments of the plug-in:
void help(){
  cout << "ChiSqCube -- a Tempo2 plug-in for chi-squared mapping.\n"
       << "\n"
       << "Available command-line options:\n"
       << " -grdev 11/xs   : Determine graphics device.\n"
       << " -all           : Prints all known parameters.\n"
       << " -out filename  : Determine the output file name.\n"
       << " -v/-verb x     : sets the verbosity to level <x> \n"
       << "                  (positive integer).\n"
       << " -x KOM 0 1 360 : Determines the first parameter (KOM), to be \n"
       << "                  sampled from 0 to 360, at intervals of 1.\n"
       << " -y/-z          : Works like -x, but for the second and \n"
       << "                  third parameters.\n"
       << " -param         : Works like -x but for any arbitrary parameter.\n"
       << " -params KOM 0 1 360 KIN 0 1 180 M2 0.1 0.3 0.01 \n"
       << "                : Allows determination of all parameters at once.\n"
       << " -dist 0.1 0.01 2 : Distance range (in kpc): min, stepsize, max.\n"
       << " -plotonly      : Plot what's in the file, don't calculate \n"
       << "                  anything. (Uncalculated pixels will be defaulted\n"
       << "                  to the maximal value.)\n"
       << " -ascii         : Output ascii (instead of binary).\n"
       << " -PlotRealtime  : Update the maps in realtime (once every 20 "
       << " pixels).\n"
       << " -ResidualPlot  : Plot the postfit residuals for every pixel.\n"
       << " -cont          : Determines the contour levels to be plotted.\n"
       << "                  (Followed by the contour level values.)\n"
       << " -UseOmdot      : Hardwire Omdot to its GR prediction.\n"
       << " -nits x        : Number of fitting iterations.\n"
       << " -transform     : Transform the (output) file from ascii to binary"
       << " or vice versa.\n"
       << " -merge 1 2 3   : Merge files 1, 2 and 3 into filename (defined by "
       << " -out.\n"
       << endl;
}

void palett(){
  float GL[] = { 0.0, 1.0 };
  float GR[] = { 0.0, 1.0 };
  float GG[] = { 0.0, 1.0 };
  float GB[] = { 0.0, 1.0 };
  cpgctab( GL, GR, GG, GB, 2, 1.0, 0.5 );
}

void PlotResiduals( pulsar *psr ){
  int ct; // generic counter

  // Prepare data arrays:
  float *XX, *YY, *EE;
  XX = (float *)malloc( psr[0].nobs * sizeof( float ) );
  YY = (float *)malloc( psr[0].nobs * sizeof( float ) );
  EE = (float *)malloc( psr[0].nobs * sizeof( float ) );

  // Enter data
  for( ct = 0; ct < psr[0].nobs; ct++ ){
    XX[ct] = (float)psr[0].obsn[ct].bat;
    YY[ct] = (float)psr[0].obsn[ct].residual;
    EE[ct] = (float)psr[0].obsn[ct].toaErr*1e-6;
  }

  // Find extremes
  float xmin, xmax, ymin, ymax;
  xmin = XX[0];
  xmax = XX[0];
  ymin = YY[0] - EE[0];
  ymax = YY[0] + EE[0];
  for( ct = 1; ct < psr[0].nobs; ct++ ){

    if( XX[ct] < xmin ) 
      xmin = XX[ct];

    if( XX[ct] > xmax )
      xmax = XX[ct];
    
    if( YY[ct] - EE[ct] < ymin )
      ymin = YY[ct] - EE[ct];
    
    if( YY[ct] + EE[ct] > ymax )
      ymax = YY[ct] + EE[ct];
  }

  cpgbeg( 0, "1982/xs", 1, 1 );
  cpgscir( 0, 15 );
  cpgenv( xmin, xmax, ymin, ymax, 0, 2 );
  cpgpt( psr[0].nobs, XX, YY, 1 );
  cpgerrb( 6, psr[0].nobs, XX, YY, EE, 0.5 );
  cpgend( );

  free( XX ), free( YY ), free( EE );
}

int SetDist( pulsar *psr, long double dd ){
  // First Parallax:
  if( psr[0].param[param_px].paramSet[0] == 1 ){
    psr[0].param[param_px].val[0] = 1.0/dd;
    psr[0].param[param_px].fitFlag[0] = 0;
  }else{
    // Have not defined PX in the timing model.
    if( verb >= 0 ){
      cout << "WARNING!\n"
           << "       You have not defined parallax in your par file.\n"
           << "       It really would make sense to do this.\n"
           << "       Just so you know. Please reconsider.\n";
    }
  }

  // Now determine companion mass and inclination angle:

  long double m2=-1.0, sini = -2.0;
  if( strcmp( psr[0].binaryModel, "DDH" ) == 0 ){
    // Using DDH model (H3, stig)
    m2 = psr[0].param[param_h3].val[0] / TSUN 
      / pow( psr[0].param[param_stig].val[0], 3.0 );
    sini = 2.0 * psr[0].param[param_stig].val[0]
      / ( 1.0 + pow( psr[0].param[param_stig].val[0], 2.0 ) );
  }else if( strcmp( psr[0].binaryModel, "ELL1H" ) == 0 ){
    // Using ELL1H model (H3, H4)
    m2 = pow( psr[0].param[param_h3].val[0], 4.0 ) / TSUN
      / pow( psr[0].param[param_h4].val[0], 3.0 );
    sini = 2.0 * psr[0].param[param_h3].val[0] * psr[0].param[param_h4].val[0]
      / ( pow( psr[0].param[param_h3].val[0], 2.0 )
          + pow( psr[0].param[param_h4].val[0], 2.0 ) );
  }else{
    // Should have a model with simple Shapiro delay parameters.
  }

  if( psr[0].param[param_m2].paramSet[0] == 1 )
    m2 = psr[0].param[param_m2].val[0];
  if( psr[0].param[param_sini].paramSet[0] == 1 )
    sini = psr[0].param[param_sini].val[0];

  if( m2 == -1.0 || sini == -2.0 ){
    // Have timing model with no m2 or sini; and no DDH or ELL1H model.
    // Cannot use the distance in certain cases.
    if( verb >= 0 ){
      cout << "WARNING!\n"
           << "        You have no Shapiro Delay parameters defined.\n"
           << "        This means we cannot calculate the GR contribution to Pbdot.\n"
           << "        We will ignore this contribution, but make sure to verify\n"
           << "        that this is justified, before using any results from this plug-in.\n";
    }
  }

  // Now Pbdot:
  //
  // Kinematic (proper-motion related):
  long double Kinematic;
  if( psr[0].param[param_pmra].paramSet[0] * 
      psr[0].param[param_pmdec].paramSet[0] == 1 ){
    Kinematic = dd * 1000.0 * PARSEC // distance in meters
      / 299792458 // speed of light in meters per second
      * psr[0].param[param_pb].val[0] * 24.0 * 3600.0 // Pb in seconds
      * ( pow( psr[0].param[param_pmra].val[0] * 1e-3 / 3600.0
               * M_PI / 180.0 / 365.25 / 24.0 / 3600.0, 2.0 ) + 
          pow( psr[0].param[param_pmdec].val[0] * 1e-3 / 3600.0
               * M_PI / 180.0 / 365.25 / 24.0 / 3600.0, 2.0 ) ); // PM in rad/s
  }else{
    // No proper motion
    Kinematic = 0.0L;
  }

  long double Gal[2]; // Galactic coordinates (0: lat; 1: long)
  GetGalCoord( psr[0].name, Gal );

  // Kz (vertical Galactic potential):
  long double Kz;
  Kz = - sin( Gal[0] ) * psr[0].param[param_pb].val[0] * 24.0 * 3600.0 
    / 299792458 // speed of light in meters per second
    * GetKz( dd * 1000.0 * sin( Gal[0] ) );

  // DGR (differential Galactic rotation):
  long double DGR;
  long double beta = dd/8.34 * cos( Gal[0] ) - cos( Gal[1] );
  DGR = -pow( 240e3, 2.0 ) // Galactic rotation in m/s
    /299792458 // speed of light in m/s
    /8340. / PARSEC // Galactic radius in m
    * cos( Gal[0] ) * ( cos( Gal[1] ) + 
        beta / ( pow( sin( Gal[1] ),2.0 ) + pow( beta, 2.0 ) ) ) *
    psr[0].param[param_pb].val[0] * 24.0 * 3600.0;

  // Relativistic effect (caused by GW emission):
  long double relativistic = 0.0L;
  long double Mpsr = 0.0L;

  if( m2 != -1.0 && sini != -2.0 ){
    Mpsr = -m2 +
      sqrt( 4.925490947e-6 * 
            pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
                 24.0 * 3600.0, 2.0 ) *
            pow( m2 * sini / 
                 psr[0].param[param_a1].val[0], 3.0 ) );
    if( verb > 100 ) 
      cout << "Have Mpsr: " << Mpsr << endl; fflush( stdout );
    if( Mpsr < 0.0 || isnan(Mpsr)){
      cout << "ERROR! Used:\n"
           << " m2: " << m2 << endl
           << " pb: " << psr[0].param[param_pb].val[0] << endl
           << " a1: " << psr[0].param[param_a1].val[0] << endl
           << " sini: " << sini << endl
           << "Resulting in : " <<  -m2 
           << " + " 
           << " sqrt( " << 
            pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
                 24.0 * 3600.0, 2.0 ) 
           << " * " << 
        pow( m2 * sini / 
             psr[0].param[param_a1].val[0], 3.0 ) 
           << " = " << 
        sqrt( 4.925490947e-6 * 
              pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
                   24.0 * 3600.0, 2.0 ) *
              pow( m2 * sini / 
                   psr[0].param[param_a1].val[0], 3.0 ) )
           << endl;
      return (1);
      //exit( 0 );
    }

    relativistic = -192.0*M_PI/5.0 * 
      pow( 2.0 * M_PI * 4.925490947e-6, 5.0/3.0 ) * 
      pow( psr[0].param[param_pb].val[0] * 24.0 * 3600.0, -8.0/3.0 ) * 
      Mpsr * m2 * 
      pow( Mpsr + m2, -1.0/3.0 ) * 
      ( 1.0 + 73.0/24.0 * pow( psr[0].param[param_ecc].val[0], 2.0 ) +
        37.0/96.0 * pow( psr[0].param[param_ecc].val[0], 4.0 ) ) *
      pow( 1.0 - psr[0].param[param_ecc].val[0], -7.0/2.0 ) *
      psr[0].param[param_pb].val[0] * 24.0 * 3600.0;
  }

  if( verb > 100 ) 
    cout << "GD Contributions to Pbdot: " 
         << "Kinematic: " << setprecision( 10 ) << Kinematic
         << "\tKz: " << setprecision( 10 ) << Kz 
         << "\tDGR: " << setprecision( 10 ) << DGR
         << "\tRelativistic: " << setprecision( 10 ) << relativistic
         << endl; fflush( stdout );
  if( verb > 9 )
    cout << "Setting Pbdot to: " << setprecision( 10 )
         << Kinematic + Kz + DGR + relativistic << endl; fflush( stdout );

  
  if( psr[0].param[param_pbdot].paramSet[0] != 1 ){
    // initialise pbdot if needed. (NOT CHECKED JORIS );
    psr[0].param[param_pbdot].paramSet[0] = 1;
    strcpy( psr[0].param[param_pbdot].label[0], "PBDOT" );
    strcpy( psr[0].param[param_pbdot].shortlabel[0], "PBDOT" );
    psr[0].param[param_pbdot].fitFlag[0] = 0;
    psr[0].param[param_pbdot].err[0] = 0.0;
  }
  psr[0].param[param_pbdot].val[0] = Kinematic + Kz + DGR + relativistic;
  psr[0].param[param_pbdot].fitFlag[0] = 0;
  return 0;
}

void SetOmDot( pulsar *psr ){
  long double Mtot;
  long double OrigOmdot = psr[0].param[param_omdot].val[0];

  if( verb > 100 ){
    cout << "In SetOmDot with kin: " << psr[0].param[param_kin].val[0]
         << " and sini: " << psr[0].param[param_sini].val[0]
         << "  m2: " << psr[0].param[param_m2].val[0] << endl;
    fflush( stdout );
  }

  if( psr[0].param[param_kin].paramSet[0] == 0 ){
    Mtot = psr[0].param[param_pb].val[0] * DAY2S / ( 2.0 * M_PI ) 
      * sqrt( 4.925490947e-6 ) 
      * pow( psr[0].param[param_m2].val[0] 
             * fabs( psr[0].param[param_sini].val[0] )
             / psr[0].param[param_a1].val[0], 1.5 );
  }else{
    Mtot = psr[0].param[param_pb].val[0] * DAY2S / ( 2.0 * M_PI )
      * sqrt( 4.925490947e-6 )
      * pow( psr[0].param[param_m2].val[0]
             * fabs( sin( psr[0].param[param_kin].val[0] / 180.0 * M_PI ) )
             / psr[0].param[param_a1].val[0], 1.5 );
    if( verb > 100 )
      cout << "Setting Omdot. Mtot = " << Mtot 
           << "\t m2: " << psr[0].param[param_m2].val[0] << endl; fflush( stdout );
  }

  psr[0].param[param_omdot].val[0] = 0.19738 * pow( Mtot, 2.0 / 3.0 ) 
    * pow( psr[0].param[param_pb].val[0], - 5.0 / 3.0 ) 
    / ( 1 - pow( psr[0].param[param_ecc].val[0], 2.0 ) );
  psr[0].param[param_omdot].fitFlag[0] = 0;
  
  psr[0].param[param_pb].val[0] = 
    psr[0].param[param_pb].val[0] * 
    ( 1.0 - ( OrigOmdot - psr[0].param[param_omdot].val[0] ) / 360.0 
      / 365.25 * psr[0].param[param_pb].val[0] );

  if( verb > 9 )
    cout << "TOALAL mass: " << Mtot 
         << " OMDOT: " << psr[0].param[param_omdot].val[0] 
         << " PB: " << setprecision( 20 ) << psr[0].param[param_pb].val[0]
         << endl; fflush( stdout );
}  

void Progress( time_t start, int counter, double size ){
  time_t present;
  double hours, minutes, seconds;
  double PctDone, ToDo;
  double diff;

  // Print progress
  time( &present );
  diff = difftime( present, start );
  PctDone = (double)counter * 100.0 / size;
  ToDo = diff * ( 100.0 - PctDone ) / PctDone;
  seconds = modf( ToDo / 60.0, &minutes ) * 60.0;
  minutes = modf( minutes / 60.0, &hours ) * 60.0;
  cout << "Finished iteration " << counter 
       << " out of " << size
       << " or " << PctDone
       << "\% in " << diff << " seconds. (" 
       << hours << "h:" << minutes << "m " 
       << " left).       \r";
  fflush( stdout );
}

void WriteResults( long int counter, long double val, long double xx,
                   long double yy, long double zz, long double dd, int *map,
                   int form, FILE *fout ){
  if( form == 0 ){
    // binary output
    fwrite( &counter, sizeof( long int ), 1, fout );
    fwrite( &val, sizeof( long double ), 1, fout );
    fwrite( &dd, sizeof( long double ), 1, fout );

    // Write first parameter
    if( map[0] == 0 )
      fwrite( &xx, sizeof( long double ), 1, fout );
    else if( map[1] == 0 )
      fwrite( &yy, sizeof( long double ), 1, fout );
    else if( map[2] == 0 )
      fwrite( &zz, sizeof( long double ), 1, fout );

    // Write second parameter
    if( map[0] == 1 )
      fwrite( &xx, sizeof( long double ), 1, fout );
    else if( map[1] == 1 )
      fwrite( &yy, sizeof( long double ), 1, fout );
    else if( map[2] == 1 )
      fwrite( &zz, sizeof( long double ), 1, fout );

    // Write third parameter
    if( map[0] == 2 )
      fwrite( &xx, sizeof( long double ), 1, fout );
    else if( map[1] == 2 )
      fwrite( &yy, sizeof( long double ), 1, fout );
    else if( map[2] == 2 )
      fwrite( &zz, sizeof( long double ), 1, fout );
    
    fflush( fout );
  }else{
    // ascii output
    fprintf( fout, "%ld %Lg %Lg ", counter, val, dd );

    // Write first parameter
    if( map[0] == 0 )
      fprintf( fout, "%Lg ", xx );
    else if( map[1] == 0 )
      fprintf( fout, "%Lg ", yy );
    else if( map[2] == 0 )
      fprintf( fout, "%Lg ", zz );

    // Write second parameter
    if( map[0] == 1 )
      fprintf( fout, "%Lg ", xx );
    else if( map[1] == 1 )
      fprintf( fout, "%Lg ", yy );
    else if( map[2] == 1 )
      fprintf( fout, "%Lg ", zz );

    // Write third parameter
    if( map[0] == 2 )
      fprintf( fout, "%Lg\n", xx );
    else if( map[1] == 2 )
      fprintf( fout, "%Lg\n", yy );
    else if( map[2] == 2 )
      fprintf( fout, "%Lg\n", zz );
    fflush( fout );
  }
}

void Transform( InputParams *IP ){
  FILE *fin;
  char read[1000];

  // Now read header and verify format
  char *names[3]; // Parameter names
  for( int ct = 0; ct < 3; ct++ )
    names[ct] = (char *)malloc( 50* sizeof( char ) );
  char psrname[20]; // pulsar name
  char par[100], tim[100]; // par and tim file
  char format[35]; // binary or ascii

  int testval = ReadHeader( IP->form(), IP->fout(), names, psrname, par, tim );
  if( testval != IP->form() )
    IP->ReForm();

  // Now reopen the file and get rid of the header
  if( IP->form() == 0 ){
    // binary
    fin = fopen( IP->fout(), "rb" );
    fread( read, 100, 1, fin );
  }else{
    fin = fopen( IP->fout(), "r" );
    fgets( read, 1000, fin );
  }

  // Now read file and write out in different format
  FILE *fout; // output file
  char outname[100]; 
  if( IP->form() == 1 ){
    // binary output file
    sprintf( outname, "%s.bin", IP->fout() );
    fout = fopen( outname, "wb" );
    char header[1000];
    sprintf( header, "%s %s %s %s %s %s %s\n", psrname, par, tim, 
             names[0], names[1], names[2], "binary" );
    fwrite( header, 100, 1, fout );
  }else{
    // ascii output file
    sprintf( outname, "%s.asc", IP->fout() );
    if( verb > 10 )
      cout << "Opening output file: " << outname << endl; fflush( stdout );
    fout = fopen( outname, "w" );
    fprintf( fout, "%s %s %s %s %s %s %s\n", psrname, par, tim, names[0],
             names[1], names[2], "ascii" );
  }

  long int counter;
  long double Mapval, vals[3];
  while( !feof( fin ) ){
    if( IP->form() == 1 ){
      // ascii input, binary output
      fscanf( fin, "%ld %Lg %Lg %Lg %Lg", &counter, &Mapval, &vals[0], &vals[1],
              &vals[2] );
      // now write:
      fwrite( &counter, sizeof( long int ), 1, fout );
      fwrite( &Mapval, sizeof( long double ), 1, fout );
      fwrite( &vals[0], sizeof( long double ), 1, fout );
      fwrite( &vals[1], sizeof( long double ), 1, fout );
      fwrite( &vals[2], sizeof( long double ), 1, fout );
      fflush( fout );
    }else{
      // binary input, ascii output
      fread( &counter, sizeof( long int ), 1, fin );
      fread( &Mapval, sizeof( long double ), 1, fin );
      fread( &vals[0], sizeof( long double ), 1, fin );
      fread( &vals[1], sizeof( long double ), 1, fin );
      fread( &vals[2], sizeof( long double ), 1, fin );
      // now write:
      fprintf( fout, "%ld %Lg %Lg %Lg %Lg\n", counter, Mapval, vals[0], vals[1],
               vals[2] );
      fflush( fout );
    }
  }
  fclose( fin );
  fclose( fout );
}

void Merge( InputParams *IP, int argc, char *argv[] ){
  char **files; // Files to merge
  int ct; // generic counter
  char tmpname[100];
  int nread;
  
  int Nfiles = 0;
  for( ct = 1; ct < argc; ct++ ){
    if( strcmp( argv[ct], "-merge" ) == 0 ){
      int testID = ct;
      while( ct < argc - 1 && argv[ct+1][0] != '-' ){
        Nfiles++;
        ct++;
      }
      ct = testID;
      files = (char **)malloc( Nfiles * sizeof( char *) );
      for( int ct2 = 0; ct2 < Nfiles; ct2++ )
        files[ct2] = (char *)malloc( 100 * sizeof( char ) );

      while( ct < argc - 1 && argv[ct+1][0] != '-' ){
        nread = sscanf( argv[++ct], "%s", tmpname );
        Feedback( "file to merge", tmpname );
        if( nread > 0 )
          strcpy( files[ct-testID-1], tmpname );
        //files.push_back( tmpname );
      }// Finished reading all files to be merged.
    }// End "if -merge"
  }// End loop over all arguments

  if( verb > 5 ){
    cout << "Will merge >>" << Nfiles << "<< files.\n"; 
    fflush( stdout );
    for( ct = 0; ct < Nfiles; ct++ )
      cout << "File " << ct << " :: " << files[ct] << endl; fflush( stdout );
  }
  
  // First read headers of all files 
  char ***names; // parameter names
  char **psr; // pulsar names
  char **par; // par-file names
  char **tim; // tim-file names
  // Three dimensions. First: number of input files:
  names = (char ***)malloc( Nfiles * sizeof( char ** ) );
  psr = (char **)malloc( Nfiles * sizeof( char * ) );
  par = (char **)malloc( Nfiles * sizeof( char * ) );
  tim = (char **)malloc( Nfiles * sizeof( char * ) );
  for( ct = 0; ct < Nfiles; ct++ ){
    // Second dimension: 3 parameter names for each file:
    names[ct] = (char **)malloc( 3 * sizeof( char * ) );
    psr[ct] = (char *)malloc( 20 * sizeof( char ) );
    par[ct] = (char *)malloc( 100 * sizeof( char ) );
    tim[ct] = (char *)malloc( 100 * sizeof( char ) );
    for( int ct2 = 0; ct2 < 3; ct2++ ){
      // Third dimension: 20 characters for each parameter name:
      names[ct][ct2] = (char *)malloc( 20 * sizeof( char ) );
    }
  }
  int formvals[Nfiles];
  for( ct = 0; ct < Nfiles; ct++ )
    formvals[ct] = ReadHeader( IP->form(), files[ct], names[ct], psr[ct], par[ct],
                              tim[ct] );

  // Now compare header information
  // Compare essentials: parameter names:
  int parmtest = 0;
  for( ct = 1; ct < Nfiles; ct++ )
    for( int ct2 = 0; ct2 < 3; ct2++ )
      parmtest += fabs( strcmp( names[0][ct2], names[ct][ct2] ) );

  if( parmtest != 0 ){
    cout << "ERROR! You are trying to merge files for different parameters.\n";
    cout << "       Maybe you just need to reorder the columns. Try converting\n";
    cout << "       to ascii and resorting the columns with awk.\n";
    cout << "EXITING...\n";
    fflush( stdout );
    exit( 1 );
  }

  // Compare non-essentials: pulsar, par and tim files:
  int psrtest = 0;
  int partest = 0;
  int timtest = 0;
  for( ct = 1; ct < Nfiles; ct++ ){
    psrtest += fabs( strcmp( psr[0], psr[ct] ) );
    partest += fabs( strcmp( par[0], par[ct] ) );
    timtest += fabs( strcmp( tim[0], tim[ct] ) );
  }
  if( psrtest != 0 ){
    cout << "WARNING! You seem to be merging files for different pulsars.\n"
         << "         Normally this is a problem, but maybe it's just a typo.\n"
         << "         We will proceed, assuming all well, but maybe you'd like\n"
         << "         to check your files.\n";
    GiveMeABreak();
  }

  if( partest != 0 ){
    cout << "WARNING! You seem to be merging files from different par files.\n"
         << "         This may be a problem, but maybe you're lucky.\n"
         << "         We will proceed, assuming all well, but maybe you'd like\n"
         << "         to check your files.\n";
    GiveMeABreak();
  }

  if( timtest != 0 ){
    cout << "WARNING! You seem to be merging files from different tim files.\n"
         << "         This may be a problem, but maybe you're lucky.\n"
         << "         We will proceed, assuming all well, but maybe you'd like\n"
         << "         to check your files.\n";
    GiveMeABreak();
  }

  // If we made it this far, all is well and we can start reading and writing.
  // First write header into new output file.
  
  // Open output file
  // Test whether file already exists or not
  FILE *fout = fopen( IP->fout(), "r" );
  if( fout == NULL ){
    // File does not yet exist. All is well.
  }else{
    // The file does already exist. We need to rename.
    fclose( fout );
    IP->Rename();
  }
  if( IP->form() == 0 ){
    fout = fopen( IP->fout(), "wb" );
    char header[1000];
    sprintf( header, "%s %s %s %s %s %s %s\n", psr[0], par[0], tim[0], names[0][0],
             names[0][1], names[0][2], "binary" );
    fwrite( header, 100, 1, fout );
  }else{
    fout = fopen( IP->fout(), "w" );
    fprintf( fout, "%s %s %s %s %s %s %s\n", psr[0], par[0], tim[0], names[0][0],
             names[0][1], names[0][2], "ascii" );
  }
  fflush( fout );
  if( fout != NULL && verb > 0 )
    cout << "Opened output file (2).\n"; fflush( stdout );

  // Now read and write data:
  FILE *fin;
  char read[1000];
  long int counter;
  long double Mapval, vals[3];
  for( ct = 0; ct < Nfiles; ct++ ){ // Foreach file
    // first open file and get rid of header:
    if( formvals[ct] == 0 ){
      // binary input file
      fin = fopen( files[ct], "rb" );
      fread( read, 100, 1, fin );
    }else{
      // ascii file
      fin = fopen( files[ct], "r" );
      fgets( read, 1000, fin );
    }
    while( !feof( fin ) ){
      // Read data
      if( formvals[ct] == 0 ){
        // binary input file
        fread( &counter, sizeof( long int ), 1, fin );
        fread( &Mapval, sizeof( long double ), 1, fin );
        fread( &vals[0], sizeof( long double ), 1, fin );
        fread( &vals[1], sizeof( long double ), 1, fin );
        fread( &vals[2], sizeof( long double ), 1, fin );
      }else{
        fscanf( fin, "%ld %Lg %Lg %Lg %Lg", &counter, &Mapval, &vals[0], &vals[1],
                &vals[2] );
      }

      // Write data
      if( IP->form() == 0 ){
        // binary output format
        fwrite( &counter, sizeof( long int ), 1, fout );
        fwrite( &Mapval, sizeof( long double ), 1, fout );
        fwrite( &vals[0], sizeof( long double ), 1, fout );
        fwrite( &vals[1], sizeof( long double ), 1, fout );
        fwrite( &vals[2], sizeof( long double ), 1, fout );
        fflush( fout );
      }else{
        // ascii output format
        fprintf( fout, "%ld %Lg %Lg %Lg %Lg\n", counter, Mapval, vals[0], vals[1],
                 vals[2] );
        fflush( fout );
      }
    }
    fclose( fin );
  }// end "foreach file"
  fclose( fout );
}

int ReadHeader( int form, char *fname, char *names[], char *psr, char *par,
                 char *tim ){
  // This function reads the header information from file "fname" and
  // puts the header information into the containers passed as
  // arguments. The file format of the input files is given as return
  // value.
  int returnval = -100;

  FILE *fin;
  char read[1000];

  if( form == 0 ){
    // binary
    fin = fopen( fname, "rb" );
    fread( read, 100, 1, fin );
  }else if( form == 1 ){
    // ascii
    fin = fopen( fname, "r" );
    fgets( read, 1000, fin );
  }else{
    cout << "ERROR! Don't know input format >>" << form << "<<.\n";
    cout << "       Will assume whatever is in input file.\n";
  }

  if( verb > 3 ) 
    cout << "Read input file header: >>" << read << "<<\n"; fflush( stdout );

  // Now read header information and verify format
  char format[100];
  sscanf( read, "%s %s %s %s %s %s %s", psr, par, tim, names[0], names[1], names[2],
          format );

  // Check format:
  if( strcmp( format, "ascii" ) == 0 && form != 1 ){
    cout << "Note you indicated the input format is binary, while it really "
         << "is ascii.\n"; fflush( stdout );
    returnval = 1;
  }else if( strcmp( format, "binary" ) == 0 && form != 0 ){
    cout << "Note you indicated the input format is ascii, while it really "
         << "is binary.\n"; fflush( stdout );
    returnval = 0;
  }else if( strcmp( format, "ascii" ) == 0 && form == 1 ){
    // ascii read and ascii expected. All is well.
    returnval = 1;
  }else if( strcmp( format, "binary" ) == 0 && form == 0 ){
    // binary expected and binary read. All is well.
    returnval = 0;
  }else{
    cout << "Don't know what's going on! Expected format >>" << form 
         << "<< and read format >>" << format << "<<. Whatever that means."
         << "\n     Panicking.... Quitting.\n";
    fflush( stdout );
    exit( 1 );
  }
  if( verb > 100 ) 
    cout << "ReadHeader finished.\n"; fflush( stdout );
  fclose( fin );
  return( returnval );
}

void GiveMeABreak(){
  for( int ct = 5; ct > 0; ct-- ){
    cout << "Continuing in " << ct << " sec...\r"; fflush( stdout );
    system( "sleep 1" );
  }
  cout << "Continuing...                     \n"; fflush( stdout );
}

/*
  The following function is copied from the LKB code provided with
  Verbiest, Weisberg, Chael, Lee & Lorimer, ApJ 2012.
 */
void GetGalCoord( char name[], long double *Gal ){
  // Determines the Galactic coordinates of the pulsar
  int lengte = strlen( name );

  int xtra; 
  if( name[0] == 'J' || name[0] == 'B' )
    xtra = 1;
  else
    xtra = 0;

  long double hours = 0.0;
  long double temp = 0.0;
  long double RA = 0.0;

  // Determine RA
  temp = (long double)( name[xtra] - '0' );
  hours += 10.0 * temp;
  temp = (long double)( name[++xtra] - '0' );
  hours += temp;

  temp = (long double)( name[++xtra] - '0' );
  hours += temp / 6.0;
  temp = (long double)( name[++xtra] - '0' );
  hours += temp / 60.0;
  RA = hours / 12.0 * M_PI;
  if( verb > 9 )
    cout << "RA: " << RA << " hours: " << hours << endl;

  // Declination
  int PosNeg = 1;
  if( name[++xtra] == '-' )
    PosNeg = -1;
  else if( name[xtra] == '+' )
    PosNeg = 1;
  else 
    printf( "ERROR in Dec sign?!: >>%c<<\n", name[xtra] );

  temp = 0.0;
  long double degrees = 0.0;
  long double Dec = 0.0; 

  temp = (long double)( name[++xtra] - '0' );
  degrees += 10.0 * temp;
  temp = (long double)( name[++xtra] - '0' );
  degrees += temp;

  if( lengte > 8 ){
    temp = (long double)( name[++xtra] - '0' );
    degrees += temp / 6.0;
    temp = (long double)( name[++xtra] - '0' );
    degrees += temp / 60.0;
  }

  degrees *= (long double) PosNeg;

  Dec = degrees / 180.0 * M_PI;
  if( verb > 9 ) 
    cout << "Dec: " << Dec << " degrees: " << degrees << endl;

  // Determine Galactic Coordinates
  // Latitude
  Gal[0] = asin( cos( Dec ) * cos( 27.4 / 180.0 * M_PI ) *
                 cos( RA - 192.25 / 180.0 * M_PI ) +
                 sin( Dec ) * sin( 27.4 / 180.0 * M_PI ));
  long double SinGL, CosGL;
  SinGL = ( sin( Dec ) * sin( 62.6 / 180.0 * M_PI ) +
            cos( Dec ) * sin( RA - 282.25 / 180.0 * M_PI ) *
            cos( 62.6 / 180.0 * M_PI )) / cos( Gal[0] );
  CosGL = ( cos( Dec ) * cos( RA - 282.25 / 180.0 * M_PI )) / cos( Gal[0] );
  
  // Longitude
  Gal[1] = atan2( SinGL, CosGL );
  Gal[1] += 33.0 / 180.0 * M_PI;
  if( verb > 9 )
    cout << "Sin: " << SinGL << " Cos: " << CosGL << endl;
  if( verb > 0 )
    cout << "Galactic coordinates: \n\t GB: " << Gal[0] / M_PI * 180.0 
         << "\n\t GL: " << Gal[1] / M_PI * 180.0 << endl;
}

long double GetKz( long double zz ){
  // The following is a linear interpolation of the Holmberg & Flynn (MNRAS 2004)
  // Figure 8 Kz-z relationship.
  //
  // Input is the height above the plane, in pc.
  // Output is the Kz force in m/s^2.

  long double convFactor = 1e6/PARSEC; // km^2/pc to m.
  
  if( zz < 22.5L )
    return( ( (0.15-0.03125)/(22.5-2.5)*zz+(0.03125*22.5-0.15*2.5)/(22.5-2.5) ) * convFactor );
  else if( zz < 60.0L )
    return( ( (zz-22.5)/(60-22.5)*0.325 + (60-zz)/(60-22.5)*0.15 ) * convFactor );
  else if( zz < 130.0L )
    return( ( (zz-60)/(130-60)*0.60625 + (130-zz)/(130-60)*0.325 ) * convFactor );
  else if( zz < 202.5L )
    return( ( (zz-130)/(202.5-130)*0.8125 + (202.5-zz)/(202.5-130)*0.60625 ) * convFactor );
  else if( zz < 280.0L )
    return( ( (zz-202.5)/(280-202.5)*1.0 + (280-zz)/(280-202.5)*0.8125 ) * convFactor );
  else if( zz < 365.0L )
    return( ( (zz-280)/(365-280)*1.15 + (365-zz)/(365-280)*1.0 ) * convFactor );
  else if( zz < 452.5 )
    return( ( (zz-365)/(452.5-365)*1.2875 + (452.5-zz)/(452.5-365)*1.15 ) * convFactor );
  else if( zz < 542.5 )
    return( ( (zz-452.5)/(542.5-452.5)*1.4125 + (542.5-zz)/(542.5-452.5)*1.2875 ) * convFactor );
  else if( zz < 690 )
    return( ( (zz-542.5)/(690-542.5)*1.56875 + (690-zz)/(690-542.5)*1.4125 ) * convFactor );
  else if( zz < 805 )
    return( ( (zz-690)/(805-690)*1.675 + (805-zz)/(805-690)*1.56875 ) * convFactor );
  else if( zz < 905 )
    return( ( (zz-805)/(905-805)*1.75625 + (905-zz)/(905-805)*1.675 ) * convFactor );
  else if( zz < 1027.5 )
    return( ( (zz-905)/(1027.5-905)*1.85 + (1027.5-zz)/(1027.5-905)*1.75625 ) * convFactor );
  else if( zz < 1145 )
    return( ( (zz-1027.5)/(1145-1027.5)*1.94375 + (1145-zz)/(1145-1027.5)*1.85 ) * convFactor );
  else if( zz < 1227.5 )
    return( ( (zz-1145)/(1227.5-1145)*1.99375 + (1227.5-zz)/(1227.5-1145)*1.94375 ) * convFactor );
  else if( zz < 1287.5 )
    return( ( (zz-1227.5)/(1287.5-1227.5)*2.0375 + (1287.5-zz)/(1287.5-1227.5)*1.99375 ) * convFactor );
  else if( zz < 1385 )
    return( ( (zz-1287.5)/(1385-1287.5)*2.09375 + (1385-zz)/(1385-1287.5)*2.0375 ) * convFactor );
  else 
    return( ( (zz-1385)/(1482.5-1385)*2.1625 + (1482.5-zz)/(1482.5-1385)*2.09375 ) * convFactor );
  
}

void CheckRange( InputParams IP, Param *XYZ, pulsar *psr ){
  int ct; 
  int parDef = IP.pd();
  // Check if we are in a positive M_psr range:
  long double kinMin = -1.0;
  long double kinMax = -1.0;
  long double M2Min = -1.0;
  int H3ID = 0;
  int Defined = 0;
  for( ct = 0; ct < parDef; ct++ ){
    if( strcmp( XYZ[ct].name(), "KIN" ) == 0 ){
      // Have KIN
      kinMin = XYZ[ct].I();
      kinMax = XYZ[ct].A();
      Defined += 1;
    }else if( strcmp( XYZ[ct].name(), "M2" ) == 0 ){
      // Have M2
      M2Min = XYZ[ct].I();
      Defined += 2;
    }else if( strcmp( XYZ[ct].name(), "COSI" ) == 0 ){
      // Have COSI
      kinMin = acos( XYZ[ct].A() ) / M_PI * 180.0;
      kinMax = acos( XYZ[ct].I() ) / M_PI * 180.0;
      Defined += 4;
    }else if( strcmp( XYZ[ct].name(), "H3" ) == 0 ){
      // Have H3
      // Will have to calculate the minimum M2 value based on H3 and COSI
      // ergo: we'll do that after this loop.
      H3ID = ct;
      Defined += 8;
    }
  }

  if( Defined == 12 ){
    long double stigMax;
    // Have COSI and H3
    // Get minimum M2 value: 
    // minimum H3 and maximum stig:
    if( kinMin < 90.0 && kinMax > 90.0 ){
      // maximum stig value is 1, at KIN = 90.
      stigMax = 1.0;
    }else if( kinMax < 90.0 ){
      // only in first quadrant -- maximum stig value at kinMax
      stigMax = sin( kinMax / 180.0 * M_PI ) 
        / ( 1 + sqrt( 1 - sin( kinMax / 180.0 * M_PI ) ) );
    }else if( kinMin > 90.0 ){ 
      // only in second quadrant -- maximum stig value at kinMin
      stigMax = sin( kinMin / 180.0 * M_PI )
        / ( 1 + sqrt( 1 - sin( kinMax / 180.0 * M_PI ) ) );
    }else{
      cout << "ERROR! Unphysical inclination angle range!\n";
      exit( 1 );
    }
    M2Min = XYZ[ct].I() / TSUN / pow( stigMax, 3.0 );
  }

  if( Defined == 0 ){
    // Have none of M2, KIN, COSI, H3
    // Ergo: nothing to check for.
  }else if( Defined == 3 || Defined == 12 ){
    // Have M2 and KIN or H3 and COSI

    // Check minimum inclination angle:
    if( M2Min * 4.925490947e-6 * 
        pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
             24.0 * 3600.0, 2.0 ) *
        pow( sin( kinMin * M_PI / 180.0 ) / 
             psr[0].param[param_a1].val[0], 3.0 )  <= 1 ){
      long double minM2 = 1.0 /( 4.925490947e-6 * 
                                 pow( psr[0].param[param_pb].val[0] /2.0 / M_PI
                                      * 24.0 * 3600.0, 2.0 ) *
                                 pow( sin( kinMin * M_PI / 180.0 ) / 
                                      psr[0].param[param_a1].val[0], 3.0 ));
      long double minKin = 180.0 / M_PI * 
        asin( pow( 4.925490947e-6 * 
                   pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
                        24.0 * 3600.0, 2.0 ) *
                   pow( psr[0].param[param_a1].val[0], -3.0 ) * 
                   M2Min, -1.0/3.0 ) );
      
      if( Defined == 3 ){
        cout << "ERROR! Your mass range doesn't work out!" << endl
             << "       Please increase your companion mass " << endl
             << "       to at least " << minM2
             << " Or increase your minimum KIN value to: " 
             << minKin << endl;
      }else if( Defined == 12 ){
        long double minStig = sin( minKin / 180.0 * M_PI )
          / ( 1 + sqrt( 1 - sin( minKin / 180.0 * M_PI ) ) );
        cout << "ERROR! Your mass range doesn't work out!\n"
             << "       Please decrease your maximum COSI value to at least "
             << cos( minKin / 180.0 * M_PI )
             << "      Or increase your minimum H3 value to " 
             << minM2 * TSUN * pow( minStig, 3.0 ) << endl;
        exit( 1 );
      }else{
        cout << "ERROR! This situation should not occur (0).\n";
        exit( 1 );
      }
    } 
    
    // Check maximum inclination angle:
    if( M2Min * 4.925490947e-6 * 
        pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
             24.0 * 3600.0, 2.0 ) *
        pow( sin( kinMax * M_PI / 180.0 ) / 
             psr[0].param[param_a1].val[0], 3.0 )  <= 1 ){
      long double minM2 = 1.0 /( 4.925490947e-6 * 
                                 pow( psr[0].param[param_pb].val[0] /2.0 / M_PI
                                      * 24.0 * 3600.0, 2.0 ) *
                                 pow( sin( kinMax * M_PI / 180.0 ) / 
                                      psr[0].param[param_a1].val[0], 3.0 ));
      long double maxKin = 180.0 - 180.0 / M_PI * 
        asin( pow( 4.925490947e-6 * 
                   pow( psr[0].param[param_pb].val[0] /2.0 / M_PI *
                        24.0 * 3600.0, 2.0 ) *
                   pow( psr[0].param[param_a1].val[0], -3.0 ) * 
                   M2Min, -1.0/3.0 ) );
      if( Defined == 3 ){
        cout << "ERROR! Your mass range doesn't work out!." << endl
             << "       Please increase your companion mass " << endl
             << "       to at least " << minM2
             << " Or decrease your maximum KIN value to:"
             << maxKin << endl;
        exit( 1 );
      }else if( Defined == 12 ){
        long double minStig = sin( maxKin / 180.0 * M_PI )
          / ( 1 + sqrt( 1 - sin( maxKin / 180.0 * M_PI ) ) );
        cout << "ERROR! Your mass range doesn't work out!\n"
             << "       Please increase your minimum COSI value to at least "
             << cos( maxKin / 180.0 * M_PI )
             << "      Or increase your minimum H3 value to " 
             << minM2 * TSUN * pow( minStig, 3.0 ) << endl;

      }else{
        cout << "ERROR! This situation should not occur (1).\n";
        exit( 1 );
      }
    }
  }else{
    cout << "WARNING! You seem to have a mixture of COSI/KIN/M2/H3.\n"
         << "         Normally you use either M2-KIN or COSI-H3.\n"
         << "         I suggest you press Ctrl-c and reconsider NOW.\n"
         << "         Alternatively, you can proceed at your own risk.\n";
    GiveMeABreak();
  }
}
