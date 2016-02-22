#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include <stdlib.h>

#include <iostream>
#include <iomanip>

/* Timing model      */
/* Based on bnrydd.f */

// A Tempo2 plug-in to map the chi-squared space of pulsar timing data.
// Joris P.W. Verbiest, MPIfR Bonn, 28 March 2012
// Transformed to a timing model by Gregory Desvignes MPIfR Bonn, Feb. 2016

using namespace std;

#define DAY2S 86400

int verb = 0; // generic verbosity flag
int eq2gal(float ra, float dec, double *glon,double *glat);
int SetPbDot( pulsar *psr, long double dist);
void SetOmDot( pulsar *psr );

longdouble D3model(pulsar *psr,int p,int ipos,int param)
{
  long double an;
  long double pb,k;
  long double rad2deg = 180.0/M_PI;
  long double SUNMASS = 4.925490947e-6;
  long double m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,am2,ct;
  long double pbdot,xpbdot,phase,u,du,gamma;
  long double orbits;
  int norbits;
  long double  cu,onemecu,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
  long double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb;
  long double csigma,ce,cx,comega,cgamma,cdth,cm2,csi;
  long double dist;
  const char *CVS_verNum = "$Revision: 1.7 $";

  if (displayCVSversion == 1) CVSdisplayVersion("D3model.C","D3model()",CVS_verNum);

  if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
  else si = 0.0;

  if (si > 1.0)
    {
      displayMsg(1,"BIN1","SIN I > 1.0, setting to 1: should probably use D3S model","",psr[p].noWarnings);
      si = 1.0;
      psr[p].param[param_sini].val[0] = 1.0L;
    }

  if (psr[p].param[param_dist].paramSet[0]==1) {
      dist = psr[p].param[param_dist].val[0];
      // Reset the value of parallax
      psr[p].param[param_px].paramSet[0] = 1;
      psr[p].param[param_px].fitFlag[0] = 0;
      psr[p].param[param_px].val[0] = 1000./ dist;
  } else {
      printf("ERROR! For the D3 model, DIST must be set.\n");
      exit( 1 );
  }

  SetPbDot(psr, dist);

  if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
  else am2 = 0.0;

  pb = psr[p].param[param_pb].val[0]*SECDAY;
  an = 2.0*M_PI/pb;
  k = psr[p].param[param_omdot].val[0]/(rad2deg*365.25*86400.0*an);

  m2 = am2*SUNMASS;
  t0 = psr[p].param[param_t0].val[0];
  ct = psr[p].obsn[ipos].bbat;    
  
  tt0 = (ct-t0)*SECDAY;

  if (psr[p].param[param_gamma].paramSet[0]==1)
    gamma = psr[p].param[param_gamma].val[0];
  else
    gamma = 0.0;

  if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0];
  else omz = 0.0;

  if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
  else xdot  = 0.0;

  if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
  else pbdot = 0.0;

  if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
  else edot = 0.0;

  if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
  else xpbdot = 0.0;

  if (psr[p].param[param_dr].paramSet[0] == 1) dr = psr[p].param[param_dr].val[0];
  else dr = 0.0;

  if (psr[p].param[param_dtheta].paramSet[0] == 1) dth = psr[p].param[param_dtheta].val[0];
  else dth = 0.0;

  if (psr[p].param[param_a0].paramSet[0] == 1) a0 = psr[p].param[param_a0].val[0];
  else a0 = 0.0;

  if (psr[p].param[param_b0].paramSet[0] == 1) b0 = psr[p].param[param_b0].val[0];
  else b0 = 0.0;


  x = psr[p].param[param_a1].val[0]+xdot*tt0;
  ecc = psr[p].param[param_ecc].val[0]+edot*tt0;
  er = ecc*(1.0+dr);
  eth = ecc*(1.0+dth);

  if (ecc < 0.0 || ecc > 1.0)
    {
      printf("D3model: problem with eccentricity = %Lg [%s]\n",psr[p].param[param_ecc].val[0],psr[p].name);
      exit(1);
    }

  orbits = tt0/pb - 0.5L*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb);
  norbits = (int)orbits;
  if (orbits<0.0) norbits--;
  phase=2.0*M_PI*(orbits-norbits);
  /*  Compute eccentric anomaly u by iterating Kepler's equation. */
  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
  do {
    du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
    u=u+du;
  } while (fabs(du)>1.0e-12);
  
  /*  D3 equations 17b, 17c, 29, and 46 through 52 */
  su=sin(u);
  cu=cos(u);
  onemecu=1.0-ecc*cu;
  cae=(cu-ecc)/onemecu;
  sae=sqrt(1.0-pow(ecc,2))*su/onemecu;
  ae=atan2(sae,cae);
  if(ae<0.0) ae=ae+2.0*M_PI;
  ae=2.0*M_PI*orbits + ae - phase;
  omega=omz/rad2deg + k*ae;
  sw=sin(omega);
  cw=cos(omega);
  alpha=x*sw;
  beta=x*sqrt(1-pow(eth,2))*cw;
  bg=beta+gamma;
  dre=alpha*(cu-er) + bg*su;
  drep=-alpha*su + bg*cu;
  drepp=-alpha*cu - bg*su;
  anhat=an/onemecu;

  /* D3 equations 26, 27, 57: */
  sqr1me2=sqrt(1-pow(ecc,2));
  cume=cu-ecc;
  brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
  //  printf("GEORGE: si = %g, brace = %g\n",(double)si,(double)brace);
  dlogbr=log(brace);
  ds=-2*m2*dlogbr;
  da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw);

  /*  Now compute d2bar, the orbital time correction in D3 equation 42. */
  d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
					  0.5*ecc*su*dre*drep/onemecu)) + ds + da;
  torb=-d2bar;

  if (param==-1) return torb;
  
  /*  Now we need the partial derivatives. Use D3 equations 62a - 62k. */
  csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
  ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
  cx=sw*cume+sqr1me2*cw*su;
  comega=x*(cw*cume-sqr1me2*sw*su);
  cgamma=su;
  cdth=-ecc*ecc*x*cw*su/sqr1me2;
  cm2=-2*dlogbr;
  csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 
  if (param==param_pb)
    return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
  else if (param==param_a1)
    return cx;
  else if (param==param_ecc)
    return ce;
  else if (param==param_edot)
    return ce*tt0;
  else if (param==param_om)
    return comega;
  else if (param==param_omdot)
    return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
  else if (param==param_t0)
    return -csigma*an*SECDAY;
  else if (param==param_pbdot)
    return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
  else if (param==param_sini)
    return csi;
  else if (param==param_gamma)
    return cgamma;
  else if (param==param_m2)
    return cm2*SUNMASS;
  else if (param==param_a1dot) /* Also known as xdot */
    return cx*tt0;
  else if (param==param_dtheta)
    return cdth;

  return 0.0;
}


void updateD3(pulsar *psr,double val,double err,int pos)
{
  if (pos==param_pb)
    {
      psr->param[param_pb].val[0] += val/SECDAY;
      psr->param[param_pb].err[0]  = err/SECDAY;
    }
  else if (pos==param_a1 || pos==param_ecc || pos==param_t0 || pos==param_sini || pos==param_m2
	   || pos == param_gamma || pos==param_edot || pos==param_dist)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_om)
    {
      psr->param[pos].val[0] += val*180.0/M_PI;
      psr->param[pos].err[0]  = err*180.0/M_PI;
    }
  else if (pos==param_pbdot || pos==param_dtheta)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_a1dot)
    {
      psr->param[pos].val[0] += val;
      psr->param[pos].err[0]  = err;
    }
  else if (pos==param_omdot)
    {
      psr->param[pos].val[0] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
      psr->param[pos].err[0]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
    }
}

int SetPbDot( pulsar *psr, long double dist ){

  // Now determine companion mass and inclination angle:
  dist /= 1000.;

  long double m2=-1.0, sini = -2.0;
  if(psr[0].param[param_h3].paramSet[0] == 1 && psr[0].param[param_stig].paramSet[0] == 1 ){
    // Using DDH model (H3, stig)
    m2 = psr[0].param[param_h3].val[0] / TSUN 
      / pow( psr[0].param[param_stig].val[0], 3.0 );
    sini = 2.0 * psr[0].param[param_stig].val[0]
      / ( 1.0 + pow( psr[0].param[param_stig].val[0], 2.0 ) );
  }else if(psr[0].param[param_h3].paramSet[0] == 1 && psr[0].param[param_h4].paramSet[0] == 1 ){
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
    Kinematic = dist * 1000.0 * PCM // distance in meters
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

  double gl, gb; // Galactic coordinates (0: lat; 1: long)
  //GetGalCoord( psr[0].name, Gal );
  // TODO
  eq2gal(psr[0].param[param_raj].val[0], psr[0].param[param_decj].val[0], &gl, &gb);

  // Kz (vertical Galactic potential):
  long double Kz;
  Kz = - sin( gb ) * psr[0].param[param_pb].val[0] * 24.0 * 3600.0 
    / 299792458 // speed of light in meters per second
    * Kz_HF04( dist * 1000.0 * sin( gb ) );

  // DGR (differential Galactic rotation):
  long double DGR;
  long double beta = dist/8.34 * cos( gb ) - cos( gl );
  DGR = -pow( 240e3, 2.0 ) // Galactic rotation in m/s
    /299792458 // speed of light in m/s
    /8340. / PCM // Galactic radius in m
    * cos( gb ) * ( cos( gl ) + 
        beta / ( pow( sin( gl ),2.0 ) + pow( beta, 2.0 ) ) ) *
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
      std::cout << "Have Mpsr: " << Mpsr << std::endl; fflush( stdout );
    if( Mpsr < 0.0 || isnan(Mpsr)){
      std::cout << "ERROR! Used:\n"
           << " m2: " << m2 << std::endl
           << " pb: " << psr[0].param[param_pb].val[0] << std::endl
           << " a1: " << psr[0].param[param_a1].val[0] << std::endl
           << " sini: " << sini << std::endl
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
           << std::endl;
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
    std::cout << "GD Contributions to Pbdot: " 
         << "Kinematic: " << std::setprecision( 10 ) << Kinematic
         << "\tKz: " << std::setprecision( 10 ) << Kz 
         << "\tDGR: " << std::setprecision( 10 ) << DGR
         << "\tRelativistic: " << std::setprecision( 10 ) << relativistic
         << std::endl; fflush( stdout );
  if( verb > 9 )
    std::cout << "Setting Pbdot to: " << std::setprecision( 10 )
         << Kinematic + Kz + DGR + relativistic << std::endl; fflush( stdout );

  
  if( psr[0].param[param_pbdot].paramSet[0] != 1 ){
    // initialise pbdot if needed. (NOT CHECKED JORIS );
    psr[0].param[param_pbdot].paramSet[0] = 1;
    strcpy( psr[0].param[param_pbdot].label[0], "PBDOT" );
    strcpy( psr[0].param[param_pbdot].shortlabel[0], "PBDOT" );
  }
  psr[0].param[param_pbdot].val[0] = Kinematic + Kz + DGR + relativistic;
  psr[0].param[param_pbdot].fitFlag[0] = 0;
  psr[0].param[param_pbdot].err[0] = 0.0;
  return 0;
}

/* J2000 coordinates of galactic north pole and x-axis */
#define N_GAL_POLE_RA    (12 + (51 + (26.2754/60))/60)*15
#define N_GAL_POLE_DEC   (27 + (07 + (41.705/60))/60)
#define GAL_CENTER_DEC  -(28 + (56 + (10.219/60))/60)

int eq2gal(float ra, float dec, double *glon,double *glat) { 
  double raddeg, spdec, cpdec, glong_off, ra_off;
  double sra, cra, sdec, cdec;

  raddeg = M_PI/180.0;

  spdec = sin((90 - N_GAL_POLE_DEC)*raddeg);
  cpdec = cos((90 - N_GAL_POLE_DEC)*raddeg);
  glong_off = -asin(sin(GAL_CENTER_DEC*raddeg)/spdec)/raddeg;
  ra_off = 90 + N_GAL_POLE_RA;

  sra = sin((ra - ra_off)*raddeg);
  cra = cos((ra - ra_off)*raddeg);
  sdec = sin(dec*raddeg);
  cdec = cos(dec*raddeg);

  *glon = glong_off;
  *glon = glong_off + atan2((cdec*sra*cpdec + sdec*spdec), cdec*cra)/raddeg;
  if (*glon < 0) *glon = *glon + 360;
  *glat = asin(sdec*cpdec - cdec*sra*spdec)/raddeg;

  return(0);
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



#if 0
void cel2gal(int rah, int ram, double ras,int decd,int decm,double decs,double *glon,double *glat) /*          includefile*/
{
double ra, dec, gl, gb;
int isign=1;
double twopi, raddeg, spdec, cpdec, glong_off, ra_off;
double sra, cra, sdec, cdec;

    ra = (rah + (ram + (ras/60))/60)*15;
    if (decd != 0) isign = decd/abs(decd);
    dec = decd + isign*(decm + (decs/60))/60;

  twopi = 4*acos(0.0);
  raddeg = twopi/360.0;

  spdec = sin((90 - N_GAL_POLE_DEC)*raddeg);
  cpdec = cos((90 - N_GAL_POLE_DEC)*raddeg);
  glong_off = -asin(sin(GAL_CENTER_DEC*raddeg)/spdec)/raddeg;
  ra_off = 90 + N_GAL_POLE_RA;

  sra = sin((ra - ra_off)*raddeg);
  cra = cos((ra - ra_off)*raddeg);
  sdec = sin(dec*raddeg);
  cdec = cos(dec*raddeg);

  *glon = glong_off;
  *glon = glong_off + atan2((cdec*sra*cpdec + sdec*spdec), cdec*cra)/raddeg;
  if (*glon < 0)
    *glon = *glon + 360;
  *glat = asin(sdec*cpdec - cdec*sra*spdec)/raddeg;

}
#endif

