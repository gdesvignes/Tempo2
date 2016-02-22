#include "tempo2.h"

//  Implemented by Joris Verbiest

long double Kz_HF04( long double zz ){
  // The following is a linear interpolation of the Holmberg & Flynn (MNRAS 2004)
  // Figure 8 Kz-z relationship.
  //
  // Input is the height above the plane, in pc.
  // Output is the Kz force in m/s^2.

  long double convFactor = 1e6/PCM; // km^2/pc to m.
  
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

