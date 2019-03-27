#include <math.h>
#include <cmath>
#include <iostream>

using namespace std;

int main(){
  long double ro=3430.18601252541;  ro=3430.18601252541d0	! equatorial component in km - Irbene
  long double albypi=0.121415123729; albypi=0.121415123729d0	! along/pi                   - Irbene
  long double vhor=2*M_PI*ro/(24*3600)*1.002737909350795; vhor=2*pi*ro/(24*3600)*1.002737909350795d0	!km/s
  long double t = 6987.8545138887130;
  long double aLMST=M_PI*fmod(5.55811454652+fmod(t+t,(long double) 2.0)+ t*(.547581870159-2.0+t*(1.61549-15-t*1.473-24))+albypi, (long double) 2.0);
              aLMST=pi*dmod(5.55811454652d0+dmod(t+t,2d0)+ t*(.547581870159d0-2+t*(1.61549d0-15-t*1.473d-24))+albypi,2d0)
  long double ravhor=aLMST+M_PI/2.0; ravhor=aLMST+pi/2
  long double dec = 62.0 + 1.0/60.0 + 46.7/3600.0;  
  long double cdec=cos(dec* M_PI / 180.0 ); cdec=dcos(dec)
  long double ra = 22.0 + 56.0/60.0 + 17.90/3600.0;
  long double vobs=vhor*cdec*cos((ra*M_PI/180.0)-ravhor); vobs=vhor*cdec*dcos(ra-ravhor)
  cout<<"aLMST "<<(long double) aLMST<<endl;
  cout<<"ravhor "<<(long double) ravhor<<endl;
  cout<<"cdec "<<(long double) cdec<<endl;
  cout<<"vobs "<<(long double) vobs<<endl;  

  return 0;
}

	
	
	
	
	
	
	
