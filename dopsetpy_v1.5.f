c dopset for Irbene 32m antenna

	implicit real*8 (a-h,o-z)
	character*(*) argument(12)*20


        pi=3.1415926535898d0
        c=2.99792458d+5 

       l_arg=iargc()
        
       if (l_arg.eq.12) then 

 
       do i=1,l_arg
          call getarg(i,argument(i))
       enddo

        read(argument(1),*)date7
        read(argument(2),*)im
        read(argument(3),*)date5

        read(argument(4),*)iUThh
        read(argument(5),*)iUTmm
        read(argument(6),*)iUTss

        read(argument(7),*)rahh
        read(argument(8),*)ramm
        read(argument(9),*)rass

        read(argument(10),*)decdd
        read(argument(11),*)decmm
        read(argument(12),*)decss




ccc inputs:
c        date7=2014               !real year
c        im=12                    !int month
c        date5=21                !real day
        
c        iUTss=0
c        iUTmm=10
c        iUThh=10

c        rahh=12
c        ramm=0
c        rass=0
        
c        decdd=50
c        decmm=0
c        decss=0
        



        Epoch=2000.0


cccccccc Compute the Doppler shift cccccccccccccccccccccc
c Julian Day at given UT
	DJ = JD(int(date7),im,int(date5),0) - .5d0
     &       + ((iUTss/60d0 + iUTmm)/60d0 + iUThh)/24d0
	rar0  = (rahh + ramm/60d0 + rass/3600d0)/12*pi
	decr0 = (dabs(decdd)+dabs(decmm)/60d0+dabs(decss)/36d2)/180*pi
	if(decdd.lt.0d0.or.decmm.lt.0d0.or.decss.lt.0d0) decr0=-decr0

c       write(*,*)'rar0= ',rar0,'  decr0= ',decr0

c Find 'epoch' for precession and nutation
		if(int(Epoch*100).eq.195000) then
	   if(DJ.gt.2454476d0) then
	DJ0 = 2433282.42345905d0	! That is what B1950.0 means
	   else   ! Before Jan 10, 2008 our 'epoch' was fake and
	DJ0 = DJ  ! the coordinates already precessed to date
	Epoch = 2000 + (DJ0 - 2451545d0)/365.25d0
	   endif
		else
	DJ0 = JD(int(Epoch),1,1,0) - .5d0 + dmod(Epoch,1d0)*365.25d0
		endif
c Precess and nutate the equatorial coordinates (precessing over
c 50 years makes some 0.5 km/s difference in velocity)

	call PREnew(rar0,decr0,rarp,decrp,DJ0,DJ)
c        write(*,*)'rarp= ',rarp,'  decrp= ',decrp


	call nutate(DJ,rarp,decrp,rar,decr)

c        write(*,*)'rar= ',rar,'  decr= ',decr



c The following subroutine computes the Earth velocity that agrees
c with the JPL ephemeris within 0.0006 km/s (checked for 1990 - 2019)
	call VdRT32(dj, rar, decr, aLMST, Vsun, Vobs, Vdop)
c         Units:  d   rad  rad   rad    km/s  km/s  km/s
c Vdop is the quantity of interest - the projection of the location
c (RT32) velocity with respect to the LSR onto the source direction
cccccccccccccccccccccccccccccccccccccccccccccccccc
        freq0=6668.5192
        f_shift=(Vdop/c)*freq0
    
        write(*,*)'----------------------------------'
        write(*,*)'Julian Day at given UT: ', dj
        write(*,*)'V Sun: ', Vsun
        write(*,*)'V obs: ', Vobs
        write(*,*)'V tot: ', Vdop
        write(*,*)'f shift at ', freq0,'MHz  = ', f_shift,' MHz'
        write(*,*)'----------------------------------'

        open(11,file="V_tot.dat")
        write(11,*)Vdop,0.0
        close(11)

        open(11,file="lsrShift.dat")
        write(11,995) "DateTime",int(date7), im, int(date5),
     & iUThh,iUTmm,iUTss
c     MJD=JD-2400000.5d0   ! this formula gives Modified Julian Day number
        write(11,996) "LSRshift",f_shift,"MHz"
        write(11,997) "MJD",dj-2400000.5d0
        write(11,997)'Vobs', Vobs,"km/s"
        write(11,996) "AtFreq",freq0,"MHz"
        write(11,996) "FreqShift",f_shift,"MHz"
        write(11,996) "VelTotal",Vdop,"km/s"
        close(11)

        else  
        write(*,*)'Input error.'
        write(*,*)'Syntax: dopsetpy yyyy mm dd  UT_h UT_m UT_s RA_h RA_m
     &RA_s DEC_d DEC_m DEC_s'

995     format (a,";"i4,"-",i0.2,"-",i0.2,1x,i0.2,":",i0.2,":",i0.2,";")
996     format (a,";",f15.10,";",a,";")
997     format (a,";",f10.4,";",a,";")


        endif
        
        END

c====================================================================
c====================================================================
c====================================================================


	SUBROUTINE PREnew(dra,d,dra1,d1,Dje1,Dje2)

C CALCULATES GENERAL PRECESSION FROM Dje1 TO Dje2 (new IAU1976 theory)

	implicit real*8 (d,r)
	DATA dpi/3.141592653589793d0/
	CALL PRE(Dje1,Dje2,dzeta,dz,dth)
	dc=cos(dra+dzeta)
	dra1=dz+datan2(sin(dra+dzeta),-tan(d)*sin(dth)+cos(dth)*dc)
	d1=asin(sin(d)*cos(dth)+cos(d)*sin(dth)*dc)
	dra1=dmod(dra1+2*dpi,2*dpi)
	end


	SUBROUTINE nutate(DJ,ra,d,ra1,d1)
c Poprawia rektascensje (ra) i deklinacje (d) na efekt nutacji.
c Na wyjsciu ra1 i d1 sa wspolrzednymi poprawionymi [rad].
	implicit real*8 (d,e,r)
	call nutation(dj,dp,de,e)
	dce=dcos(e)
	dse=dsin(e)
	dsd=dsin(d)
	dcd=dcos(d)
	dsa=dsin(ra)
	dsb=dce*dsd-dse*dcd*dsa
	dcb=dsqrt(1-dsb*dsb)
	dl=datan2(dse*dsd+dce*dcd*dsa,dcd*dcos(ra))+dp
	dce=dcos(e+de)
	dse=dsin(e+de)
	dsl=dsin(dl)
	d1=dasin(dce*dsb+dse*dcb*dsl)
	ra1=datan2(-dse*dsb+dce*dcb*dsl,dcb*dcos(dl))
	if(ra1.lt.0.) ra1=ra1+2.*3.141592653589793d0
	end

	subroutine nutation(dj,dpsi,deps,eps)
c oblicza nutacje w dlugosci (dpsi) i nachyleniu ekliptyki do rownika
c (deps) oraz nachylenie (eps, w radianach) na date julianska dj.
	implicit real*8 (d,e,t)
	integer*2 k,k1,k2
	integer*4 ide
	dimension da(5),k(5,106),k1(5,61),k2(5,45),ide(2,106),bd(15)
     *,be(15)
	equivalence (k(1,1),k1(1,1)),(k(1,62),k2(1,1))
	data bd/-174.2,.2,-1.6,-3.4,1.2,-.5,.1,-.1,.1,-.2,.1,-.4,0.,.1,-.1
     */,be/8.9,.5,-3.1,-.1,-.6,.3,0.,0.,0.,-.5,0.,0.,-.1,2*0./,ide/
     *-171996,92025,2062,-895,-13187,5736,1426,54,-517,224,217,-95,
     *129,-70,17,0,-16,7,-2274,977,712,-7,-386,200,-301,129,63,-33,
     *-58,32,     46,-24,11,0,-3,1,-3,0,-2,1,1,0,48,1,-22,0,-15,9,
     *-12,6,-6,3,-5,3,4,-2,4,-2,-4,0,1,0,1,0,-1,0,1,0,1,0,-1,0,
     *-158,-1,123,-53,63,-2,-59,26,-51,27,-38,16,29,-1,29,-12,-31,13,
     *26,-1,21,-10,16,-8,-13,7,-10,5,-7,0,7,-3,-7,3,-8,3,6,0,6,-3,-6,3,
     *-7,3,6,-3,-5,3,5,0,-5,3,-4,0,4,0,-4,0,-3,0,3,0,-3,1,-3,1,-2,1,
     *-3,1,-3,1,2,-1,-2,1,2,-1,-2,1,2,0,2,-1,1,-1,-1,0,1,-1,-2,1,-1,0,
     *1,-1,-1,1,-1,1,1,0,1,0,1,-1,-1,0,-1,0,1,0,1,0,-1,0,1,0,1,0,
     *-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,-1,0,1,0/
	data k1/4*0,1, 4*0,2, 0,0,2,-2,2, 0,1,4*0,1,2,-2,2 ,0,-1,2,-2,2,
     c0,0,2,-2,1, 0,2,4*0,2,2,-2,2, 0,0,2,0,2, 1,6*0,2,0,1, 1,0,2,0,2,
     c1,3*0,1,-1,3*0,1,    -2,0,2,0,1, 2,0,-2,0,0, -2,0,2,0,2, 1,-1,0,-1
     c,0,0,-2,2,-2,1, 2,0,-2,0,1, 2,0,0,-2,3*0,2,-2,0, 0,1,0,0,1,0,-1,0,
     c0,1,-2,0,0,2,1,0,-1,2,-2,1,2,0,0,-2,1,0,1,2,-2,1,1,0,0,-1,0,2,1,0,
     c-2,0,0,0,-2,2,1,0,1,-2,2,0, 0,1,0,0,2, -1,0,0,1,1, 0,1,2,-2,0,
     c1,0,0,-2,0,-1,0,2,0,2,3*0,2,0,  -1,0,3*2,1,0,2,0,1,0,0,4*2,4*0,
     c1,0,2,-2,2,2,0,2,0,2,0,0,2,0,0, -1,0,2,0,1, -1,0,0,2,1, 1,0,0,-2,1
     c, -1,0,2,2,1, 1,1,0,-2,0, 0,1,2,0,2, 0,-1,2,0,2, 1,0,3*2, 1,0,0,2,
     c0, 2,0,2,-2,2, 3*0,2,1, 0,0,2,2,1, 1,0,2,-2,1, 3*0,-2,1, 1,-1,3*0/
c62
     *, k2/2,0,2,0,1, 0,1,0,-2,0, 1,0,-2,5*0,1,0, 1,1,3*0, 1,0,2,0,0,
     c1,-1,2,0,2, 2*-1,3*2, -2,3*0,1, 3,0,2,0,2, 0,-1,3*2, 1,1,2,0,2,
     c-1,0,2,-2,1, 2,3*0,1, 1,3*0,2, 3,6*0,2,1,2, -1,3*0,2, 1,0,0,-4,0,
     c-2,0,3*2,-1,0,2,4,2, 2,0,0,-4,0, 1,1,2,-2,2, 1,0,2,2,1, -2,0,2,4,2
     c,-1,0,4,0,2, 1,-1,0,-2,0, 2,0,2,-2,1,2,0,3*2, 1,0,0,2,1, 0,0,4,-2,
     c2, 3,0,2,-2,2, 1,0,2,-2,0, 0,1,2,0,1, 2*-1,0,2,1, 0,0,-2,0,1,
     c0,0,2,-1,2, 0,1,0,2,0, 1,0,2*-2,0, 0,-1,2,0,3*1,0,-2,1, 1,0,-2,2,0
     c, 2,0,0,2,3*0,2,4,2, 0,1,0,1,0/
	df(d1,i2,d2,d3,d4)=(d1+dmod(i2*t,1d0)*360*3600+
     *t*(d2+t*(d3+t*d4)))*dpi/(3600*180)
	dpi=3.141592653589793d0
	t=(dj-2451545d0)/36525d0
	da(1)=df( 485866.733d0,1325, 715922.633d0, 31.310d0,.064d0)
	da(2)=df(1287099.804d0,  99,1292581.224d0,  -.577d0,-.012d0)
	da(3)=df( 335778.877d0,1342, 295263.137d0,-13.257d0,.011d0)
	da(4)=df(1072261.307d0,1236,1105601.328d0, -6.891d0,.019d0)
	da(5)=df( 450160.28d0,   -5,-482890.539d0,  7.455d0,.008d0)
	dpsi=0
	deps=0
	cd=0
	ce=0
	do 3 i=106,1,-1
	dfaz=0
	do 1 j=1,5
1	dfaz=dfaz+k(j,i)*da(j)
	if(i.gt.15) go to 2
	cd=t*bd(i)
	ce=t*be(i)
2	dpsi=dpsi+(ide(1,i)+cd)*dsin(dfaz)
3	deps=deps+(ide(2,i)+ce)*dcos(dfaz)
	dpsi=dpsi/(10000d0*3600*180)*dpi
	deps=deps/(10000d0*3600*180)*dpi
	eps=df(84381.448d0,0,-46.815d0,-.00059d0,.001813d0)
	end





	subroutine VdRT32(dj, ra, dec, aLMST, Vsun, Vobs, Vtot)
c  Jednostki:            d rad   rad  rad    km/s  km/s km/s

c Oblicza rzut predkosci, Vtot [km/s], teleskopu RT32 (Piwnice
c k. Torunia) na kierunek obserwowanego zrodla o wspolrzednych
c     ra, dec - rektascensja i deklinacja [rad] przeprecesowane 
c na zadany moment czasu:
c     dj - data julianska (UTC)
c     aLMST - sredni czas gwiazdowy dla polozenia RT32 [rad]
c     Vsun - predkosc Slonca (ku apeksowi) jest uwzgledniona w Vtot
c dlatego, jesli nie chcemy jej uwzgledniac, trzeba uzyc Vtot-Vsun.
c Maksymalna roznica Vtot w stosunku do obliczen scislych z
c wykorzystaniem numerycznych efemeryd JPL w latach 1990 - 2019
c jest mniejsza niz 0.57 m/s

	implicit real*8 (a-h,o-z)
	real*8 VHelio(3),ve(3)
	real*4 dtuts,y
	data pi/3.141592653589793d0/

	ras=18*pi/12d0 ! Wspolrzedne apeksu Slonca  rekt. (epoch J1900)
	d=30*pi/180d0  ! 			    dekl.
	vSun0=20       ! km/s
	call PREold(ras,d,rav,dv,2415020.d0,dj) ! precesja tych wpolrz.
	cdec=dcos(dec)
	sdec=dsin(dec)
	vSun=vSun0*(dsin(dv)*sdec+dcos(dv)*cdec*dcos(ra-rav))

c x=3638.55851  y=1221.96972  z=5077.03676 <- wspolrz. RT32 [km]
c x=3183.661    y= 1276.902   z=5359.291   <- Irbene   RT32 coordinates [km]
c (wg.: P.Charlot et al., 2001, Proc. 15th Working Meet. Europ. VLBI
c for Geodesy and Astrometry, 194-200)
c ro=sqrt(x^2+y^2), along=atan2(y,x)
c ro=3838.27018685375, along=0.32400392423rad=pi*0.103133652246  for Torun telescope
c ro=3430.18601252541, along=0.381436860742553rad=pi*0.121415123729264  for Irbene telescope
c	ro=3838.2701869d0	! skladowa rownikowa w km    - Torun
c	albypi=0.103133652246d0	! along/pi                   - Torun
	ro=3430.18601252541d0	! equatorial component in km - Irbene
	albypi=0.121415123729d0	! along/pi                   - Irbene
	vhor=2*pi*ro/(24*3600)*1.002737909350795d0	!km/s
c vhor ma kierunek wschodni, tj. decl=0, ra=almst+pi/2
c sredni czas gwiazdowy miejsca
	t=dj-2451545.D0
	aLMST=pi*dmod(5.55811454652d0+dmod(t+t,2d0)+
     * t*(.547581870159d-2+t*(1.61549d-15-t*1.473d-24))+albypi,2d0)
	ravhor=aLMST+pi/2
	vobs=vhor*cdec*dcos(ra-ravhor)

	y=2000+t/365.25		! przyblizony rok
	dje=dj+dtuts(y)/86400		! czas efemeryd
	call BARVEL(DJE,0.d0,VHelio,ve)	
c ve zawiera predkosc wzgledem barycentrum Ukladu Slonecznego
	vEarth=(ve(1)*dcos(ra)+ve(2)*dsin(ra))*cdec+ve(3)*sdec

	Vtot = vSun + vEarth + vobs
	end


	SUBROUTINE PREold(dra,d,dra1,d1,Dje1,Dje2)

C CALCULATES GENERAL PRECESSION FROM Dje1 TO Dje2 (old IAU theory)

	IMPLICIT REAL*8 (D,t)
	DATA DCSAR/4.848136812D-6/,dpi/3.141592653589793d0/
     * ,DC1/2304.25D0/,DC2/1.396D0/,DC3/.302D0/,DC4/.018D0/,DC5/.791D0/
     * ,DC6/2004.683D0/,DC7/-.853D0/,DC8/-.426D0/,DC9/-.042D0/
	DT0=(Dje1-2415020d0)/36525
	DT=(Dje2-Dje1)/36525
	DTS=DT*DT
	DTC=DTS*DT
	DZETA=((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT=DZETA+DC5*DTS*DCSAR
	DTHET=((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSTHET=DSIN(DTHET)
	DCTHET=DCOS(DTHET)
	dra1=dzett+datan2(sin(dra+dzeta),-tan(d)*dsthet+
     + dcthet*cos(dra+dzeta))
	d1=asin(sin(d)*dcthet+cos(d)*dsthet*cos(dra+dzeta))
c correct the RA at the epoch of observations
	t=(dje2-2451545.)/36525
	dra1=dra1+(.0775+(.0851+.00002*t)*t)*dpi/(12*3600)
	if(dra1.lt.0d0) dra1=dra1+2*dpi
	END


	SUBROUTINE BARVEL(DJE,DEQ,DVELH,DVELB)

C Wg: P. Stumpff, 1980, Astron. Astrophys. Suppl. Ser. 41, 1-8.

C CALCULATES HELIOCENTRIC AND BARYCENTRIC VELOCITY COMPONENTS OF THE EARTH.
C THE LARGEST DEVIATIONS FROM THE JPL-DE96 ARE 42 CM/S.
c W latach 1999 -2020 stwierdzono odchylki od JPL DE405 <60 cm/s (KB)

C GIVEN  DJE = JULIAN EPHEMERIS DATE
C  	DEQ = EPOCH OF MEAN EQUATOR AND EQUINOX. IF DEQ=0, BOTH VECTORS
C	    DVELH AND DVELB ARE REFERRED TO MEAN EQUATOR AND EQUINOX 
C	    OF DATE (DJE).

C RESULT DVELH(K) = HELIOCENTRIC, DVELB(K) = BARYCENTRIC VELOCITY COMPONENTS.
C 	    (K=1,2,3 --- DX/DT,DY/DT,DZ/DT; UNIT=KM/S [Bar] i A.U./S [Hel])

	IMPLICIT REAL*8 (D)
	DIMENSION DVELH(3),DVELB(3),SN(4)
	DIMENSION DCFEL(3,8),DCEPS(3),CCSEL(3,17),DCARGS(2,15),
     *	CCAMPS(5,15),CCSEC(3,4),DCARGM(2,3),CCAMPM(4,3),CCPAMV(4)
  	EQUIVALENCE (SORBEL(1),E),(FORBEL(1),G)

	COMMON/BARXYZ/DPREMA(3,3),DPSI,D1PDRO,DSINLS,DCOSLS,DSINEP,DCOSEP,
     * FORBEL(7),SORBEL(17),SINLP(4),COSLP(4),SINLM,COSLM,SIGMA,IDEQ

	DATA DC2PI/6.2831853071796D0/,CC2PI/6.283185/,
     * DC1/1.0D0/,DCT0/2415020.D0/,DCJUL/36525D0/,DCBES/.313D0/,
     * DCTROP/365.24219572D0/,DC1900/1900.D0/

C CONSTANTS DCFEL(I,K) OF FAST CHANGING ELEMENTS

	DATA DCFEL/1.7400353D+00, 6.2833195099091D+02, 5.2796D-6,
     *		 6.2565836D0,   6.2830194572674D2,  -2.6180D-6,
     *		 4.7199666D0,   8.3997091449254D3,  -1.9780D-5,
     *		 1.9636505D-1,  8.4334662911720D3,  -5.6044D-5,
     *		4.1547339D0,    5.2993466764997D1,   5.8845D-6,
     *		4.6524223D0,    2.1354275911213D1,   5.6797D-6,
     *		4.2620486D0,    7.5025342197656D0,   5.5317D-6,
     *		1.4740694D0,    3.8377331909193D0,   5.6093D-6/

	DATA DCEPS/4.093198D-1,-2.271110D-4,-2.860401D-8/
	DATA CCSEL/1.675104E-2,-4.179579E-5,-1.260516E-7,
     *		2.220221E-1, 2.809917E-2, 1.852532E-5,
     *		1.589963E0,  3.418075E-2, 1.4302E-5,
     *		2.994089E0,  2.590824E-2, 4.155840E-6,
     *		8.155457E-1, 2.486352E-2, 6.836840E-6,
     *		1.735614E0,  1.763719E-2, 6.370440E-6,
     *		1.968564E0,  1.524020E-2,-2.517152E-6,
     *		1.282417E0,  8.703393E-3, 2.289292E-5,
     *		2.280820E0,  1.918010E-2, 4.484520E-6,
     *		4.833473E-2, 1.641773E-4,-4.654200E-7,
     *		5.589232E-2,-3.455092E-4,-7.388560E-7,
     *		4.634443E-2,-2.658234E-5, 7.757E-8,
     *		8.997041E-3, 6.329728E-6,-1.939256E-9,
     *		2.284178E-2,-9.941590E-5, 6.787400E-8,
     *		4.350267E-2,-6.839749E-5,-2.714956E-7,
     *		1.348204E-2, 1.091504E-5, 6.903760E-7,
     *		3.106570E-2,-1.665665E-4,-1.590188E-7/

	DATA DCARGS/5.0974222D0,-7.8604195454652D2,
     *3.9584962D0,-5.7533848094674D2,
     *1.6338070D0,-1.1506769618935D3,
     *2.5487111D0,-3.9302097727326D2,
     *4.9255514D0,-5.8849265665348D2,
     *1.3363463D0,-5.5076098609303D2,
     *1.6072053D0,-5.2237501616674D2,
     *1.3629480D0,-1.1790629318198D3,
     *5.5657014D0,-1.0977134971135D3,
     *5.0708205D0,-1.5774000881978D2,
     *3.9318944D0, 5.296346478D1,
     *4.8989497D0, 3.9809289073258D1,
     *1.3097446D0, 7.7540959633708D1,
     *3.5147141D0, 7.9618578146517D1,
     *3.5413158D0,-5.4868336758022D2/

	DATA CCAMPS/
     *-2.279594E-5, 1.407414E-5, 8.273188E-6, 1.340565E-5,-2.490817E-7,
     *-3.494537E-5, 2.860401E-7, 1.289448E-7, 1.627237E-5,-1.823138E-7,
     * 6.593466E-7, 1.322572E-5, 9.258695E-6,-4.674248E-7,-3.646275E-7,
     * 1.140767E-5,-2.049792E-5,-4.747930E-6,-2.638763E-6,-1.245408E-7,
     * 9.516893E-6,-2.748894E-6,-1.319381E-6,-4.549908E-6,-1.864821E-7,
     * 7.310990E-6,-1.924710E-6,-8.772849E-7,-3.334143E-6,-1.745256E-7,
     *-2.603449E-6, 7.359472E-6, 3.168357E-6, 1.119056E-6,-1.655307E-7,
     *-3.228859E-6, 1.308997E-7, 1.013137E-7, 2.403899E-6,-3.736225E-7,
     * 3.442177E-7, 2.671323E-6, 1.832858E-6,-2.394688E-7,-3.478444E-7,
     * 8.702406E-6,-8.421214E-6,-1.372341E-6,-1.455234E-6,-4.998479E-8,
     *-1.488378E-6,-1.251789E-5, 5.226868E-7,-2.049301E-7,0.E0,
     *-8.043059E-6,-2.991300E-6, 1.473654E-7,-3.154542E-7,0.E0,
     * 3.699128E-6,-3.316126E-6, 2.901257E-7, 3.407826E-7,0.e0,
     * 2.550120E-6,-1.241123E-6, 9.901116E-8, 2.210482E-7,0.e0,
     *-6.351059E-7, 2.341650E-6, 1.061492E-6, 2.878231E-7,0.e0/

	DATA CCSEC3/-7.757020E-8/
	DATA CCSEC/ 1.289600E-6,5.550147E-1,2.076942E0,
     *3.102810E-5,4.035027E0, 3.525565E-1,
     *9.124190E-6,9.990265E-1,2.622706E0,
     *9.793240E-7,5.508259E0, 1.559103E1/
     
	DATA DCSLD/1.990987D-7/, CCSGD/1.990969E-7/
	DATA CCKM/3.122140E-5/,CCMLD/2.661699E-6/,CCFDI/2.399485E-7/
	DATA DCARGM/5.1679830D0,8.3286911095275D3,
     *5.4913150D0,-7.2140632838100D3,
     *5.9598530D0, 1.5542754389685D4/
	DATA CCAMPM/
     * 1.097594E-1,2.896773E-7,5.450474E-2, 1.438491E-7,
     *-2.223581E-2,5.083103E-8,1.002548E-2,-2.291823E-8,
     * 1.148966E-2,5.658888E-8,8.249439E-3, 4.063015E-8/
	DATA CCPAMV/8.326827E-11,1.843484E-11,1.988712E-12,1.881276E-12/,
     *     DC1MME/.99999696D0/
     *     ,DAUKM/1.4959787D8/ !dolozyl KB (2003.03.28)

C  EXECUTION

	IDEQ=DEQ
	DT=(DJE-DCT0)/DCJUL
	T=DT
	DTSQ=DT*DT
	TSQ=DTSQ
	DO 100 K=1,8
	DLOCAL=DMOD(DCFEL(1,K)+DT*DCFEL(2,K)+DTSQ*DCFEL(3,K),DC2PI)
	IF(K.EQ.1) DML=DLOCAL
100   IF(K.NE.1) FORBEL(K-1)=DLOCAL
	DEPS=DMOD(DCEPS(1)+DT*DCEPS(2)+DTSQ*DCEPS(3),DC2PI)
	DO 200 K=1,17
200	SORBEL(K)=AMOD(CCSEL(1,K)+T*CCSEL(2,K)+TSQ*CCSEL(3,K),CC2PI)
	DO 300 K=1,4
300	SN(K)=SIN(AMOD(CCSEC(2,K)+T*CCSEC(3,K),CC2PI))
	PERTL=   CCSEC(1,1)*SN(1)+CCSEC(1,2)*SN(2)+
     +(CCSEC(1,3)+T*CCSEC3)*SN(3)+CCSEC(1,4)*SN(4)
	PERTLD=0.
	PERTR =0.
	PERTRD=0.
	DO 400 K=1,15
	A=DMOD(DCARGS(1,K)+DT*DCARGS(2,K),DC2PI)
	COSA=COS(A)
	SINA=SIN(A)
	PERTL=PERTL+CCAMPS(1,K)*COSA+CCAMPS(2,K)*SINA
	PERTR=PERTR+CCAMPS(3,K)*COSA+CCAMPS(4,K)*SINA
	IF(K.GE.11) GO TO 400
	PERTLD=PERTLD+(CCAMPS(2,K)*COSA-CCAMPS(1,K)*SINA)*CCAMPS(5,K)
	PERTRD=PERTRD+(CCAMPS(4,K)*COSA-CCAMPS(3,K)*SINA)*CCAMPS(5,K)
400	CONTINUE
	ESQ=E*E
	DPARAM=DC1-ESQ
	PARAM=DPARAM
	TWOE=E+E
	TWOG=G+G
	PHI=TWOE*((1.-ESQ*.125)*SIN(G)+E*.625*SIN(TWOG)+
     +     ESQ*.5416667*SIN(G+TWOG))
	F=G+PHI
	SINF=SIN(F)
	COSF=COS(F)
	DPSI=DPARAM/(DC1+E*COSF)
	PHID=E*CCSGD*((2.+ESQ*3.)*COSF+E*(2.5-SINF*SINF))
	PSID=E*CCSGD*SINF/SQRT(PARAM)

	D1PDRO=DC1+PERTR
	DRD=D1PDRO*(PSID+DPSI*PERTRD)
	DRLD=D1PDRO*DPSI*(DCSLD+PHID+PERTLD)
	DTL=DMOD(DML+PHI+PERTL,DC2PI)
	DSINLS=DSIN(DTL)
	DCOSLS=DCOS(DTL)
	DXHD=DRD*DCOSLS-DRLD*DSINLS
	DYHD=DRD*DSINLS+DRLD*DCOSLS

	PERTL=0.
	PERTLD=0.
	PERTP=0.
	PERTPD=0.
	DO 500 K=1,3
	A=DMOD(DCARGM(1,K)+DT*DCARGM(2,K),DC2PI)
	SINA=SIN(A)
	COSA=COS(A)
	PERTL=PERTL+CCAMPM(1,K)*SINA
	PERTLD=PERTLD+CCAMPM(2,K)*COSA
	PERTP=PERTP+CCAMPM(3,K)*COSA
500	PERTPD=PERTPD-CCAMPM(4,K)*SINA

	TL=FORBEL(2)+PERTL
	SINLM=SIN(TL)
	COSLM=COS(TL)
	SIGMA=CCKM/(1.+PERTP)
	A=SIGMA*(CCMLD+PERTLD)
	B=SIGMA*PERTPD
	DXHD=DXHD+A*SINLM+B*COSLM
	DYHD=DYHD-A*COSLM+B*SINLM
	DZHD=-SIGMA*CCFDI*COS(FORBEL(3))

	DXBD=DXHD*DC1MME
	DYBD=DYHD*DC1MME
	DZBD=DZHD*DC1MME
	DO 600 K=1,4
	PLON=FORBEL(K+3)
	POMG=SORBEL(K+1)
	PECC=SORBEL(K+9)
	TL=AMOD(PLON+2.*PECC*SIN(PLON-POMG),CC2PI)
	SINLP(K)=SIN(TL)
	COSLP(K)=COS(TL)
	DXBD=DXBD+CCPAMV(K)*(SINLP(K)+PECC*SIN(POMG))
	DYBD=DYBD-CCPAMV(K)*(COSLP(K)+PECC*COS(POMG))
600	DZBD=DZBD-CCPAMV(K)*SORBEL(K+13)*COS(PLON-SORBEL(K+5))

	DCOSEP=DCOS(DEPS)
	DSINEP=DSIN(DEPS)
	DYAHD=DCOSEP*DYHD-DSINEP*DZHD
	DZAHD=DSINEP*DYHD+DCOSEP*DZHD
	DYABD=DCOSEP*DYBD-DSINEP*DZBD
	DZABD=DSINEP*DYBD+DCOSEP*DZBD
C Ponizej czynnik DAUKM dolozony przez KB
	IF(IDEQ.NE.0) GO TO 700
	DVELH(1)=DXHD
	DVELH(2)=DYAHD
	DVELH(3)=DZAHD
	DVELB(1)=DXBD*DAUKM
	DVELB(2)=DYABD*DAUKM
	DVELB(3)=DZABD*DAUKM
	RETURN

700	DEQDAT=(DJE-DCT0-DCBES)/DCTROP+DC1900
	CALL PREc(DEQDAT,DEQ,DPREMA)
	DO 710 N=1,3
	DVELH(N)=DXHD*DPREMA(N,1)+DYAHD*DPREMA(N,2)+DZAHD*DPREMA(N,3)
710	DVELB(N)=
     *(DXBD*DPREMA(N,1)+DYABD*DPREMA(N,2)+DZABD*DPREMA(N,3))*DAUKM
	RETURN
	END




	SUBROUTINE PREc(DEQ1,DEQ2,DPREMA)

C CALCULATES THE MATRIX OF GENERAL PRECESSION FROM DEQ1 TO DEQ2

	IMPLICIT REAL*8 (D)
	DIMENSION DPREMA(3,3)
	DATA DCSAR/4.848136812D-6/,DC1900/1900.D0/,DC1M2/.01D0/,
     * DC1/2304.25D0/,DC2/1.396D0/,DC3/.302D0/,DC4/.018D0/,DC5/.791D0/
     * ,DC6/2004.683D0/,DC7/-.853D0/,DC8/-.426D0/,DC9/-.042D0/
	DT0=(DEQ1-DC1900)*DC1M2
	DT=(DEQ2-DEQ1)*DC1M2
	DTS=DT*DT
	DTC=DTS*DT
	DZETA=((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT=DZETA+DC5*DTS*DCSAR
	DTHET=((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSZETA=DSIN(DZETA)
	DCZETA=DCOS(DZETA)
	DSZETT=DSIN(DZETT)
	DCZETT=DCOS(DZETT)
	DSTHET=DSIN(DTHET)
	DCTHET=DCOS(DTHET)
	DA=DSZETA*DSZETT
	DB=DCZETA*DSZETT
	DC=DSZETA*DCZETT
	DD=DCZETA*DCZETT
	DPREMA(1,1)= DD*DCTHET-DA
	DPREMA(1,2)=-DC*DCTHET-DB
	DPREMA(1,3)=-DSTHET*DCZETT
	DPREMA(2,1)= DB*DCTHET+DC
	DPREMA(2,2)=-DA*DCTHET+DD
	DPREMA(2,3)=-DSTHET*DSZETT
	DPREMA(3,1)= DCZETA*DSTHET
	DPREMA(3,2)=-DSZETA*DSTHET
	DPREMA(3,3)= DCTHET
	RETURN
	END

	FUNCTION JD(L,M,N,J1G0)

c Input:  L - calendar year (years BC numbered 0, -1, -2, ...)
c         M - calendar month (for January M=1, February M=2, ..., M=12)
c         N - calendar day of the month M (1 to 28/29/30/31)
c      J1G0 - to be set to 1 for Julian and to 0 for Gregorian calendar
c Output: JD - Julian Day number
c  Calculates the Julian Day number (JD) from Gregorian or Julian
c  calendar dates. This integer number corresponds to the noon of 
c  the date (i.e. 12 hours of Universal Time).
c  The procedure was tested to be good since 1 March, -100100 (of both 
c  the calendars) up to a few millions (10**6) years into the future.
c  The algorithm is based on D.A. Hatcher, Q.Jl.R.Astron.Soc. 25(1984), 53-55
c  slightly modified by K.M. Borkowski (Post.Astron. 25(1987), 275-279).

      JD=(L+(M-8)/6+100100)*1461/4+(153*MOD(M+9,12)+2)/5+N-34840408
      IF(J1G0.LE.0)          JD = JD-(L+100100+(M-8)/6)/100*3/4+752
c     MJD=JD-2400000.5d0   ! this formula gives Modified Julian Day number
      END
	

	function dtuts(y)
	dimension idt(32),dt(235)		! dt(235=2014-1780+1)
c wg: FR Stephenson, LV Morrison, 1984, Phil Trans R Soc Lond A313, 47
c iy 1630 - 1775 (1625 dodalem wg rownania; 1780 z nast.)
	Data idt/78,    85,72,62,54,48,43,37,32,26,21,16,12,10,9,9,9,
     * 10,10,11,11,11,12,12,13,13,14,15,16,16,17,    17/
c iy 1780 - 1859
	Data dt/2*16.9,17.,6*17.1,17.,16.9,16.7,16.5,16.2,15.9,15.6,15.2
     *,14.8,14.4,14.1,13.7,13.4,13.1,12.9,12.7,12.6,11*12.5,12.4,12.3,
     *12.2,12.,11.7,11.4,11.1,10.6,10.2,9.6,9.1,8.6,8.,7.5,7.,6.6,6.3,6.
     *,5.8,5.7,3*5.6,5.7,5.8,5.9,6.1,6.2,6.3,6.5,6.6,6.8,6.9,7.1,7.2,
     *7.3,7.4,7.5,7.6,2*7.7,2*7.8,
c iy 1860 - 1980
     *7.88,7.82,7.54,6.97,6.4,6.02,5.41,4.1,2.92,1.82,1.61,
     *.1,-1.02,-1.28,-2.69,-3.24,-3.64,-4.54,-4.71,-5.11,-5.4,-5.42,
     *-5.2,2*-5.46,-5.79,-5.63,-5.64,-5.8,-5.66,-5.87,-6.01,-6.19,
     *-6.64,-6.44,-6.47,-6.09,-5.76,-4.66,-3.74,-2.72,-1.54,-.02,1.24,
     *2.64,3.86,5.37,6.14,7.75,9.13,10.46,11.53,13.36,14.65,16.01,17.2,
     *18.24,19.06,20.25,20.95,21.16,22.25,22.41,23.03,23.49,23.62,
     *23.86,24.49,24.34,24.08,24.02,24.,23.87,23.95,23.86,23.93,23.73,
     *23.92,23.96,24.02,24.33,24.83,25.3,25.7,26.24,26.77,27.28,27.78,
     *28.25,28.71,29.15,29.57,29.97,30.36,30.72,31.07,31.35,31.68,
     *32.18,32.68,33.15,33.59,34.,34.47,35.03,35.73,36.54,37.43,38.29,
     *39.2,40.18,41.17,42.22,43.37,44.48,45.47,46.46,47.52,48.53,49.59,
     *50.54,
c iy 1981 - 2014 (uaktualniono 2010.08.19 wg. Astron.Almanac 2011)
c iy:  1981    82    83     84    85   86    87    88    89    90   91
     *51.38,52.17,52.96,53.79,54.34,54.87,55.32,55.82,56.30,56.86,57.57
cc     1992   93    94    95    96     97    98    99   2000   1
     *,58.31,59.12,59.98,60.78,61.63,62.29,62.97,63.47,63.83,64.09,
c                 (AA data extrapolated to 2010 - 2014)
cc    2002   3     4     5     6    7     8     9    10   
     *64.3,64.47,65.57,64.69,64.85,65.15,65.46,65.78,66.1,
cc    2011 12  13 14
     *66.4,67.,67,67./

c Uaktualnienie tabeli dt: dopisanie N nowych wartosci powiazac ze
c zwiekszeniem wymiaru tabeli dt o N oraz w nastepnym wierszu
c powiekszyc graniczny rok po "if(iy.gt." takze o N.

	iytab=2014	!!! -> wymiar tab. dt(ii): ii=iytab+1-1780

	iy=y
	f=y-iy
	if(iy.gt.iytab-1.or.iy.lt.1637) go to 2

	if(iy.lt.1780) go to 1
	dtuts=dt(iy-1779)*(1-f)+dt(iy-1778)*f
	return

1	i=(iy-1620)/5
cc	dtuts=idt(i)+(iy-i*5-1619.5)*(idt(i+1)-idt(i))/5.
	dtuts=idt(i)+(y-i*5-1620)*(idt(i+1)-idt(i))/5.
	return
2	t=(y-1800)/100.
c w zasadzie wstecz tylko do ok. 700 BC
	if(iy.lt.948) dtuts=(44.3*t+320)*t+1360.
	if(iy.ge.948) dtuts=25.5*t*t 	! 44 s dla ciaglosci; nizej

c musze odjac wartosc skoku dla ciaglosci z danymi obserw. w XXI w.
		if(iy.ge.iytab) then
	dtutii=dt(iytab+1-1780)		! ostatnia wartosc w tabeli
	dtutform=25.5*((iytab-1800)/100.)**2 ! odpwiada jej ta wg wzoru
	dtuts= dtuts-(dtutform-dtutii)
		endif
	end

	SUBROUTINE PRE(DJE1,DJE2,dzeta,dzet,dthet)

C CALCULATES GENERAL PRECESSION FROM DJE1 TO DJE2

	IMPLICIT REAL*8 (D)
c Stale wg J.H.Lieske 1979, Astron. Astrophys. 73, 282.
	DATA DCSAR/4.848136812D-6/,J2000/2451545/,DC1/2306.2181D0/,
     * DC2/1.39656D0/,DC3/.30188D0/,DC4/.017998D0/,DC5/1.09468D0/
     * ,DC6/2004.3109D0/,DC7/-.8533D0/,DC8/-.42665D0/,DC9/-.041833D0/
	DT0=(Dje1-j2000)/36525
	DT=(Dje2-Dje1)/36525
	DTS=DT*DT
	DTC=DTS*DT
	DZETA=((DC1+(DC2-1.39d-4*dt0)*DT0)*DT+(DC3-3.44d-4*dt0)*DTS+
     * DC4*DTC)*DCSAR
	dzet =DZETA+((DC5-dc3+4.1d-4*dt0)*DTS+.000205*dtc)*DCSAR
	DTHET=((DC6+(DC7-2.17d-4*dt0)*DT0)*DT+(DC8-2.17d-4*dt0)*DTS+
     * DC9*DTC)*DCSAR
	END

