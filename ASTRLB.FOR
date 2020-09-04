C
C Library for software astrolabes
C
C Copyright (c) 2020, Michael Robinson
C
	BLOCK DATA ASTRLB
	INTEGER MONTHS(12)
	COMMON MONTHS
	DATA MONTHS/31,28,31,30,31,30,31,31,30,31,30,31/
	END
C
	FUNCTION RADC2E(LAT,RA,DEC,ST)
C Convert latitude (deg),right ascension (hr), declination (deg), 
C and local sideral time (rad) to elevation (radians)
	REAL LAT,RA,DEC,ST,RADC2E
	REAL LATRD,RARD,DECRD,LHA,NUM,DEN
C
	LATRD=LAT*3.14159/180.0
	RARD=RA*6.28319/24
	DECRD=DEC*3.14159/180.0
C       
	LHA=ST-RARD
C Bit of a hack below since PDP-11 fortran does have ASIN()
	NUM=SIN(LATRD)*SIN(DECRD)+COS(LATRD)*COS(DECRD)*COS(LHA)
	DEN=SQRT(1-NUM**2)
	RADC2E=ATAN2(NUM,DEN)
	END
C
	FUNCTION RADC2A(LAT,RA,DEC,ST)
C Convert latitude (deg),right ascension (hr), declination (deg), 
C local sideral time (rad) to azimuth (radians)
	REAL LAT,RA,DEC,ST,RADC2A,EL,NUM,DEN
	REAL LATRD,RARD,DECRD,LHA
C
	LATRD=LAT*3.14159/180.0
	RARD=RA*6.28319/24
	DECRD=DEC*3.14159/180.0
C
	EL=RADC2E(LAT,RA,DEC,ST)
	LHA=ST-RARD
C	
	NUM=-SIN(LHA)*COS(DECRD)/COS(EL)
	DEN=(SIN(DECRD)-SIN(EL)*SIN(LATRD))/(COS(EL)*COS(LATRD))
	RADC2A=ATAN2(NUM,DEN)
C
	IF (RADC2A.LT.0) RADC2A=RADC2A+6.28319
	END
C
	FUNCTION UTC2S(LON,D,H,M)
C Convert Longitude (degrees), Julian date past 1 Jan 2000,
C UTC (hours, minutes) to local sidereal time in radians
C Source https://en.wikipedia.org/wiki/Sidereal_time
	REAL LON,UTC2S
	INTEGER D,H,M
	UTC2S=D+(H-12)/24.0+M/1440.0
	UTC2S=0.77905727+1.0027378*UTC2S
	UTC2S=UTC2S+LON/360.0
	UTC2S=(UTC2S-INT(UTC2S))*6.28319
	END
C                                     
	FUNCTION JDATE(YEAR,MONTH,DAY)
C Convert month and day to Julian date past 1 Jan 2000
	INTEGER YEAR,MONTH,DAY
	INTEGER I,Y,JDATE,MONTHS(12)
	COMMON MONTHS	
	JDATE=DAY
C Accumulate days for all years up to this one
C Account for leap years as well
	DO 310 Y=2000,YEAR-1
	JDATE=JDATE+365
	IF (MOD(Y,4).EQ.0) JDATE=JDATE+1
	IF (MOD(Y,100).EQ.0) JDATE=JDATE-1
	IF (MOD(Y,400).EQ.0) JDATE=JDATE+1
310	CONTINUE
C Accumulate days from the current year
C Leap year accounting only happens in month 2 (February)
	DO 300 I=1,MONTH-1
	JDATE=JDATE+MONTHS(I)
	IF (I.NE.2) GOTO 300
	IF (MOD(YEAR,4).EQ.0) JDATE=JDATE+1
	IF (MOD(YEAR,100).EQ.0) JDATE=JDATE-1
	IF (MOD(YEAR,400).EQ.0) JDATE=JDATE+1
300	CONTINUE
	END
C
	FUNCTION UTC2JD(YEAR,MONTH,DAY,HOUR,MINUTE)
C Compute Julian date with UTC hour and minute
	INTEGER YEAR,MONTH,DAY,HOUR,MINUTE
	INTEGER JD
	REAL UTC2JD
	JD=JDATE(YEAR,MONTH,DAY)
	UTC2JD=JD+HOUR/24.0+MINUTE/1440.0
	END
C
	FUNCTION SUNRA(N)
C Compute right ascension (in hours) of the sun from its day past 1 Jan 2000
C Source https://en.wikipedia.org/wiki/Position_of_the_Sun
	REAL L,G,LAM,PI,EPS,SUNRA
	INTEGER N
	PI=3.14159
	L=280.460+0.9856003*N
	G=357.528+0.9856003*N
	LAM=L+1.915*SIN(G*PI/180)+0.02*SIN(2*PI/180*G)
	EPS=23.439+0.00000004*N
C       
	SUNRA=ATAN2(COS(EPS*PI/180)*SIN(LAM*PI/180),COS(LAM*PI/180))
	SUNRA=SUNRA*24/6.28319
	END
C
	FUNCTION SUNDC(N)
C Compute declination of the sun (in degrees) of the sun from its day
C past 1 Jan 2000
C Source https://en.wikipedia.org/wiki/Position_of_the_Sun
	REAL L,G,LAM,PI,EPS,SUNDC,NUM,DEN
	INTEGER N
	PI=3.14159
	L=280.460+0.9856003*N
	G=357.528+0.9856003*N
	LAM=L+1.915*SIN(G*PI/180)+0.02*SIN(2*PI/180*G)
	EPS=23.439+0.00000004*N
C Hack below because PDP-11 fortran does not have ASIN()
	NUM=SIN(EPS*PI/180)*SIN(LAM*PI/180)
	DEN=SQRT(1-NUM**2)
	SUNDC=ATAN2(NUM,DEN)
	SUNDC=SUNDC*180/PI
	END
C
	FUNCTION EQNTIM(D)
C Compute equation of time (minutes) for Julian day D
C Source https://en.wikipedia.org/wiki/Equation_of_time
	REAL M,EQNTIM
	INTEGER D
	M=6.24004077+0.01720197*D
	EQNTIM=-7.659*SIN(M)+9.863*SIN(2*M+3.5932)
	END
C
	SUBROUTINE MOONPOS(D,MOONRA,MOONDC)
C Compute the right ascention (hours) and declination (degrees)
C of the moon for Julian day D
	REAL MOONRA,MOONDC
	REAL D
	REAL NN,II,WW,AA,EE,MM,PI
	REAL E,XV,YV,V,R
	REAL XG,YG,ZG
	REAL ECL
	REAL XE,YE,ZE
	PI=3.14159
C Oribital elements of the moon
C NN = longitude of ascending node (radians)
C II = inclination to the ecliptic (radians)
C WW = Argument of perigee (radians)
C AA = mean distance to earth (Earth radii)
C EE = Eccentricity (dimensionless)
C MM = Mean anomaly (radians)
	NN=(125.1228-0.0529538083*D)*PI/180
	II=5.1454*PI/180
	WW=(318.0634+0.1643573223*D)*PI/180
	AA=60.2666
	EE=0.054900
	MM=(115.3654+13.0649929509*D)*PI/180
C Moon's true distance R (Earth radii) and 
C true anomaly (radians)
	E=MM+EE*SIN(MM)*(1.0+EE*COS(MM))
        XV=AA*(COS(E)-EE)
	YV=AA*SQRT(1.0-EE*EE)*SIN(E)
	V=ATAN2(YV,XV)
	R=SQRT(XV*XV+YV*YV)
C Moon's position in geocentric rectangular coordinates
	XG=R*(COS(NN)*COS(V+WW)-SIN(NN)*SIN(V+WW)*COS(II))
	YG=R*(SIN(NN)*COS(V+WW)+COS(NN)*SIN(V+WW)*COS(II))
        ZG=R*SIN(V+WW)*SIN(II)
C Ecliptic coordinates
	ECL=(23.4393-3.563E-7*D)*PI/180
	XE=XG
	YE=YG*COS(ECL)-ZG*SIN(ECL)
	ZE=YG*SIN(ECL)+ZG*COS(ECL)
C Final output in equatorial coordinates
	MOONRA=ATAN2(YE,XE)*24/6.28319
	MOONDC=ATAN2(ZE,SQRT(XE*XE+YE*YE))*180/PI
	END
C
	SUBROUTINE JUPPOS(D,JUPRA,JUPDC)
C Compute the right ascension (hours) and declination (degrees)
C of Jupiter for Julian day D
	REAL JUPRA,JUPDC
	REAL D
	REAL NN,II,WW,AA,EE,MM,PI
	REAL E,XV,YV,R
	REAL XH,YH,ZH
	REAL XG,YG,ZG
        REAL LATECL,LONECL
	REAL EES,MMS,WWS
	REAL ES,XVS,YVS,VS,RS,LONSUN,XS,YS
	REAL ECL
	REAL XE,YE,ZE
	PI=3.14159
C Orbital elements of Jupiter (about the sun)
C NN = longitude of ascending node (radians)
C II = inclination to the ecliptic (radians)
C WW = Argument of perihelion (radians)
C AA = mean distance to the sun (AU)
C EE = eccentricity (dimensionless)
C MM = Mean anomaly (radians)
	NN=(100.4542+2.76854E-5*D)*PI/180
	II=(1.3030-1.557E-7*D)*PI/180
	WW=(273.8777+1.64505E-5*D)*PI/180
	AA=5.20256
	EE=0.048498+4.469E-9*D
	MM=(19.8950+0.0830853001*D)*PI/180
C Jupiter's distance R (AU) from the sun
C true anomaly V (radians)
	E=MM+EE*SIN(MM)*(1.0+EE*COS(MM))
	XV=AA*(COS(E)-EE)
	YV=AA*SQRT(1.0-EE*EE)*SIN(E)
	V=ATAN2(YV,XV)
	R=SQRT(XV*XV+YV*YV)
C Jupiter's position in space (centered on the sun
C units in AU)
	XH=R*(COS(NN)*COS(V+WW)-SIN(NN)*SIN(V+WW)*COS(II))
	YH=R*(SIN(NN)*COS(V+WW)+COS(NN)*SIN(V+WW)*COS(II))
	ZH=R*(SIN(V+WW)*SIN(II))
C Ecliptic coordinates (radians)
	LONECL=ATAN2(YH,XH)
	LATECL=ATAN2(ZH,SQRT(XH*XH+YH*YH))
C Orbital elements of the sun
	EES=(0.016709-1.151E-9*D)*PI/180
	MMS=(356.0470+0.9856002585*D)*PI/180
	WWS=(282.9404+4.70935E-5*D)*PI/180
C Position of the sun
	ES=MMS+EES*SIN(MMS)*(1.0+EES*COS(MMS))
	XVS=COS(ES)-EES
	YVS=SQRT(1.0-EES*EES)*SIN(ES)
	VS=ATAN2(YVS,XVS)	
	RS=SQRT(XVS*XVS+YVS*YVS)
	LONSUN=VS+WWS
	XS=RS*COS(LONSUN)
	YS=RS*SIN(LONSUN)
C Geocentric coodinates of Jupiter
	XG=XH+XS
	YG=YH+YS
	ZG=ZH
C Ecliptic coordinates
	ECL=(23.4393-3.563E-7*D)*PI/180
	XE=XG
	YE=YG*COS(ECL)-ZG*SIN(ECL)
	ZE=YG*SIN(ECL)+ZG*COS(ECL)
C Final output in equatorial coordinates
	JUPRA=ATAN2(YE,XE)*24/6.28319
	JUPDC=ATAN2(ZE,SQRT(XE*XE+YE*YE))*180/PI
	END
