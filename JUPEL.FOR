	PROGRAM MOONEL
C Compute Jupiter elevation
C
C Copyright (c) 2020, Michael Robinson
C
	INTEGER YEAR,MONTH,DAY,HOUR,MINUTE,JD
	REAL JD2
	REAL LAT,LON,ST,SRA,SDC,SEL
C Default home location
	LAT=39
	LON=-77
C
	TYPE *, 'Enter Year Month Day Hour Minute (UTC) '
 	ACCEPT *, YEAR,MONTH,DAY,HOUR,MINUTE
C Compute Julian date
	JD=JDATE(YEAR,MONTH,DAY)
	JD2=UTC2JD(YEAR,MONTH,DAY,HOUR,MINUTE)
C Compute Sidereal time
	ST=UTC2S(LON,JD,HOUR,MINUTE)
C Compute Jupiter's position
	CALL JUPPOS(JD2,SRA,SDC)
	SEL=RADC2E(LAT,SRA,SDC,ST)
	TYPE 60, JD2
	TYPE 50, SRA,SDC
	TYPE 40, SEL*180/3.14159
60	FORMAT (' Julian date : ',F20.5)
40	FORMAT (' Jupiter elevation : ',F9.5,' deg')
50	FORMAT (' Jupiter RA,DEC : ',F9.5,' hr, ',F9.5,' deg')
20	END
