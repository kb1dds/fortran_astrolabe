	PROGRAM SUNPOS
C Compute position of the sun
C
C Copyright (c) 2020, Michael Robinson
	INTEGER YEAR,MONTH,DAY,HOUR,MINUTE,JD
	REAL LAT,LON,ST,SRA,SDC,SEL
C Default home location
	LAT=39
	LON=-77
C
	TYPE *,'Enter Year Month Day Hour Minute (UTC) '
 	ACCEPT *,YEAR,MONTH,DAY,HOUR,MINUTE
C Compute Julian date
	JD=JDATE(YEAR,MONTH,DAY)
C Compute Sidereal time
	ST=UTC2S(LON,JD,HOUR,MINUTE)
C Locate the sun
	SRA=SUNRA(JD)
	SDC=SUNDC(JD)
C Elevation of the sun
	SEL=RADC2E(LAT,SRA,SDC,ST)*180/6.28319
	TYPE *,'Sun elevation: ',SEL,' degrees'
	END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                         