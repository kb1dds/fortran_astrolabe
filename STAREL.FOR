	PROGRAM STAREL
C Compute star elevations
C
C Copyright (c) 2020, Michael Robinson
C
	INTEGER YEAR,MONTH,DAY,HOUR,MINUTE,JD
	REAL LAT,LON,ST,SRA,SDC,SEL
	BYTE STAR(10)
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
C Open star chart
	OPEN(UNIT=1,NAME='STARS.DAT',TYPE='OLD')
10	READ (1,30,END=20) STAR,SRA,SDC
30	FORMAT (10A,' ',F9.5,F9.5)
	SEL=RADC2E(LAT,SRA,SDC,ST)
	TYPE 40,STAR,SEL*180/3.14159
40	FORMAT (' ',10A,' elevation : ',F9.5,' deg')
	GOTO 10
20	END
                                                                                                                                                                                                                                                                                                                                                              