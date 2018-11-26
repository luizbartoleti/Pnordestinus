#!/usr/bin/python

#import modules
import random
import os
import math

#draw from prior
SimNum = 300000
Counter = 0


#while loop
while Counter < SimNum:
	Model = random.randint(4,6)

#variables
	Thetamit = random.uniform(1,40.0)
	Thetanuc = 4*Thetamit
	DivNS = random.uniform(0, 0.008)
	DivWS = random.uniform(0, 0.008)
	GrowthRatio = random.uniform(0.25, 0.9)
	GrowthRateN = -(1/DivNS) * math.log(GrowthRatio)
	GrowthRateW = -(1/DivWS) * math.log(GrowthRatio)
	GrowthRateS = -(1/0.008) * math.log(GrowthRatio)
	DivWS2 = random.uniform(0, 0.0014)
	GrowthRateN2= -(1/0.0014) * math.log(GrowthRatio)
	GrowthRateW2 = -(1/DivWS2) * math.log(GrowthRatio)
	GrowthRateS2 = -(1/0.0014) * math.log(GrowthRatio)
	NBotN = random.uniform(5, 10)
	NBotW = random.uniform(5, 10) 
	DivNS2 = random.uniform(0.0014, 0.008)
	DivWN = random.uniform(0, 0.0014)
	GrowthRateW3 = -(1/DivWN) * math.log(GrowthRatio)

#ms's commands
	if Model == 4:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -eg 0.008 2 0.0 -ej %f 3 2 -ej %f 1 2 | perl msSS.pl >> ./scenario2_ND2.txt" % (Thetamit, GrowthRateN, GrowthRateS, GrowthRateW, DivNS, NBotN, DivWS, NBotW, DivWS, DivNS))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -eg 0.008 2 0.0 -ej %f 3 2 -ej %f 1 2 | perl msSS.pl >> ./scenario2_16S.txt" % (Thetamit, GrowthRateN, GrowthRateS, GrowthRateW, DivNS, NBotN, DivWS, NBotW, DivWS, DivNS))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -eg 0.008 2 0.0 -ej %f 3 2 -ej %f 1 2 | perl msSS.pl >> ./scenario2_SiaH.txt" % (Thetanuc, GrowthRateN, GrowthRateS, GrowthRateW, DivNS, NBotN, DivWS, NBotW, DivWS, DivNS))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -g 1 %f -g 2 %f -g 3 %f -en %f 1 %f -en %f 3 %f -eg 0.008 2 0.0 -ej %f 3 2 -ej %f 1 2 | perl msSS.pl >> ./scenario2_Rhodo.txt" % (Thetanuc, GrowthRateN, GrowthRateS, GrowthRateW, DivNS, NBotN, DivWS, NBotW, DivWS, DivNS))
		Counter += 1
	elif Model == 5:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario2_ND2.txt" % (Thetamit, GrowthRateN2, GrowthRateS2, GrowthRateW2, DivWS2, NBotW, DivWS2, DivNS2))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario2_16S.txt" % (Thetamit, GrowthRateN2, GrowthRateS2, GrowthRateW2, DivWS2, NBotW, DivWS2, DivNS2))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario2_SiaH.txt" % (Thetanuc, GrowthRateN2, GrowthRateS2, GrowthRateW2, DivWS2, NBotW, DivWS2, DivNS2))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario2_Rhodo.txt" % (Thetanuc, GrowthRateN2, GrowthRateS2, GrowthRateW2, DivWS2, NBotW, DivWS2, DivNS2))
		Counter += 1
	elif Model == 6:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario2_ND2.txt" % (Thetamit, GrowthRateN2, GrowthRateS2, GrowthRateW3, DivWN, NBotW, DivWN, DivNS2))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario2_16S.txt" % (Thetamit, GrowthRateN2, GrowthRateS2, GrowthRateW3, DivWN, NBotW, DivWN, DivNS2))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario2_SiaH.txt" % (Thetanuc, GrowthRateN2, GrowthRateS2, GrowthRateW3, DivWN, NBotW, DivWN, DivNS2))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -g 1 %f -g 2 %f -g 3 %f -en %f 3 %f -eg 0.0014 1 0.0 -eg 0.0014 2 0.0 -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario2_Rhodo.txt" % (Thetanuc, GrowthRateN2, GrowthRateS2, GrowthRateW3, DivWN, NBotW, DivWN, DivNS2))
		Counter += 1
		
#print prior values
	print '_%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f' % (Model, Thetamit, Thetanuc, DivNS, DivWS, GrowthRateN, GrowthRateW, GrowthRateS, DivWS2, GrowthRateN2, GrowthRateW2, GrowthRateS2, NBotN, NBotW, DivNS2, DivWN, GrowthRateW3)
