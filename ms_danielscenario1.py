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
	Model = random.randint(1,3)

#variables
	Thetamit = random.uniform(1,40.0)
	Thetanuc = 4*Thetamit
	DivTime = random.uniform(0.008, 0.36)
	DivTime2 = random.uniform(0, DivTime)
	NBot = random.uniform(5, 10)
		

#ms's commands
	if Model == 1:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> ./scenario1_ND2.txt" % (Thetamit, DivTime, DivTime))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> ./scenario1_16S.txt" % (Thetamit, DivTime, DivTime))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> ./scenario1_SiaH.txt" % (Thetanuc, DivTime, DivTime))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -ej %f 3 1 -ej %f 2 1 | perl msSS.pl >> ./scenario1_Rhodo.txt" % (Thetanuc, DivTime, DivTime))
		Counter += 1
	elif Model == 2:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario1_ND2.txt" % (Thetamit, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario1_16S.txt" % (Thetamit, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario1_SiaH.txt" % (Thetanuc, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -en %f 3 %f -ej %f 3 2 -ej %f 2 1 | perl msSS.pl >> ./scenario1_Rhodo.txt" % (Thetanuc, DivTime2, NBot, DivTime2, DivTime))
		Counter += 1
	elif Model == 3:
		os.system("./ms 168 1 -t %f -I 3 131 34 3 -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario1_ND2.txt" % (Thetamit, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 151 1 -t %f -I 3 129 18 4 -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario1_16S.txt" % (Thetamit, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 310 1 -t %f -I 3 260 46 4 -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario1_SiaH.txt" % (Thetanuc, DivTime2, NBot, DivTime2, DivTime))
		os.system("./ms 378 1 -t %f -I 3 292 82 4 -en %f 3 %f -ej %f 3 1 -ej %f 1 2 | perl msSS.pl >> ./scenario1_Rhodo.txt" % (Thetanuc, DivTime2, NBot, DivTime2, DivTime))
		Counter += 1
		
#print prior values
	print '_%d\t%f\t%f\t%f\t%f\t%f' % (Model, Thetamit, Thetanuc, DivTime, DivTime2, NBot)
