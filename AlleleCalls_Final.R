#########################################################Raw Data Preparation############################################################
#Dependencies#
require(MsatAllele) #Archived version of MSatAllele is available at:  https://cran.r-project.org/src/contrib/Archive/MsatAllele/

##########Load Data##########
###Note on data format###
#Data is loaded from a .csv file (preparable in Microsoft Excel).
#Dataframe should have a CODE column that indicates individual ID followed by two columns for each locus, one for each allele.  
#See included example file for formatting reference.  
#############################

###########Format Data#########
#The following code reformats the data for use in microsatellite binning procedure
gen.data<-read.csv("example_samples.csv")  
row<-as.vector(gen.data$CODE)
marker<-cbind("RT9","BM4107","P","Cervid1","Q","RT7", "BM6438", "INRA011",
              "RT5", "OarFCB193", "BL42")
n=as.numeric(length(gen.data[,1]))
mx1<-rep(1,n)
mx2<-rep(2,n)
mx3<-rep(2,n)
gen.data[is.na(gen.data)]<-0
#############################

###########Data Preparation############
#The following code creates an object for each locus in MSatAllele Format.

RT9<-rep("RT9",n)
RT9<-cbind(row,mx1,RT9,as.numeric(gen.data$RT9a1), as.numeric(gen.data$RT9a2))

BM4107<-rep("BM4107",n)
BM4107<-cbind(row,mx1,BM4107,gen.data$BM4107a1, gen.data$BM4107a2)

P<-rep("P",n)
P<-cbind(row,mx1,P,gen.data$Pa1, gen.data$Pa2)

Cervid1<-rep("Cervid1",n)
Cervid1<-cbind(row,mx1,Cervid1,gen.data$Cervid1a1, as.character(gen.data$Cervid1a2))

Q<-rep("Q",n)
Q<-cbind(row,mx2,Q,gen.data$Qa1, gen.data$Qa2)

RT5<-rep("RT5", n)
RT5<-cbind(row, mx2, RT5, gen.data$RT5a1, gen.data$RT5a2)

BL42<-rep("BL42", n)
BL42<-cbind(row, mx2, BL42, gen.data$BL42a1, gen.data$BL42a2)


RT7<-rep("RT7",n)
RT7<-cbind(row,mx3,RT7,as.character(gen.data$RT7a1), gen.data$RT7a2)

BM6438<-rep("BM6438",n)
BM6438<-cbind(row,mx3,BM6438,gen.data$BM6438a1, gen.data$BM6438a2)

INRA011<-rep("INRA011",n)
INRA011<-cbind(row,mx3,INRA011,gen.data$INRA011a1, gen.data$INRA011a2)

OarFCB193<-rep("OarFCB193", n)
OarFCB193<-cbind(row, mx3, OarFCB193, gen.data$OarFCB193a1, gen.data$OarFCB193a2)



###########Format Full Data Frame###########
data<-rbind(RT9, BM4107, P, Cervid1, Q, RT7, BM6438, INRA011, RT5, OarFCB193, BL42)
colnames(data)<-c("Sample Name","Panel","Marker","Size1","Size2")
data<-as.data.frame(data)
data<-data[order(data[,1],data[,2],data[,3]),]
data$Size1[data$Size1==0]<-NA
data$Size2[data$Size2==0]<-NA
data$Size1[data$Size1=="NA"]<-NA
data$Size2[data$Size2=="NA"]<-NA

write.table(data,"sample_markers.txt",sep="\t",na="",quote=F, col.names=F,row.names=F)
#############################




#######Allele Calls########
##Read MSatAllele Data for Cumulative Allele plots##
gen.data<-read.frag.sizes("sample_markers.txt", date="2018",plate="BMCSamples") #Make sure this matches table name
####

##RT9 Bin##
par(mfrow=c(1,2))
hist(as.numeric(RT9[,4:5]),main="RT9",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"RT9")

RT9.1<-RT9 #Save copy of data formatted using MSatAllele
RT9<-RT9[,4:5] #Reformat for final calls 

RT9<-ifelse(RT9>=101.5 & RT9<103.2,102,RT9)
RT9<-ifelse(RT9>=103.2 & RT9<105,104,RT9)
RT9<-ifelse(RT9>=105.0 & RT9<107.2,106,RT9)
RT9<-ifelse(RT9>=107.2 & RT9<109.1,108,RT9)
RT9<-ifelse(RT9>=109.1 & RT9<111.1,110,RT9)
RT9<-ifelse(RT9>=111.1 & RT9<112.5,112,RT9)
RT9<-ifelse(RT9>=112.5 & RT9<114.4,113,RT9)
RT9<-ifelse(RT9>=114.4 & RT9<116.2,115,RT9)
RT9<-ifelse(RT9>=116.2 & RT9<118.6,117,RT9)
RT9<-ifelse(RT9>=118.6 & RT9<119.9,119,RT9)
RT9<-ifelse(RT9>=119.9 & RT9<121.8,121,RT9)
RT9<-ifelse(RT9>=121.8 & RT9<123,123,RT9)
RT9<-ifelse(RT9>=123.0 & RT9<125.7,125,RT9)
####


##BM4107 Bin##
par(mfrow=c(1,2))
hist(as.numeric(BM4107[,4:5]),main="BM4107",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"BM4107")


BM4107.1<-BM4107
BM4107<-BM4107[,4:5]
BM4107<-ifelse(BM4107>=133.7 & BM4107<135.0,134,BM4107)
BM4107<-ifelse(BM4107>=135.0 & BM4107<137.2,136,BM4107)
BM4107<-ifelse(BM4107>=137.3 & BM4107<140.0,139,BM4107)
BM4107<-ifelse(BM4107>=140.0 & BM4107<142.0,141,BM4107)
BM4107<-ifelse(BM4107>=142.0 & BM4107<144.0,143,BM4107)
BM4107<-ifelse(BM4107>=144.0 & BM4107<146.0,145,BM4107)
BM4107<-ifelse(BM4107>=146.0 & BM4107<148.0,147,BM4107)
BM4107<-ifelse(BM4107>=148.0 & BM4107<150.0,149,BM4107)
BM4107<-ifelse(BM4107>=150.0 & BM4107<152.0,151,BM4107)
BM4107<-ifelse(BM4107>=152.0 & BM4107<154.0,153,BM4107)
BM4107<-ifelse(BM4107>=154.0 & BM4107<155.4,155,BM4107)
BM4107<-ifelse(BM4107>=155.5 & BM4107<157.0,156,BM4107)
BM4107<-ifelse(BM4107>=157.0 & BM4107<159.0,158,BM4107)
BM4107<-ifelse(BM4107>=159.0 & BM4107<161.0,160,BM4107)
BM4107<-ifelse(BM4107>=161.0 & BM4107<163.0,162,BM4107)
BM4107<-ifelse(BM4107>=163.0 & BM4107<165.0,164,BM4107)
BM4107<-ifelse(BM4107>=165.0 & BM4107<167.0,166,BM4107)
BM4107<-ifelse(BM4107>=167.0 & BM4107<169.0,168,BM4107)
####

##P Bin##
par(mfrow=c(1,2))
hist(as.numeric(P[,4:5]),main="P",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"P")

P.1<-P
P<-P[,4:5]
P<-ifelse(P>=210.0 & P<213.0,210,P)
P<-ifelse(P>=215.0 & P<218.0,216,P)
P<-ifelse(P>=218.0 & P<222.0,220,P)
P<-ifelse(P>=222.0 & P<226.0,224,P)
P<-ifelse(P>=227.0 & P<229.0,228,P)
P<-ifelse(P>=230.0 & P<233.5,232,P)
P<-ifelse(P>=234.0 & P<237.5,236,P)
P<-ifelse(P>=239.0 & P<241.6,240,P)
P<-ifelse(P>=243.0 & P<246.0,244,P)
P<-ifelse(P>=247.0 & P<250.0,248,P)
####



##Cervid1 Bin##
par(mfrow=c(1,2))
hist(as.numeric(Cervid1[,4:5]),main="Cervid1",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"Cervid1")


Cervid1.1<-Cervid1
Cervid1<-Cervid1[,4:5]

Cervid1<-ifelse(Cervid1>=158.5 & Cervid1<160.0,159,Cervid1)
Cervid1<-ifelse(Cervid1>=160.0 & Cervid1<162.0,161,Cervid1)
Cervid1<-ifelse(Cervid1>=162.0 & Cervid1<164.0,163,Cervid1)
Cervid1<-ifelse(Cervid1>=164.0 & Cervid1<166.0,165,Cervid1)
Cervid1<-ifelse(Cervid1>=166.0 & Cervid1<168.0,167,Cervid1)
Cervid1<-ifelse(Cervid1>=168.0 & Cervid1<170.0,169,Cervid1)
Cervid1<-ifelse(Cervid1>=170.0 & Cervid1<172.0,171,Cervid1)
Cervid1<-ifelse(Cervid1>=172.0 & Cervid1<173.9,173,Cervid1)
Cervid1<-ifelse(Cervid1>=173.9 & Cervid1<175.5,175,Cervid1)
Cervid1<-ifelse(Cervid1>=175.5 & Cervid1<177.1,176,Cervid1)
Cervid1<-ifelse(Cervid1>=177.1 & Cervid1<179.0,178,Cervid1)
Cervid1<-ifelse(Cervid1>=179.0 & Cervid1<181.2,180,Cervid1)
Cervid1<-ifelse(Cervid1>=181.2 & Cervid1<183.0,182,Cervid1)
Cervid1<-ifelse(Cervid1>=183.0 & Cervid1<185.0,184,Cervid1)
Cervid1<-ifelse(Cervid1>=185.0 & Cervid1<187.0,186,Cervid1)
Cervid1<-ifelse(Cervid1>=187.0 & Cervid1<189.0,188,Cervid1)
Cervid1<-ifelse(Cervid1>=189.0 & Cervid1<191.0,190,Cervid1)
Cervid1<-ifelse(Cervid1>=191.0 & Cervid1<193.0,192,Cervid1)
Cervid1<-ifelse(Cervid1>=193.0 & Cervid1<195.0,194,Cervid1)
Cervid1<-ifelse(Cervid1>=195.0 & Cervid1<197.0,196,Cervid1)
Cervid1<-ifelse(Cervid1>=197.0 & Cervid1<199.0,198,Cervid1)
####


##Q Bin##
par(mfrow=c(1,2))
hist(as.numeric(Q[,4:5]),main="Q",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"Q")

Q.1<-Q
Q<-Q[,4:5]
Q<-ifelse(Q>=226.0 & Q<230.0,228,Q)
Q<-ifelse(Q>=230.0 & Q<234.0,232,Q)
Q<-ifelse(Q>=234.0 & Q<238.0,236,Q)
Q<-ifelse(Q>=238.0 & Q<242.0,240,Q)
Q<-ifelse(Q>=242.0 & Q<246.2,244,Q)
Q<-ifelse(Q>=246.2 & Q<250.0,248,Q)
Q<-ifelse(Q>=250.0 & Q<254.1,252,Q)
Q<-ifelse(Q>=254.1 & Q<258.0,256,Q)
Q<-ifelse(Q>=258.0 & Q<262.0,260,Q)
Q<-ifelse(Q>=262.0 & Q<264.3,263,Q)
Q<-ifelse(Q>=264.3 & Q<266.1,265,Q)
Q<-ifelse(Q>=266.1 & Q<268.0,267,Q)
Q<-ifelse(Q>=268.0 & Q<270.0,269,Q)
Q<-ifelse(Q>=270.0 & Q<272.0,271,Q)
Q<-ifelse(Q>=272.0 & Q<274.0,273,Q)
Q<-ifelse(Q>=274.0 & Q<277.0,275,Q)
Q<-ifelse(Q>=277.0 & Q<280.2,279,Q)
Q<-ifelse(Q>=280.2 & Q<285.0,283,Q)
Q<-ifelse(Q>=285.0 & Q<289.0,287,Q)
Q<-ifelse(Q>=289.0 & Q<293.0,291,Q)
Q<-ifelse(Q>=293.0 & Q<297.0,295,Q)
####

##RT7 Bin##
par(mfrow=c(1,2))
hist(as.numeric(RT7[,4:5]),main="RT7",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"RT7")

RT7.1<-RT7
RT7<-RT7[,4:5]
RT7<-ifelse(RT7>=204.0 & RT7<206.0,205,RT7)
RT7<-ifelse(RT7>=206.0 & RT7<208.0,207,RT7)
RT7<-ifelse(RT7>=208.0 & RT7<210.0,209,RT7)
RT7<-ifelse(RT7>=210.0 & RT7<212.4,211,RT7)
RT7<-ifelse(RT7>=212.4 & RT7<214.5,213,RT7)
RT7<-ifelse(RT7>=214.5 & RT7<216.4,215,RT7)
RT7<-ifelse(RT7>=216.4 & RT7<218.3,217,RT7)
RT7<-ifelse(RT7>=218.3 & RT7<220.4,219,RT7)
RT7<-ifelse(RT7>=220.4 & RT7<222.3,221,RT7)
RT7<-ifelse(RT7>=222.3 & RT7<224.5,223,RT7)
RT7<-ifelse(RT7>=224.5 & RT7<226.5,225,RT7)
RT7<-ifelse(RT7>=226.5 & RT7<228.5,227,RT7)
RT7<-ifelse(RT7>=228.5 & RT7<230.5,229,RT7)
RT7<-ifelse(RT7>=230.5 & RT7<232.5,231,RT7)
RT7<-ifelse(RT7>=232.5 & RT7<233.8,233,RT7)
RT7<-ifelse(RT7>=233.8 & RT7<236.0,235,RT7)
RT7<-ifelse(RT7>=236.0 & RT7<238.4,237,RT7)
RT7<-ifelse(RT7>=238.5 & RT7<240.0,239,RT7)
RT7<-ifelse(RT7>=240.0 & RT7<242.0,241,RT7)
RT7<-ifelse(RT7>=242.0 & RT7<244.0,243,RT7)
####

##BM6438 Bin##
par(mfrow=c(1,2))
hist(as.numeric(BM6438[,4:5]),main="BM6438",xlab="Raw Allele",breaks=10000)    
AlleleCum(gen.data,"BM6438")

BM6438.1<-BM6438
BM6438<-BM6438[,4:5]
BM6438<-ifelse(BM6438>=249.0 & BM6438<252.0,251,BM6438)
BM6438<-ifelse(BM6438>=252.0 & BM6438<253.9,253,BM6438)
BM6438<-ifelse(BM6438>=253.9 & BM6438<256.0,255,BM6438)
BM6438<-ifelse(BM6438>=256.0 & BM6438<257.9,257,BM6438)
BM6438<-ifelse(BM6438>=257.9 & BM6438<260.0,259,BM6438)
BM6438<-ifelse(BM6438>=260.0 & BM6438<262.0,261,BM6438)
BM6438<-ifelse(BM6438>=262.0 & BM6438<263.9,263,BM6438)
BM6438<-ifelse(BM6438>=263.9 & BM6438<265.8,265,BM6438)
BM6438<-ifelse(BM6438>=265.8 & BM6438<267.4,267,BM6438)
BM6438<-ifelse(BM6438>=267.4 & BM6438<269.7,269,BM6438)
BM6438<-ifelse(BM6438>=269.7 & BM6438<272.0,271,BM6438)
BM6438<-ifelse(BM6438>=272.0 & BM6438<273.2,272,BM6438)
BM6438<-ifelse(BM6438>=273.3 & BM6438<275.2,274,BM6438)
BM6438<-ifelse(BM6438>=275.2 & BM6438<277.0,276,BM6438)
BM6438<-ifelse(BM6438>=277.0 & BM6438<279.0,278,BM6438)
BM6438<-ifelse(BM6438>=279.0 & BM6438<281.1,280,BM6438)
####

##INRA011 Bin##
par(mfrow=c(1,2))
hist(as.numeric(INRA011[,4:5]),main="INRA011",xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data,"INRA011")

INRA011.1<-INRA011
INRA011<-INRA011[,4:5]
INRA011<-ifelse(INRA011>=189.0 & INRA011<190.0,189,INRA011)
INRA011<-ifelse(INRA011>=190.0 & INRA011<192.0,191,INRA011)
INRA011<-ifelse(INRA011>=192.0 & INRA011<194.0,193,INRA011)
INRA011<-ifelse(INRA011>=194.0 & INRA011<196.0,195,INRA011)
INRA011<-ifelse(INRA011>=196.0 & INRA011<198.0,197,INRA011)
INRA011<-ifelse(INRA011>=198.0 & INRA011<200.0,199,INRA011)
INRA011<-ifelse(INRA011>=200.0 & INRA011<202.0,201,INRA011)
INRA011<-ifelse(INRA011>=202.0 & INRA011<204.0,203,INRA011)
INRA011<-ifelse(INRA011>=204.0 & INRA011<206.0,205,INRA011)
INRA011<-ifelse(INRA011>=206.0 & INRA011<208.0,207,INRA011)
####



##BL42 Bin##
par(mfrow=c(1,1))
hist(as.numeric(BL42[,4:5]), main="BL42", xlab="Raw Allele",breaks=10000)
AlleleCum(gen.data, "BL42")


BL42.1<-BL42
BL42<-BL42[,4:5]
BL42<-ifelse(BL42>=234.5 & BL42<236.0,235,BL42)
BL42<-ifelse(BL42>=240.0 & BL42<241.5,241,BL42)
BL42<-ifelse(BL42>=242.5 & BL42<243.5,243,BL42)
BL42<-ifelse(BL42>=244.4 & BL42<245.7,245,BL42)
BL42<-ifelse(BL42>=246.5 & BL42<247.7,247,BL42)
BL42<-ifelse(BL42>=248.0 & BL42<249.7,249,BL42)
BL42<-ifelse(BL42>=249.8 & BL42<251.5,251,BL42)
BL42<-ifelse(BL42>=251.8 & BL42<253.8,253,BL42)
BL42<-ifelse(BL42>=253.8 & BL42<255.5,255,BL42)
BL42<-ifelse(BL42>=255.8 & BL42<257.0,256,BL42)
BL42<-ifelse(BL42>=257.5 & BL42<258.8,258,BL42)
BL42<-ifelse(BL42>=259.0 & BL42<260.5,260,BL42)
BL42<-ifelse(BL42>=261.5 & BL42<262.7,262,BL42)
BL42<-ifelse(BL42>=263.4 & BL42<264.5,264,BL42)
BL42<-ifelse(BL42>=265.0 & BL42<266.5,266,BL42)
######


###RT5
par(mfrow=c(1,2))
hist(as.numeric(RT5[,4:5]), main="RT5", xlab="Raw Allele", breaks=10000)
AlleleCum(gen.data, "RT5")


RT5.1<-RT5
RT5<-RT5[,4:5]
RT5<-ifelse(RT5>=97.5 & RT5<99.0,98,RT5)
RT5<-ifelse(RT5>=99.0 & RT5<101.0,100,RT5)
RT5<-ifelse(RT5>=101.0 & RT5<103.0,102,RT5)
RT5<-ifelse(RT5>=103.0 & RT5<105.0,104,RT5)
RT5<-ifelse(RT5>=106.0 & RT5<107.1,106,RT5)
RT5<-ifelse(RT5>=108.0 & RT5<109.5,109,RT5)
RT5<-ifelse(RT5>=110.0 & RT5<111.7,111,RT5)
RT5<-ifelse(RT5>=111.9 & RT5<113.5,113,RT5)
RT5<-ifelse(RT5>=114.0 & RT5<115.5,115,RT5)
RT5<-ifelse(RT5>=116.0 & RT5<117.5,117,RT5)
RT5<-ifelse(RT5>=118.0 & RT5<120.0,119,RT5)
RT5<-ifelse(RT5>=120.5 & RT5<121.5,121,RT5)
RT5<-ifelse(RT5>=122.0 & RT5<123.5,123,RT5)
RT5<-ifelse(RT5>=124.5 & RT5<125.5,125,RT5)


#####OarFCB193
par(mfrow=c(1,2))
hist(as.numeric(OarFCB193[,4:5]), main="OarFCB193", xlab="Raw Allele", breaks=10000)
AlleleCum(gen.data, "OarFCB193")

OarFCB193.1<-OarFCB193
OarFCB193<-OarFCB193[,4:5]
OarFCB193<-ifelse(OarFCB193>=89.0 & OarFCB193<90.0,90,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=91.5 & OarFCB193<92.5,92,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=93.0 & OarFCB193<94.0,94,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=95.0 & OarFCB193<96.0,96,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=97.0 & OarFCB193<98.1,98,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=99.5 & OarFCB193<101.0,100,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=101.5 & OarFCB193<102.9,102,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=103.0 & OarFCB193<105.0,104,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=105.0 & OarFCB193<107.0,106,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=107.0 & OarFCB193<109.0,108,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=109.1 & OarFCB193<111.0,110,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=111.5 & OarFCB193<113.5,113,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=114.0 & OarFCB193<115.9,115,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=116.0 & OarFCB193<118.0,117,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=118.0 & OarFCB193<119.9,119,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=120.0 & OarFCB193<121.9,121,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=122.0 & OarFCB193<123.9,123,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=124.0 & OarFCB193<125.9,125,OarFCB193)
OarFCB193<-ifelse(OarFCB193>=126.0 & OarFCB193<128.0,127,OarFCB193)

##################


#########Write Formatted Allele Calls##########
Data<-cbind(row,RT9, BM4107, P,  Cervid1,  Q, RT7, BM6438, INRA011, RT5, OarFCB193, BL42)
colnames(Data)<-c("CODE","RT9a1", "RT9a2", "BM4107a1", "BM4107a2", "Pa1", "Pa2", "Cervid1a1", "Cervid1a2", "Qa1", "Qa2", "RT7a1", "RT7a2", "BM6438a1", "BM6438a2", "INRA011a1", "INRA011a2","RT5a1", "RT5a2", "OarFCB193a1", "OarFCB193a2", "BL42a1", "BL42a2")
write.csv(Data,"example_formatted.csv")
############




