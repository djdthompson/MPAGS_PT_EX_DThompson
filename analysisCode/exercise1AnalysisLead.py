import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import math
import itertools
import sys
from scipy.stats import norm, moyal

##################
#Choose mode here: Circle or Square Mag field
#mode="Circle"
mode="Square"
mode="Squ"
##################

#choose intial momentum here
initMom=100

#file="B5%s_%sGeV.root"%(mode,initMom)
file="mu_%s_%sGeV_RandomFalse_LeadBefore.root"%(mode,initMom)
mode="Square"
# f=up.open("%s:B5"%file)
# print(f.show())
# f.close()

with up.open("%s:B5"%file) as f:
    print(f.show())
    DC1=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc1HitsVector_x","y":"Dc1HitsVector_y","z":"Dc1HitsVector_z"})
    DC2=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc2HitsVector_x","y":"Dc2HitsVector_y","z":"Dc2HitsVector_z"})


#print(sorted([len(DC1["x"][i]) for i in range(len(DC1["x"]))],reverse=True))
#sys.exit()


#z chambers are 0.5m apart
#x,y in mm x

zPosDC1=np.array([-6.25,-5.75,-5.25,-4.75,-4.25])
zPosDC2=np.array([2.25,2.75,3.25,3.75,4.25])
#first will convert cm to m to make fitting and everything easier!
sigmaX=0.0001 #in m
sigmaY=0.01 #in m

#print(DC1["x"][555])

#convert cm to m
DC1["x"]=DC1["x"]/1000
DC1["y"]=DC1["y"]/1000
DC2["x"]=DC2["x"]/1000
DC2["y"]=DC2["y"]/1000

#print(DC1["x"][555])
rB=1
B=0.5
#initMom=200#GeV
#magfield defined as circle (1m)**2=x^2+z^2
#will be staiight either side,could i reflect model curve and just get r?
#need magfield bounds

rBSquare=1 
#square/cube has side length 2m centred on 0,0,0, so "radius" is 1
#it is initiated as Dx,Dy,Dz 1mx1mx1m but this is the half length of the cube
#can see in "GetCubicVolume" in geant4, a G4Box class would return 8*Dx*Dy*Dz
#ie: 4 equations: x=-1,x=+1,z=-1,z=+1, need to check whether it is within limits of other
#x=+-1 for -1<z<+1
#z=+-1 for -1<x<+1
#determine which the same way as with circle!



def plotFits(n,momentum,p=None):

    plt.figure()
    #n=np.random.randint(1000)
    plt.plot(DC1["x"][n]*1000,zPosDC1[DC1["z"][n].astype(int)],'rx',label="DC1")
    plt.plot(DC2["x"][n]*1000,zPosDC2[DC2["z"][n].astype(int)],'bx',label="DC2")
    if mode=="Circle":
        plt.plot(np.arange(-100,101,1)/100*1000,np.sqrt(1-(np.arange(-100,101,1)/100)**2),'g',label="B-Field")
        plt.plot(np.arange(-100,101,1)/100*1000,-1*np.sqrt(1-(np.arange(-100,101,1)/100)**2),'g')
    elif mode=="Square":
        plt.gca().add_patch(pat.Rectangle((-1000,-1),2000,2,ec='g',label="B-Field"))

    xZLin1=np.linspace(np.amin(DC1["x"][n]),np.amax(DC1["x"][n]),50)
    xZLin2=np.linspace(np.amin(DC2["x"][n]),np.amax(DC2["x"][n]),50)
    plt.plot(xZLin1*1000,np.polyval(p[0],xZLin1),'k',label="DC1 Hits Fit")
    plt.plot(xZLin2*1000,np.polyval(p[1],xZLin2),'k',label="DC2 Hits Fit")
    plt.xlabel("x [mm]")
    plt.ylabel("z [m]")
    plt.legend()
    plt.title("Track in x-z plane for event %i \n with recon. momentum = %.1f."%(n,momentum))

    plt.figure()

    plt.plot(zPosDC1[DC1["z"][n].astype(int)],DC1["y"][n]*1000,'rx',label="DC1")

    plt.plot(zPosDC2[DC2["z"][n].astype(int)],DC2["y"][n]*1000,'bx',label="DC2")
    #print(np.zeros(len(DC1["z"][n])+len(DC2["z"][n])))
    plt.errorbar(np.concatenate((zPosDC1[DC1["z"][n].astype(int)],zPosDC2[DC2["z"][n].astype(int)])),np.zeros(len(DC1["z"][n])+len(DC2["z"][n])),c='g',capsize=4,yerr=10,label="Expected path")
    zYLin=np.linspace(zPosDC1[0],zPosDC2[-1],50)
    plt.plot(zYLin,np.polyval(p[2],zYLin)*1000,'k',label="Z-Y Hits Fit")
    plt.xlabel("z [m]")
    plt.ylabel("y [mm]")
    plt.legend()
    plt.title("Track in z-y plane for event %i \n with recon. momentum = %.1f."%(n,momentum))
    plt.show()





#whichever has the larger gradient magnitude SHOULD be the initial beam
#but in the circle case we know 1 is the initial
#1 should get the negative z result, 2 must be positive z to reach DC2
#so initally solve for z, gives obvious result, choose correct z point... run into line equ. to get correct x
#equ: 0=z^2(1+m^2) -2cz+c^2+r^2m^2
def getInterceptsCircle(m1,c1,m2,c2):
    #inital beam:
    if (2*c1)**2-(4*(1+m1**2)*(-1*(rB**2)*(m1**2)+c1**2))<0:
        print("why")
        return 0,0,0,0
    z1Ans=(np.array([2*c1,2*c1])+np.sqrt((2*c1)**2-(4*(1+m1**2)*(-1*rB**2*m1**2+c1**2)))*np.array([1,-1]))/(2*(1+m1**2))
    #add uncertainties with undertainty package if i want
    if len(z1Ans[z1Ans<0])==0:
        #catches an obviously failed fit, doesn't catch no interception...
        return 0,0,0,0
    z1=z1Ans[z1Ans<0][0] #zero to get it as value
    x1=(z1-c1)/m1

    if ((2*c2)**2-(4*(1+m2**2)*(-1*(rB**2)*(m2**2)+c2**2)))<0:
        return 0,0,0,0
    z2Ans=(np.array([2*c2,2*c2])+np.sqrt((2*c2)**2-(4*(1+m2**2)*(-1*rB**2*m2**2+c2**2)))*np.array([1,-1]))/(2*(1+m2**2))
    #add uncertainties with undertainty package if i want
    z2=z2Ans[z2Ans>0][0] #index to get it as value not array
    x2=(z2-c2)/m2

    return x1,z1,x2,z2
#returns of zero check for failed events/failed fits,these can be left to lead to momentum of 0 or can be manually skipped


def getInterceptsSquare(m1,c1,m2,c2):
    #inital beam should go through bottom of square always, if not fail!
    #so at z=-1...
    z1=-rBSquare
    x1=(z1-c1)/m1
    if x1>rBSquare or x1<-rBSquare:
        #beam misses mag field...fail!
        return 0,0,0,0

    #for exit beam, could intersect with z=+1, x=-1 or x=+1
    #pair x=+1 with positive gradient, -1 with negative, hmmm maybe not
    #first look at z=+1
    z2=+rBSquare
    x2=(z2-c2)/m2
    if x2<=rBSquare and x2>=-rBSquare:
        pass
    elif x2<-rBSquare:
        x2=-rBSquare
        z2=m2*x2+c2
        if z2<-rBSquare or z2>+rBSquare:
            return 0,0,0,0
    else:
        x2=+rBSquare
        z2=m2*x2+c2
        if z2<-rBSquare or z2>+rBSquare:
            return 0,0,0,0

    return x1,z1,x2,z2
#returns of zero check for failed events/failed fits,these can be left to lead to momentum of 0 or can be manually skipped



def fitLine(x,z,chamberNum):

    #only need to categorise 


    if len(x)>0:
        if chamberNum==1:
            referenceVal=[np.arange(len(z),dtype=int)[z==zPosDC1[j]] for j in range(5)] #a python list lists corresponding to the index in the full coords that correspond to hits from each chamber
        if chamberNum==2:
            referenceVal=[np.arange(len(z),dtype=int)[z==zPosDC2[j]] for j in range(5)]
        #tempZ=np.array([x[z==zPosDC1[j]] for j in range(5)])
        chiSqMin=float("inf")
        #this should mean the fit should only be attempted between 5 hits in different chambers
        if len(referenceVal[0])==0 or len(referenceVal[1])==0 or len(referenceVal[2])==0 or len(referenceVal[3])==0 or len(referenceVal[4])==0:
            #if there are no hits in any of chambers, reject event for difficulty in fitting
            return None,None,None,None
        for chamber1 in referenceVal[0]:
            for chamber2 in referenceVal[1]:
                for chamber3 in referenceVal[2]:
                    for chamber4 in referenceVal[3]:
                        for chamber5 in referenceVal[4]:
        
                            elements=np.array([chamber1,chamber2,chamber3,chamber4,chamber5])
                            #print(elements)
                            p=np.polyfit(z[elements],x[elements],1,cov=False,full=False,w=1/(np.ones(5)*sigmaX)**2)
                            chiSqTest=np.sum(((x[elements]-np.polyval(p,z[elements]))/sigmaX)**2)
                            if chiSqTest<chiSqMin:
                                chiSqMin=chiSqTest
                                bestP=p
                                bestElements=elements
    else:
        return None,None,None,None

    # if chamberNum==1:
    #     chiSqAxis=np.sum(((x[bestElements])/sigmaX)**2)/len(bestElements)
    #     if chiSqAxis>1.1 or chiSqAxis<0.01:
    #         return None,None,None,None

    #if want uncertainties add here!
    m=1/bestP[0]
    c=-bestP[1]/bestP[0]
    #remember to invert final result!
    return m,c,chiSqMin,bestElements


#first fit to lines
discardedEventsPreFit=[]

discardedEventsFromXFit=[]
discardedEventsFromYFit=[]
discardedEventsFromIntercept=[]

momentumValues=[]

for i in range(len(DC1["x"])):
    if i%(len(DC1["x"])/10)==0:
        print("%i %% complete"%int(100*i/len(DC1["x"])))
    #first check for events with far too many hits to fit!
    eventX1,eventY1,eventZ1,eventX2,eventY2,eventZ2=DC1["x"][i],DC1["y"][i],zPosDC1[DC1["z"][i].astype(int)],DC2["x"][i],DC2["y"][i],zPosDC2[DC2["z"][i].astype(int)]
    if len(eventX1)>30 or len(eventX2)>30:
        #discard these events as way too many points to iterate through to find best 5! These are interesting multiple scattering events so store!
        #15 hits is about 3000 different combinations to fit,possible
        discardedEventsPreFit.append(i)
        continue

    #instead i should re-categorise hits into sections,will do it in fitLine function


    m1,c1,chiSqMin1,bestElements1=fitLine(eventX1,eventZ1,1)
    m2,c2,chiSqMin2,bestElements2=fitLine(eventX2,eventZ2,2)
    if m1==None or m2==None:
        discardedEventsFromXFit.append(i)
        continue
    #print(chiSqMin1,chiSqMin2)

    #the chisq are so small, this doesn't do much but will limit the chisq to be less than
    if chiSqMin2/(len(bestElements2)-2)>1.5 or chiSqMin1/(len(bestElements1)-2)>1.5:
        #continue without adding this event (discard event for bad fit on straight lines)
        #potentially due to missed hits or multiple scattering?
        discardedEventsFromXFit.append(i)
        continue

    #now check y fit, ideally want a reasonable straight line: errors may be too large, so only restrict reduced chisq to be <1.5
    #print(eventZ2[bestElements2])
    joinZ=np.concatenate((eventZ1[bestElements1],eventZ2[bestElements2]))
    joinY=np.concatenate((eventY1[bestElements1],eventY2[bestElements2]))
    pY,covY=np.polyfit(joinZ,joinY,1,cov=True,full=False,w=1/(np.ones(len(bestElements1)+len(bestElements2))*sigmaY)**2)
    #redChiSqY=np.sum(((joinY-np.polyval(pY,joinZ))/sigmaY)**2)/(len(bestElements1)+len(bestElements2)-2)
    redChiSqY=np.sum(((joinY-0)/sigmaY)**2)/(len(bestElements1)+len(bestElements2))

    #print(redChiSqY)
    if redChiSqY>1:
        #fit on Y is bad, potentially a back scatter that has some how passed in x-z plane, will add to separate list
        discardedEventsFromYFit.append(i)
        continue

    #now we know we want to keep these lines.. continue with calculation

    #get intercept positions
    x1,z1,x2,z2=eval("getIntercepts%s"%mode)(m1,c1,m2,c2)
    if x1==0 and z1==0 and x2==0 and z2==0:
        #no intercept with one of the lines, discard event
        discardedEventsFromIntercept.append(i)
        continue


    #will do the full calculation of r (and hence B) without approximations for 2rB=L, sintheta=theta etc
    r=0.5*np.sqrt(2*((x1-x2)**2+(z1-z2)**2)/(1-(1+((m2-m1)/(1+m2*m1))**2)**-0.5))
    #so p=0.3Ber# in units of GeV
    p=r*B*3e-1

    #if want to look for and plot certain points, look for them here
    # if p>180:
    #     plotFits(i,p,p=[[m1,c1],[m2,c2],[pY[0],pY[1]]])
    momentumValues.append(p)



#print events that have been discarded
print("Events discarded from having too many points to reasonably fit: ", discardedEventsPreFit)
print("Events discarded in z-x fit of two lines: ", discardedEventsFromXFit)
print("Events discarded in z-y fit of two lines: ", discardedEventsFromYFit)
print("Events discarded from having no intercept with B-Field: ", discardedEventsFromIntercept)


plt.figure()
h=plt.hist(momentumValues,bins="auto",fill=False,ec="b",label="Entries: %i"%len(momentumValues))

momentumValues=np.array(momentumValues)
xmin, xmax = initMom-initMom/100*35,initMom+initMom/100*70
#mu,std=norm.fit(momentumValues[(momentumValues>xmin)&(momentumValues<xmax)])
mu,std=moyal.fit(momentumValues[(momentumValues>xmin)&(momentumValues<xmax)])
xVals = np.linspace(xmin, xmax, 100)
# fitNumbers = norm.pdf(xVals, mu, std)
fitNumbers = moyal.pdf(xVals, mu, std)
boolChoose=(h[1]>xmin)&(h[1]<xmax)#bool to choose bins within fitting area
boolChooseBins=(boolChoose[:-1].astype(int)+boolChoose[1:].astype(int))==2
#boolChooseBins=[]
print(boolChooseBins)
normaliser=len(momentumValues[(momentumValues>xmin)&(momentumValues<xmax)])*np.abs(h[1][1]-h[1][0])
#chiSqGauss=np.sum((h[0][boolChooseBins]-normaliser*norm.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2/(normaliser*norm.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2)/(len(boolChooseBins)-2)

chiSqGauss=np.sum((h[0][boolChooseBins]-normaliser*moyal.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2/(normaliser*moyal.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2)/(len(boolChooseBins)-2)
#plt.plot(xVals, fitNumbers*normaliser, 'k', linewidth=2,label="Norm. Dist. Fit\n"+r'$\mu=$'+"%.1f , "%mu+r'$\sigma=$'+"%.1f\n"%std+r"$\chi^2=$"+"%.1f"%chiSqGauss)
plt.plot(xVals, fitNumbers*normaliser, 'k', linewidth=2,label="Landau Dist. Fit\n"+r'$\mu=$'+"%.1f , "%mu+r'$\sigma=$'+"%.1f\n"%std+r"$\chi^2=$"+"%.1f"%chiSqGauss)

plt.legend()
print((h[1][-1]-h[1][0])/len(h[0]))
plt.ylabel("Number of Events / (%.1f GeV)"%((h[1][-1]-h[1][0])/len(h[0])))
plt.xlabel("Reconstructed momentum [GeV]")
plt.title('Reconstructed momentum from fit of deviation of an antimuon beam\nat %i GeV in a %s B-field of B= %.2f \nwith 5cm of lead placed before mag field'%(initMom,mode,B))
plt.show()


expecAngle=0.3*0.5*2/100
plt.figure()
#vals=np.random.uniform(0,2*expecAngle,10000)
vals=np.random.uniform(expecAngle-2*4.4e-4,expecAngle+2*4.4e-4,10000)

pVals=0.3/vals
hBG=plt.hist(pVals,bins=h[1])
plt.ylabel("Number of Events / (%.1f GeV)"%((h[1][-1]-h[1][0])/len(h[0])))
plt.xlabel("momentum [GeV]")
plt.title("Reconstructed momentum for random sample of angles")
plt.show()

cast=(hBG[1][:-1]>125)&(hBG[1][1:]<150)
normBG=np.sum(h[0][cast])/np.sum(hBG[0][cast])

newHist=h[0]-hBG[0]*normBG

centre=(hBG[1][1:]+hBG[1][:-1])/2
width = hBG[1][1:] - hBG[1][:-1]
newHist[newHist<0]=np.zeros(len(newHist[newHist<0])) 
plt.figure()
plt.bar(centre,newHist, align='center', width=width*0.9,label="1/theta BG subtracted plot")
plt.ylabel("Number of Events-BG (RemoveNegatives) / (%.1f GeV)"%((h[1][-1]-h[1][0])/len(h[0])))
plt.xlabel("momentum [GeV]")
plt.title('Reconstructed momentum from fit of deviation of an antimuon beam\nat %i GeV in a %s B-field of B= %.2f \nwith 5cm of lead placed before mag field'%(initMom,mode,B))
# xmin, xmax = mu-50,mu+50

# fitVals=[]
# i=0
# for a in newHist[(newHist>xmin)&(newHist<xmax)]:
#     for b in range(np.floor(a).astype(int)):
#         fitVals.append(centre[(newHist>xmin)&(newHist<xmax)][i])
#     i+=1
# mu,std=norm.fit(fitVals)
# xVals = np.linspace(xmin, xmax, 100)
# fitNumbers = norm.pdf(xVals, mu, std)
# boolChoose=(h[1]>xmin)&(h[1]<xmax)#bool to choose bins within fitting area
# boolChooseBins=(boolChoose[:-1].astype(int)+boolChoose[1:].astype(int))==2
# normaliser=np.sum(newHist[(newHist>xmin)&(newHist<xmax)])*np.abs(h[1][1]-h[1][0])
# #chiSqShift=np.sum((h[0][boolChooseBins]-normaliser*moyal.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2/(normaliser*moyal.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2)/(len(boolChooseBins)-2)
# plt.plot(xVals, fitNumbers*normaliser, 'k', linewidth=2,label="Gauss. Fit\n"+r'$\mu=$'+"%.1f , "%mu+r'$\sigma=$'+"%.1f\n"%std)#+r"$\chi^2=$"+"%.1f"%chiSqShift)

plt.legend()

plt.show()




