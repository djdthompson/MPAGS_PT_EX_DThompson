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
mode="Squ"
##################

#choose intial momentum here
initMom=100

#choose particle: prot,mu,elec
particle="elec"
#mu_Squ_100GeV_RandomFalse
file="%s_%s_%sGeV_RandomFalse.root"%(particle,mode,initMom)

names={"prot":"protons","elec":"positrons","mu":"antimuons"}

# f=up.open("%s:B5"%file)
# print(f.show())
# f.close()

with up.open("%s:B5"%file) as f:
    print(f.show())
    DC1=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc1HitsVector_x","y":"Dc1HitsVector_y","z":"Dc1HitsVector_z"})
    DC2=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc2HitsVector_x","y":"Dc2HitsVector_y","z":"Dc2HitsVector_z"})
    EC=f.arrays(["E","V"],library="np",aliases={"E":"ECEnergy","V":"ECEnergyVector"})
    HC=f.arrays(["E","V"],library="np",aliases={"E":"HCEnergy","V":"HCEnergyVector"})    

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

#measured in MeV
#V is vector corresponding to cells in calorimeter, 80 for ECAL, 20 for HCAL
# for i in range(0,1000):
#     #a=np.random.randint(1000)
#     #print(EC["E"][a]/np.sum(EC["V"][a]))
#     if(np.abs(HC["E"][i]-np.sum(HC["V"][i]))<1):
#         print(i)
#         #print
        
#     # print(i)
#     # print(HC["E"][i])
#     # print(HC["V"][i])
#     # print(np.sum(HC["V"][i]))

# for i in range(0,1000):
#     #a=np.random.randint(1000)
#     #print(EC["E"][a]/np.sum(EC["V"][a]))
#     # if(np.abs(EC["E"][i]-np.sum(EC["V"][i]))<1):
#     #     print(i)
#     #     #print
#     #     break
#     # print(i)
#     # print(HC["E"][i])
#     # print(HC["V"][i])
#     # print(np.sum(EC["V"][i]))

#     for j in range(20):

#         if EC["V"][i][j]<0:
#             print(i)
#             print(EC["V"][i])

#going to just go ahead with using the vector of hits

#plt.figure()
#hitsECal=EC["V"][0]
#split into 20x4 to rep x,y
EC["V"]=np.array(EC["V"])/1000
HC["V"]=np.array(HC["V"])/1000
#print(EC["V"][0])
hitsSplitECal=np.zeros((1000,4,20))
maxIndicesECal=np.zeros((1000,2),dtype=int)
hitsSplitHCal=np.zeros((1000,2,10))
maxIndicesHCal=np.zeros((1000,2),dtype=int)
for i in range(len(EC["V"])):
    hitsSplitECal[i]=np.split(EC["V"][i],4)
    maxIndicesECal[i]=np.unravel_index(np.argmax(hitsSplitECal[i]),(4,20))
    hitsSplitHCal[i]=np.split(HC["V"][i],2)
    maxIndicesHCal[i]=np.unravel_index(np.argmax(hitsSplitHCal[i]),(2,10))
#maxIndices=np.unravel_index(np.argmax(hitsSplitECal,axis=1),hitsSplitECal[0].shape)
# print(hitsSplitECal[0])
# print(maxIndices[0,0])
# print(hitsSplitECal[0,maxIndices[0,0],maxIndices[0,1]])
# print(maxIndices[:,1])
if 0:
    plt.figure()
    hitsHistECAL=plt.hist2d(maxIndicesECal[:,1],maxIndicesECal[:,0],bins=[np.arange(21),np.arange(5)],cmap="viridis")
    plt.xticks(np.arange(21).astype(int))
    plt.yticks(np.arange(5).astype(int))
    plt.xlabel("Cells along x")
    plt.ylabel("Cells along y")
    cB=plt.colorbar()
    cB.set_label("Number of Events")
    plt.title("Plot showing distribution of Maximum energy deposition in ECAL for 1000 %s\nAverage maximal energy deposition: %.1f GeV"%(names[particle],np.mean(hitsSplitECal[:,maxIndicesECal[0,0],maxIndicesECal[0,1]])))
    plt.show()

    plt.figure()
    hitsHistHCAL=plt.hist2d(maxIndicesHCal[:,1],maxIndicesHCal[:,0],bins=[np.arange(11),np.arange(3)],cmap="viridis")
    plt.xticks(np.arange(11).astype(int))
    plt.yticks(np.arange(3).astype(int))
    plt.xlabel("Cells along x")
    plt.ylabel("Cells along y")
    cB=plt.colorbar()
    cB.set_label("Number of Events")
    plt.title("Plot showing distribution of Maximum energy deposition in HCAL for 1000 %s\nAverage maximal energy deposition: %.1f GeV"%(names[particle],np.mean(hitsSplitHCal[:,maxIndicesHCal[0,0],maxIndicesHCal[0,1]])))
    plt.show()



    #hitsSplitECal=np.array(hitsSplitECal)
    plt.figure()
    a=np.random.randint(1000)
    plt.imshow(hitsSplitECal[a])
    plt.xlabel("Cells along x")
    plt.ylabel("Cells along y")
    plt.title("Plot showing distribution of hits in ECAL for a proton event %i"%a)
    cB=plt.colorbar()
    cB.set_label("Energy Deposition GeV")
    plt.show()
    #hitsSplitECal=np.array(hitsSplitECal)
    plt.figure()
    #a=np.random.randint(1000)
    plt.imshow(hitsSplitHCal[a])
    plt.xlabel("Cells along x")
    plt.ylabel("Cells along y")
    plt.title("Plot showing distribution of hits in HCAL for a proton event %i"%a)
    cB=plt.colorbar()
    cB.set_label("Energy Deposition GeV")
    plt.show()



totalEnergy=np.array([np.sum(EC["V"][i])+np.sum(HC["V"][i]) for i in range(len(EC["V"]))])
totalECAL=np.array([np.sum(EC["V"][i]) for i in range(len(EC["V"]))])
totalHCAL=np.array([np.sum(HC["V"][i]) for i in range(len(EC["V"]))])

#ignore any with total energy more than 1.1e5, definitely wrong!
fig=plt.figure()
ax1=fig.add_subplot((211))
ax2=fig.add_subplot((212))

hECAL=ax1.hist(totalECAL[totalECAL<1.1e5],bins="auto",fill=False,ec="b",alpha=0.5,label="ECAL\nEntries: %i"%len(totalECAL[totalECAL<1.1e5]))
hHCAL=ax2.hist(totalHCAL[totalHCAL<1.1e5],bins="auto",fill=False,ec="r",alpha=0.5,label="HCAL\nEntries: %i"%len(totalHCAL[totalHCAL<1.1e5]))
ax1.set_xlabel("Total Energy Deposited in Electromagnetic Calorimeter [GeV]")
ax1.set_ylabel("Number of Events")# / (%i MeV)"%(int(hEnergy[1][-1]-hEnergy[1][0])/len(hEnergy[0])))
ax2.set_xlabel("Total Energy Deposited in Hadronic Calorimeter [GeV]")
ax2.set_ylabel("Number of Events")
mu1,std1=moyal.fit(100-totalECAL[totalECAL<1.1e5]) #this for electron
#mu1,std1=moyal.fit(totalECAL[totalECAL<1.1e5]) #this for proton and muon
mu2,std2=moyal.fit(totalHCAL[totalHCAL<1.1e5]) #this for electron and muon
#mu2,std2=moyal.fit(100-totalHCAL[totalHCAL<1.1e5]) #this for proton
xVals1=100-((hECAL[1][1:]+hECAL[1][:-1])/2) #this for electron
#xVals1=((hECAL[1][1:]+hECAL[1][:-1])/2) #this for proton and muon
xVals2=(hHCAL[1][1:]+hHCAL[1][:-1])/2 #this for electron and muon
#xVals2=100-(hHCAL[1][1:]+hHCAL[1][:-1])/2 #this for proton

normaliser1=len(totalECAL[totalECAL<1.1e5])*np.abs(hECAL[1][1]-hECAL[1][0])
fitVals1=moyal.pdf(xVals1,mu1,std1)*normaliser1

normaliser2=len(totalHCAL[totalHCAL<1.1e5])*np.abs(hHCAL[1][1]-hHCAL[1][0])
fitVals2=moyal.pdf(xVals2,mu2,std2)*normaliser2
#chiSqGauss=np.sum((h[0][boolChooseBins]-normaliser*norm.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2/(normaliser*norm.pdf(((h[1][boolChoose])[:-1]+(h[1][boolChoose])[1:])/2, mu, std))**2)/(len(boolChooseBins)-2)

chiSq1=np.sum(((hECAL[0][hECAL[0]>0]-fitVals1[hECAL[0]>0])/np.sqrt(hECAL[0][hECAL[0]>0]))**2)/(len(hECAL[0])-2)
chiSq2=np.sum(((hHCAL[0][hHCAL[0]>0]-fitVals2[hHCAL[0]>0])/np.sqrt(hHCAL[0][hHCAL[0]>0]))**2)/(len(hHCAL[0])-2)

ax1.plot(np.linspace(hECAL[1][0],hECAL[1][-1],1000),normaliser1*moyal.pdf(100-np.linspace(hECAL[1][0],hECAL[1][-1],1000),mu1,std1),'k', linewidth=2,label="Landau Dist. Fit\n"+r'$\mu=$'+"%.1f , "%(100-mu1)+r'$\sigma=$'+"%.1f\n"%std1+r"$\chi^2=$"+"%.1f"%chiSq1)
#this for electron

#ax1.plot(np.linspace(hECAL[1][0],hECAL[1][-1],1000),normaliser1*moyal.pdf(np.linspace(hECAL[1][0],hECAL[1][-1],1000),mu1,std1),'k', linewidth=2,label="Landau Dist. Fit\n"+r'$\mu=$'+"%.1f , "%(mu1)+r'$\sigma=$'+"%.1f\n"%std1+r"$\chi^2=$"+"%.1f"%chiSq1)
#this for proton and muon

ax2.plot(np.linspace(hHCAL[1][0],hHCAL[1][-1],1000),normaliser2*moyal.pdf(np.linspace(hHCAL[1][0],hHCAL[1][-1],1000),mu2,std2),'k', linewidth=2,label="Landau Dist. Fit\n"+r'$\mu=$'+"%.2f , "%mu2+r'$\sigma=$'+"%.2f\n"%std2+r"$\chi^2=$"+"%.1f"%chiSq2+"\n$f_{samp}$ adjusted "+r"$\mu=$"+" %.2f"%(mu2*5))
#this for electron and muon

#ax2.plot(np.linspace(hHCAL[1][0],hHCAL[1][-1],1000),normaliser2*moyal.pdf(100-np.linspace(hHCAL[1][0],hHCAL[1][-1],1000),mu2,std2),'k', linewidth=2,label="Landau Dist. Fit\n"+r'$\mu=$'+"%.2f , "%(100-mu2)+r'$\sigma=$'+"%.2f\n"%std2+r"$\chi^2=$"+"%.1f"%chiSq2+"\n$f_{samp}$ adjusted "+r"$\mu=$"+" %.2f"%((100-mu2)*5))
#this for proton


ax1.legend()
ax2.legend()
fig.suptitle("Total Energy Deposited in Calorimeters for 1000, 100 GeV %s"%names[particle])
fig.show()


# plt.figure()
# hitsHCal=HC["V"][0]
# #split into 10x4 to rep x,y
# hitsSplitHCal=np.split(hitsHCal,2)
# plt.imshow(hitsSplitHCal)
# plt.show()


#plt.figure()

#plt.hist2d(np.sum(hitsSplitECal,axis=0),np.sum(hitsSplitECal,axis=1),cmap="viridis")

#print(np.unravel_index(np.argmax(hitsSplitECal[0]),hitsSplitECal[0].shape))
#print(hitsSplitECal[np.argmax(hitsSplitECal)%len(hitsSplitECal),np.argmax(hitsSplitECal)%len(hitsSplitECal)])
#plt.colormap()
#plt.show()
