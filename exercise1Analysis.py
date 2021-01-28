import uproot as up
import numpy as np
import matplotlib.pyplot as plt

file="B5.root"


# f=up.open("%s:B5"%file)
# print(f.show())
# f.close()

with up.open("%s:B5"%file) as f:
    print(f.show())
    DC1=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc1HitsVector_x","y":"Dc1HitsVector_y","z":"Dc1HitsVector_z"})
    DC2=f.arrays(["x","y","z"],library="np",aliases={"x":"Dc2HitsVector_x","y":"Dc2HitsVector_y","z":"Dc2HitsVector_z"})


print(DC1["x"][34][1])
#z chambers are 0.5m apart
#x,y in cm x
sigmaX=0.01 #in cm
sigmaY=1 #in cm
zPosDC1=np.array([-6.25,-5.75,-5.25,-4.75,-4.25])
zPosDC2=np.array([2.25,2.75,3.25,3.75,4.25])

#magfield defined as circle radius 1m center (0,0) in x-z plane

#gab between chamber centres is 10m

#will be staiight either side,could i reflect model curve and just get r?
#need magfield bounds


plt.figure()
n=np.random.randint(1000)
plt.plot(DC1["x"][n],zPosDC1[DC1["z"][n].astype(int)],'rx',label="DC1")
plt.plot(DC2["x"][n],zPosDC2[DC2["z"][n].astype(int)],'bx',label="DC2")
plt.plot(np.arange(-100,101,1),np.sqrt(100**2-np.arange(-100,101,1)**2)/100,'g',label="B-Field")
plt.plot(np.arange(-100,101,1),-1*np.sqrt(100**2-np.arange(-100,101,1)**2)/100,'g')
plt.xlabel("x [cm]")
plt.ylabel("z [m]")
plt.legend()
plt.title("Track in x-z plane for event %i."%n)

plt.figure()

plt.plot(zPosDC1[DC1["z"][n].astype(int)],DC1["y"][n],'rx',label="DC1")

plt.plot(zPosDC2[DC2["z"][n].astype(int)],DC2["y"][n],'bx',label="DC2")
print(np.zeros(len(DC1["z"][n])+len(DC2["z"][n])))
plt.errorbar(np.concatenate((zPosDC1[DC1["z"][n].astype(int)],zPosDC2[DC2["z"][n].astype(int)])),np.zeros(len(DC1["z"][n])+len(DC2["z"][n])),c='g',capsize=4,yerr=1,label="Expected path")

plt.xlabel("z [m]")
plt.ylabel("y [cm]")
plt.legend()
plt.title("Track in z-y plane for event %i."%n)
plt.show()
