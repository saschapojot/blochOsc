from pbcFuncs import *

import inspect


tStart=datetime.now()
for q in range(0,Q):
    Aq=AAll[q]
    Bq=BAll[q]
    Aq=reNormalization(Aq)
    Bq=reNormalization(Bq)
    AqNext,BqNext=RK4(q,Aq,Bq)
    AAll.append(AqNext)
    BAll.append(BqNext)

tEnd=datetime.now()
print("computation time: ", tEnd - tStart)
outDir = "./fig2/"
dataAll=[]
for q in range(0,Q+1):
    psiq=[]
    Aq=AAll[q]
    Bq=BAll[q]
    for j in range(0,len(Aq)):
        psiq.append(Aq[j])
        psiq.append(Bq[j])
    dataAll.append(psiq)
outCsvName=outDir+"wv.csv"
dataAll=np.asarray(dataAll)
# dtFrame=pd.DataFrame(data=dataAll)
# dtFrame.to_csv(outCsvName,index=False)
np.savetxt(outCsvName,dataAll,delimiter=",")
# # plotting wavefunctions
# plotStart=datetime.now()
# for q in range(0,Q):
#     plt.figure()
#     nmTmp=[np.abs(elem)**2 for elem in dataAll[q]]
#     plt.plot(range(0,L),nmTmp,color="black")
#     plt.xlabel("position"
#                )
#     plt.title("time = "+str(q*dt)+", g = "+str(g))
#     plt.ylabel("magnitude")
#     # outFile=outDir+"omegaF"+str(omegaF)+"omega"+str(omega)+"g"+str(g)+"q"+str(q)+".png"
#     outFile=outDir+"q"+str(q)+".png"
#     plt.savefig(outFile)
#     plt.close()
# plotEnd=datetime.now()
# print("plotting time: ",plotEnd-plotStart)
###### end plotting wave functions
xPos = []
# wdAll=[]
for q in range(0, Q+1):
    vecTmp = dataAll[q]
    # vecTmp=reNormalization(vecTmp)
    xTmp = meanXAndXWd(vecTmp)
    xPos.append(xTmp)
    # wdAll.append(wdTmp)
drift=[elem-xc for elem in xPos]
posMax = np.max(drift)
posMin = np.min(drift)

posDiff = 0.1
tickNum = int((posMax - posMin) / posDiff)
yTicks = [posMin + j * posDiff for j in range(0, tickNum + 2)]
tAll = [dt * q for q in range(0, Q+1)]
plt.figure(figsize=(20,20))
# plt.yticks(yTicks)
plt.plot(tAll, drift, color="black")
plt.xlabel("time")
plt.ylabel("avg position")
plt.title("g = " + str(g))
plt.savefig(outDir + "g" + str(g) + "position.png")
plt.close()


outTxt = outDir + "info.txt"

fptr = open(outTxt, "w+")
fptr.write("Total drift = "+str(drift[-1]-drift[0])+"\n")
fptr.write("g=" + str(g) + "\n")
fptr.write("omega=" + str(omega) + "\n")
fptr.write("omegaF=" + str(omegaF) + "\n")
fptr.write(inspect.getsource(x))
fptr.write(inspect.getsource(y))
fptr.write(inspect.getsource(u))
fptr.write(inspect.getsource(v))
fptr.write(inspect.getsource(w))
fptr.close()