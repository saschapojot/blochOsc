from funcs import *
import matplotlib.pyplot as plt
import inspect
tStart = datetime.now()
dataAll=[]
dataAll.append(psi0)

for q in range(0,Q):
    psiCurr=dataAll[q]
    psiNext=S2(psiCurr,q)
    dataAll.append(psiNext)


tEnd = datetime.now()
print("computation time: ", tEnd - tStart)
outDir = "./pump0/"

plotStart=datetime.now()
for q in range(0,Q):
    plt.figure()
    nmTmp=[np.abs(elem)**2 for elem in dataAll[q]]
    plt.plot(range(0,L),nmTmp,color="black")
    plt.xlabel("position"
               )
    plt.title("time = "+str(q*dt)+", g = "+str(g))
    plt.ylabel("magnitude")
    # outFile=outDir+"omegaF"+str(omegaF)+"omega"+str(omega)+"g"+str(g)+"q"+str(q)+".png"
    outFile=outDir+"q"+str(q)+".png"
    plt.savefig(outFile)
    plt.close()
plotEnd=datetime.now()
print("plotting time: ",plotEnd-plotStart)
xPos = []
for q in range(0, Q):
    vecTmp = dataAll[q]
    xTmp = meanX(vecTmp)
    xPos.append(xTmp)

posMax = np.max(xPos)
posMin = np.min(xPos)

posDiff = 0.1
tickNum = int((posMax - posMin) / posDiff)
yTicks = [posMin + j * posDiff for j in range(0, tickNum + 2)]
tAll = [dt * q for q in range(0, Q)]
plt.figure(figsize=(20,20))
plt.yticks(yTicks)
plt.plot(tAll, xPos, color="black")
plt.xlabel("time")
plt.ylabel("ave position")
plt.title("g = " + str(g))
plt.savefig(outDir + "g" + str(g) + "position.png")
plt.close()

# write params info
outTxt = outDir + "info.txt"

fptr = open(outTxt, "w+")
fptr.write("g=" + str(g) + "\n")
fptr.write("omega=" + str(omega) + "\n")
fptr.write("omegaF=" + str(omegaF) + "\n")
fptr.write(inspect.getsource(x))
fptr.write(inspect.getsource(y))
fptr.write(inspect.getsource(u))
fptr.write(inspect.getsource(v))
fptr.write(inspect.getsource(w))
fptr.close()
