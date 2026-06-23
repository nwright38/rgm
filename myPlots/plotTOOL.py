import matplotlib.pyplot as plt
import ROOT
import numpy as np
import math
from matplotlib import ticker, cm
from scipy.ndimage.filters import gaussian_filter1d


def plotTGEStep(ax,myFile,myName,mycolor):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    e=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        e.append( myGraph.GetErrorY(i))    
    ax.errorbar(x,y,e,color=mycolor,linestyle='')
    ax.step(x,y,color=mycolor,where='mid')

def plotTGEStepEx(ax,myFile,myName,mycolor,factor,shift,mymarker):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    e=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i]+shift)
        y.append(factor*yP[i])
        e.append(factor*myGraph.GetErrorY(i))    
    ax.errorbar(x,y,e,color=mycolor,linestyle='',marker=mymarker)
    
def plotTGEStepQ2(ax,myFile,myName,mycolor):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    e=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        e.append( myGraph.GetErrorY(i))    

    x.insert(0,0)
    y.insert(0,y[0])
    e.insert(0,e[0])
    x.append(10)
    y.append(y[-1])
    e.append(e[-1])
    ax.errorbar(x,y,e,color=mycolor,linestyle='')
    ax.step(x,y,color=mycolor,where='mid')

def plotTGELine(ax,myFile,myName,mycolor,factor):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i]*factor)
        ylow.append((yP[i]-myGraph.GetErrorYlow(i))*factor)
        yhigh.append((yP[i]+myGraph.GetErrorYhigh(i))*factor)
    ax.plot(x,y,color=mycolor)
    ax.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=mycolor)

def plotTGELineQ2(ax,myFile,myName,mycolor):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    x.insert(0,0)
    y.insert(0,y[0])
    ylow.insert(0,ylow[0])
    yhigh.insert(0,yhigh[0])
    x.append(10)
    y.append(y[-1])
    ylow.append(ylow[-1])
    yhigh.append(yhigh[-1])
    ax.plot(x,y,color=mycolor)
    ax.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=mycolor)
    
def plotTGraphAsymmErrors(myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)

def plotTGraphAsymmErrorsSmoothQ5(myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    
    y = gaussian_filter1d(y, sigma=1.5)
    ylow = gaussian_filter1d(ylow, sigma=1.5)
    yhigh = gaussian_filter1d(yhigh, sigma=1.5)
    x = np.append(x,5)
    y = np.append(y,y[-1])
    ylow = np.append(ylow,ylow[-1])
    yhigh = np.append(yhigh,yhigh[-1])
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)

def plotTGraphAsymmErrorsSmoothQ5sig(myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    
    #y = gaussian_filter1d(y, sigma=1.5)
    #ylow = gaussian_filter1d(ylow, sigma=1.5)
    #yhigh = gaussian_filter1d(yhigh, sigma=1.5)
    x = np.append(x,5)
    y = np.append(y,y[-1])
    ylow = np.append(ylow,ylow[-1])
    yhigh = np.append(yhigh,yhigh[-1])
    y*=1000*135/150
    ylow*=1000*135/150
    yhigh*=1000*135/150
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)

def plotTGraphAsymmErrorsSmooth(myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
        
    #y = gaussian_filter1d(y, sigma=1.5)
    #ylow = gaussian_filter1d(ylow, sigma=1.5)
    #yhigh = gaussian_filter1d(yhigh, sigma=1.5)
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)

def plotTGraphAsymmErrorssig(myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-2*myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+2*myGraph.GetErrorYhigh(i))

    x=np.array(x)
    y=np.array(y)
    ylow=np.array(ylow)
    yhigh=np.array(yhigh)
    x*=1000
    y*=7000
    ylow*=7000
    yhigh*=7000
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)
    
def plotTGraphAsymmErrorsData(myFile,myName,lc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    xe=[]
    ye=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        xe.append(myGraph.GetErrorX(i))
        ye.append(myGraph.GetErrorY(i))
    plt.errorbar(x,y,xerr=xe,yerr=ye,linestyle='',fmt=lc)

def plotTGraphAsymmErrorsAlpha(myFile,myName,lc,fc,myalpha):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    plt.plot(x,y,lc)
    plt.fill_between(x,ylow,yhigh,alpha=myalpha,facecolor=fc)

    
def FindMin(myFile,myName,delChiSq,scale):
    myGraph=myFile.Get(myName)
    x=myGraph.GetX()
    y=myGraph.GetY()
    xMin=0
    yMin=10
    for i in range(myGraph.GetN()):
        if(y[i]<yMin):
            xMin=x[i]
            yMin=y[i]
    xLow=0
    yLow=0
    dLow=100
    xHigh=0
    yHigh=0
    dHigh=100
    for j in range(myGraph.GetN()):
        diff=math.pow((yMin+delChiSq)-y[j],2)
        if( (diff<dLow) and (x[j]<xMin) ):
            xLow=x[j]
            yLow=y[j]
            dLow=diff
        if( (diff<dHigh) and (x[j]>xMin) ):
            xHigh=x[j]
            yHigh=y[j]
            dHigh=diff
    one=round(xMin*scale,3)
    two=round((xMin-xLow)*scale,3)
    three=round((xHigh-xMin)*scale,3)
    print(one,"-",two,"+",three)
        
def getXforGraph(myGraph):
    result=[]
    axis=myGraph.GetXaxis()
    for i in range(1,(axis.GetNbins()+1)):
        result.append(axis.GetBinCenter(i))
    return result

def getXforGraphSpacing(myGraph,spacing):
    result=[]
    axis=myGraph.GetXaxis()
    smalldiff=(axis.GetBinWidth(2))/(spacing)
    minR=axis.GetXmin()+smalldiff
    maxR=axis.GetXmax()-smalldiff
    stepR=(maxR-minR)/spacing
    for i in range(spacing+1):
        result.append(minR + (stepR * i))
    return result

def getYforGraph(myGraph):
    result=[]
    axis=myGraph.GetYaxis()
    for i in range(1,(axis.GetNbins()+1)):
        result.append(axis.GetBinCenter(i))
    return result

def getYforGraphSpacing(myGraph,spacing):
    result=[]
    axis=myGraph.GetYaxis()
    smalldiff=(axis.GetBinWidth(2))/(spacing)
    minR=axis.GetXmin()+smalldiff
    maxR=axis.GetXmax()-smalldiff
    stepR=(maxR-minR)/spacing
    for i in range(spacing+1):
        result.append(minR + (stepR * i))
    return result

def getZforGraphWithArrays(myGraph,xArray,yArray):
    result=[]
    for y in yArray:
        result.append([])
        for x in xArray:
            z=myGraph.Interpolate(x,y)
            result[yArray.index(y)].append(z)
    return result

def getZListforGraphWithArrays(myGraph,xArray,yArray):
    result=[]
    for y in yArray:
        for x in xArray:
            z=myGraph.Interpolate(x,y)
            result.append(z)
    return result

def getMinZ(myGraph):
    x=getXforGraph(myGraph)
    y=getYforGraph(myGraph)
    zList = getZListforGraphWithArrays(myGraph,x,y)
    return min(zList)

def getMinZSpacing(myGraph,xSpacing,ySpacing):
    x=getXforGraphSpacing(myGraph,xSpacing)
    y=getYforGraphSpacing(myGraph,ySpacing)
    zList = getZListforGraphWithArrays(myGraph,x,y)
    return min(zList)

def getZforGraph(myGraph):
    x=getXforGraph(myGraph)
    y=getYforGraph(myGraph)
    return getZforGraphWithArrays(myGraph,x,y)

def plotContourTest(fig1,ax1,myFile,name,levels,mycolors,mylinestyle,xFactor,yFactor):
    print(name)
    graph = myFile.Get(name)
    x = getXforGraph(graph)
    y = getYforGraph(graph)
    z = np.array(getZforGraph(graph))
    x = np.array(x)*xFactor
    y = np.array(y)*yFactor
    z[:] -= getMinZ(graph)
    CS = ax1.contour(x, y, z,levels,colors=mycolors,linestyles=mylinestyle)
    #newZ=fixZ(z)
    z[:]=np.exp(-z[:])
    CS = ax1.contourf(x,y,z,15,cmap=plt.cm.jet)

def plotNC(fig1,ax1,myFile,name):
    print(name)
    graph = myFile.Get(name)
    x = getXforGraph(graph)
    y = getYforGraph(graph)
    z = np.array(getZforGraph(graph))
    z[:] -= getMinZ(graph)
    z[:]=np.exp(-z[:])
    plt.imshow(z,extent=[0.0,0.06,0.04,0.17])

def plotVecPoint(myFile,name,x,xfactor,y,yfactor):
    vec = myFile.Get(name)
    plt.plot(vec[x]*xfactor,vec[y]*yfactor,"bo")
    
def fixZ(z):
    rows = z.shape[0]
    cols = z.shape[1]
    for i in range(rows):
        for j in range(cols):
            if( z[i][j] < 2.3 ):
                z[i][j] = 2.3-z[i][j]
            else:
                z[i][j]=0
    return z
    
def plotContourSpacing(fig1,ax1,myFile,name,levels,mycolors,xSpacing,ySpacing,xFactor,yFactor):
    graph = myFile.Get(name)
    x = getXforGraphSpacing(graph,xSpacing)
    y = getYforGraphSpacing(graph,ySpacing)
    z = np.array(getZforGraphWithArrays(graph,x,y))
    x=np.array(x)*xFactor
    y=np.array(y)*yFactor
    z[:] -= getMinZSpacing(graph,xSpacing,ySpacing)
    CS = ax1.contour(x, y, z,levels,colors=mycolors)
    z[:] += getMinZSpacing(graph,xSpacing,ySpacing)
    z[:]=np.exp(-z[:])
    CS = plt.contourf(x,y,z,15,cmap=plt.cm.jet)

#####################################################

def gXH(f1,name):
    histogram = f1.Get(name)    
    return getXforHist(histogram)

def gYH(f1,name):
    histogram = f1.Get(name)    
    return getYforHist(histogram,1)

def getXLeftforHist(histogram):
    result=[]
    for i in range(1,(histogram.GetNbinsX()+1)):
        result.append(histogram.GetBinLowEdge(i))
    return result

def getXforHist(histogram):
    result=[]
    for i in range(1,(histogram.GetNbinsX()+1)):
        result.append(histogram.GetBinCenter(i))
    return result

def getYforHist(histogram,factor):
    result=[]
    for i in range(1,(histogram.GetNbinsX()+1)):
        result.append(histogram.GetBinContent(i)*factor)
    return result

def getYforBin(histogram,mybin):
    return histogram.GetBinContent(mybin)

def getYErrforHist(histogram,factor):
    result=[]
    for i in range(1,(histogram.GetNbinsX()+1)):
        result.append(histogram.GetBinError(i)*factor)
    return result

def printHist(f1,name,factor):
    histogram = f1.Get(name)    
    for i in range(1,(histogram.GetNbinsX()+1)):
        x=histogram.GetBinCenter(i)
        y=histogram.GetBinContent(i)*factor
        e=histogram.GetBinError(i)*factor
        print(round(x,3),y,e)

def getFactorHist(histogram1,histogram2):
    y1=getYforHist(histogram1,1)
    I1=0
    for j in range(0, len(y1)):
        I1 = I1 + y1[j]

    y2=getYforHist(histogram2,1)
    I2=0
    for j in range(0, len(y2)):
        I2 = I2 + y2[j]
    return I2/I1

def getFactor(f1,name1,f2,name2):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    y1=getYforHist(histogram1,1)
    I1=0
    for j in range(0, len(y1)):
        I1 = I1 + y1[j]

    y2=getYforHist(histogram2,1)
    I2=0
    for j in range(0, len(y2)):
        I2 = I2 + y2[j]
    return I2/I1

def getFactorTGraph(f1,name1,f2,name2):
    myGraph1=f1.Get(name1)
    xP1=myGraph1.GetX()
    yP1=myGraph1.GetY()

    myGraph2=f2.Get(name2)
    xP2=myGraph2.GetX()
    yP2=myGraph2.GetY()

    I1=0
    for i in range(myGraph1.GetN()):
        I1 = I1 + yP1[i]
        
    I2=0
    for i in range(myGraph2.GetN()):
        I2 = I2 + yP2[i]

    return I2/I1

def getFactor_ep(f1,f2):
    II1=0
    II2=0
    for i in range(0,7):
        name="h_Q2_ep_SRC_Q2_"+str(i)
        histogram1 = f1.Get(name)
        histogram2 = f2.Get(name)
        y1=getYforHist(histogram1,1)
        I1=0
        for j in range(0, len(y1)):
            I1 = I1 + y1[j]
            
        y2=getYforHist(histogram2,1)
        I2=0
        for j in range(0, len(y2)):
            I2 = I2 + y2[j]

        II1=I1
        II2=I2
    return II2/II1

def getFactor_epp(f1,f2):
    II1=0
    II2=0
    for i in range(0,7):
        name="h_Q2_epp_SRC_Q2_"+str(i)
        histogram1 = f1.Get(name)
        histogram2 = f2.Get(name)
        y1=getYforHist(histogram1,1)
        I1=0
        for j in range(0, len(y1)):
            I1 = I1 + y1[j]
            
        y2=getYforHist(histogram2,1)
        I2=0
        for j in range(0, len(y2)):
            I2 = I2 + y2[j]

        II1=I1
        II2=I2
    return II2/II1

def plotFromFile(f1,name,mycolor,factor):
    histogram = f1.Get(name)    
    plt.errorbar(np.array(getXforHist(histogram)),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,marker='o',linestyle='',linewidth=0,elinewidth=3,markersize=9)

def plotFromFileAx(axs,f1,name,mycolor,factor):
    histogram = f1.Get(name)    
    axs.errorbar(np.array(getXforHist(histogram)),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,marker='o',linestyle='',linewidth=0,elinewidth=3,markersize=9)

def plotFromFileLineAx(axs,f1,name,mycolor,factor):
    histogram = f1.Get(name)    
    #axs.errorbar(np.array(getXforHist(histogram)),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,marker='',linestyle='-',linewidth=1,elinewidth=0,markersize=0)
    axs.plot(np.array(getXforHist(histogram)),getYforHist(histogram,factor),color=mycolor,marker='',linestyle='-',linewidth=1,markersize=0)
    
def plotFromFileStep(f1,name,mycolor,factor):
    histogram = f1.Get(name)    
    plt.errorbar(getXforHist(histogram),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,linestyle='')
    plt.step(getXforHist(histogram),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,where='mid')

def plotFromFileLineAx(axs,f1,name,mycolor,factor):
    histogram = f1.Get(name)    
    axs.plot(np.array(getXforHist(histogram)),getYforHist(histogram,factor),color=mycolor,marker='',linestyle='-',linewidth=1,markersize=0)

def plotFromFileLineAxRb(axs,f1,name,mycolor,factor,rebin):
    histogram = f1.Get(name)    
    f1.Rebin(rebin)
    axs.plot(np.array(getXforHist(histogram)),getYforHist(histogram,factor),color=mycolor,marker='',linestyle='-',linewidth=1,markersize=0)

def plotFromFileStepAx(ax,f1,name,mycolor,factor):
    histogram = f1.Get(name)
    ax.errorbar(getXforHist(histogram),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,linestyle='')
    ax.step(getXforHist(histogram),getYforHist(histogram,factor),color=mycolor,where='mid')
    
def plotFromFileStepAxRb(ax,f1,name,mycolor,factor,rebin):
    histogram = f1.Get(name)
    f1.Rebin(rebin)
    ax.errorbar(getXforHist(histogram),getYforHist(histogram,factor),getYErrforHist(histogram,factor),color=mycolor,linestyle='')
    ax.step(getXforHist(histogram),getYforHist(histogram,factor),color=mycolor,where='mid')
    
def plotFromFileNormStep(f1,name,mycolor):
    histogram = f1.Get(name)
    y=getYforHist(histogram,1)
    I=0
    for j in range(0, len(y)):
        I = I + y[j]
    plt.errorbar(getXforHist(histogram),getYforHist(histogram,1/I),getYErrforHist(histogram,1/I),color=mycolor,linestyle='')
    plt.step(getXforHist(histogram),getYforHist(histogram,1/I),getYErrforHist(histogram,1/I),color=mycolor,where='mid')
    
def plotFromFileNoError(f1,name,mycolor,mylinestyle,factor):
    histogram = f1.Get(name)    
    plt.plot(getXforHist(histogram),getYforHist(histogram,factor),color=mycolor,linestyle=mylinestyle)

def plotFromFileRatio(f1,name1,f2,name2,mycolor,mylinestyle):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    xe=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)        
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
            xe.append(histogram1.GetBinWidth(i)/2)
        else:
            y.append(0)
            e.append(0)
            xe.append(0)
    plt.errorbar(x,y,yerr=e,xerr=xe,color=mycolor,linestyle=mylinestyle)

def plotFromFileRatioStep(ax,f1,name1,f2,name2,mycolor):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
        else:
            y.append(0)
            e.append(0)

    x.append(1000)
    y.append(y[-1])
    e.append(e[-1])
    ax.errorbar(x,y,e,color=mycolor,linestyle='')
    ax.step(x,y,color=mycolor,where='mid')

def plotFromFileRatioLine(ax,f1,name1,f2,name2,mycolor):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
        else:
            y.append(0)
            e.append(0)
    x.append(1000)
    y.append(y[-1])
    e.append(e[-1])
    ax.plot(x,y,color=mycolor,linestyle='-')
    
    
def plotFromFileRatioForOr(f1,name1,f2,name2,mycolor):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
        else:
            y.append(0)
            e.append(0)
    plt.errorbar(x,y,e,color=mycolor,linewidth=0,elinewidth=3,zorder=3,marker='o',markersize=9)
    #plt.step(x,y,color=mycolor,where='mid')
    
def ratioGetX(f1,name1,f2,name2):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    return x

def ratioGetY(f1,name1,f2,name2):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
        else:
            y.append(0)
            e.append(0)
    return y

def ratioGetE(f1,name1,f2,name2):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    e=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        e1=histogram1.GetBinError(i)
        e2=histogram2.GetBinError(i)
        if((y2!=0) and (y1!=0)):
            y.append(y1/y2)
            e.append((y1/y2)*np.sqrt(((e1/y1)**2)+((e2/y2)**2)))
        else:
            y.append(0)
            e.append(0)
    return e


###################################################
def plotTGraphAsymmErrorsFig2(myax,myFile,myName,lc,fc):
    myGraph=myFile.Get(myName)
    xP=myGraph.GetX()
    yP=myGraph.GetY()
    x=[]
    y=[]
    ylow=[]
    yhigh=[]
    for i in range(myGraph.GetN()):
        x.append(xP[i])
        y.append(yP[i])
        ylow.append(yP[i]-myGraph.GetErrorYlow(i))
        yhigh.append(yP[i]+myGraph.GetErrorYhigh(i))
    myax.plot(x,y,lc)
    myax.fill_between(x,ylow,yhigh,alpha=0.3,facecolor=fc)

def plotFromFileRatioFig2(myax,f1,name1,f2,name2,mycolor,mylinestyle):
    histogram1 = f1.Get(name1)
    histogram2 = f2.Get(name2)
    x=getXforHist(histogram1)
    y=[]
    for i in range(1,(histogram1.GetNbinsX()+1)):
        y1=histogram1.GetBinContent(i)
        y2=histogram2.GetBinContent(i)
        y.append(y1/y2)
    myax.plot(x,y,color=mycolor,linestyle=mylinestyle)


###########################################


def plotNewSig(pdf,Q2bins,f1,name,xmin,xmax,y1label,x2label):
    fig1=plt.figure()
    plt.ylabel(y1label,fontsize=15)
    plt.xlabel(r'$Q^{2} GeV$',fontsize=15)
    plt.tight_layout()
    plt.xlim(1.5,5.0)
    plt.ylim(0.00,0.34)
    print('sigma'+name+'_Q2')
    print(f1)
    plotTGraphErrorsStep(f1,'sigma'+name+'_Q2','black')
    pdf.savefig(fig1)         
    
#    Q2_15=[0.4,0.5]
#    Q2_5=[0.07,0.02]
#    fig1=plt.figure()
#    ax_back = fig1.add_axes([0.0,0.0,1.0,1.0])
#    ax_back.set_ylim(0.0,1.0)
#    ax_back.set_xlim(0.0,1.0)
#    ax_back.plot([Q2_15[0],Q2_5[0]],[Q2_15[1],Q2_5[1]],'k',linewidth=1)
#    ax_back.text(0.1,0.35,r'$Q^{2} [GeV]$',fontsize=15)
#    ax_back.text(0.05,0.9,r'$(e,e\prime pp)$',fontsize=25)
#    #ticks
#    ticks=[2.0,2.5,3.0,3.5,4.0,4.5,5.0]
#    for t in ticks:
#        t_x=Q2_15[0]+((t-1.5)/(5.0-1.5))*(Q2_5[0]-Q2_15[0])
#        t_y=Q2_15[1]+((t-1.5)/(5.0-1.5))*(Q2_5[1]-Q2_15[1])
#        ax_back.plot([t_x,t_x-0.01],[t_y,t_y],'k',linewidth=1.0)
#        ax_back.text(t_x-0.02,t_y,str(t),horizontalalignment='right',verticalalignment='center')
#    
#    for i in range(0,7):
#        Q2=(Q2bins[i]+Q2bins[i+1])/2
#        hist_name='h_'+name+'_epp_SRC_Q2_'+str(i)
#        ax_x=Q2_15[0]+((Q2-1.5)/(5.0-1.5))*(Q2_5[0]-Q2_15[0])
#        ax_y=Q2_15[1]+((Q2-1.5)/(5.0-1.5))*(Q2_5[1]-Q2_15[1])
#        ax = fig1.add_axes([ax_x,ax_y,0.7,0.5])
#        #ax = fig1.add_axes([1.0+(i/100),1.0,1,1])
#        plotFromFileStepAx(ax,f1,hist_name,'blue',1)
#        plotFromFileLineAx(ax,f2,hist_name,'red',getFactor(f2,hist_name,f1,hist_name))
#        ax.spines['top'].set_visible(False)
#        ax.spines['right'].set_visible(False)
#        if(i!=0):
#            ax.spines['left'].set_visible(False)        
#            ax.set_yticks([])
#        if(i!=6):
#            ax.set_xticks([])        
#        ax.patch.set_alpha(0.0)    
#        ax.set_ylim(0,120)
#        ax.set_xlim(xmin,xmax)
#        if(i==0):
#            ax.set_ylabel(r'Counts',fontsize=15)
#        if(i==6):
#            ax.set_xlabel(x2label,fontsize=15)
#        #ax.title(name)
#        #ax.tight_layout()
#    pdf.savefig(fig1)         
    

def plotSig(pdf,Q2bins,f1,f2,name,xmin,xmax,y1label,x2label,factor):
    
    Q2_15=[0.4,0.5]
    Q2_5=[0.07,0.02]
    fig1=plt.figure()
    ax_back = fig1.add_axes([0.0,0.0,1.0,1.0])
    ax_back.set_ylim(0.0,1.0)
    ax_back.set_xlim(0.0,1.0)
    ax_back.plot([Q2_15[0],Q2_5[0]],[Q2_15[1],Q2_5[1]],'k',linewidth=1)
    ax_back.text(0.1,0.35,r'$Q^{2} [GeV]$',fontsize=15)
    ax_back.text(0.05,0.9,r'$(e,e^{\prime}pp)$',fontsize=25)
    #ticks
    ticks=[2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    for t in ticks:
        t_x=Q2_15[0]+((t-1.5)/(5.0-1.5))*(Q2_5[0]-Q2_15[0])
        t_y=Q2_15[1]+((t-1.5)/(5.0-1.5))*(Q2_5[1]-Q2_15[1])
        ax_back.plot([t_x,t_x-0.01],[t_y,t_y],'k',linewidth=1.0)
        ax_back.text(t_x-0.02,t_y,str(t),horizontalalignment='right',verticalalignment='center')
    
    for i in range(0,7):
        Q2=(Q2bins[i]+Q2bins[i+1])/2
        hist_name='h_'+name+'_epp_SRC_Q2_'+str(i)
        ax_x=Q2_15[0]+((Q2-1.5)/(5.0-1.5))*(Q2_5[0]-Q2_15[0])
        ax_y=Q2_15[1]+((Q2-1.5)/(5.0-1.5))*(Q2_5[1]-Q2_15[1])
        ax = fig1.add_axes([ax_x,ax_y,0.7,0.5])
        #ax = fig1.add_axes([1.0+(i/100),1.0,1,1])
        plotTGEStep(ax,f1,hist_name,'black')
        plotTGELine(ax,f2,hist_name,'blue',factor[i])
        #plotFromFileStepAx(ax,f1,hist_name,'black',1)
        #plotFromFileLineAx(ax,f2,hist_name,'red',getFactor(f2,hist_name,f1,hist_name))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if(i!=0):
            ax.spines['left'].set_visible(False)        
            ax.set_yticks([])
        if(i!=6):
            ax.set_xticks([])        
        ax.patch.set_alpha(0.0)    
        ax.set_ylim(0,120)
        ax.set_xlim(xmin,xmax)
        if(i==0):
            ax.set_ylabel(r'Counts',fontsize=15)
        if(i==6):
            ax.set_xlabel(x2label,fontsize=15)
        #ax.title(name)
        #ax.tight_layout()
    pdf.savefig(fig1)         

