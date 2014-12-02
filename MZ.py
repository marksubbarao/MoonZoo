# import packages
from numpy import *
from scipy.cluster.vq import vq,kmeans,whiten,kmeans2
import Image, ImageDraw
#import urllib
#import cStringIO


nC=0
tolerence = 10.0
w_diam = 1.0
craterList=[]
fieldCraterList=[]
fieldClusterList=[]
strayList=[]
solidCraters=[]
targetNum=0


# define Crater Class
class Crater:
  def __init__(self,class_id,annotation_id,asset_id,x,x_diam,y,y_diam,boulderyness,iden):
    self.class_id=class_id
    self.annotation_id=annotation_id
    self.asset_id=asset_id
    self.x=x
    self.x_diam=x_diam
    self.y=y
    self.y_diam=y_diam
    self.boulderyness=boulderyness
    self.iden=iden
  def circularize(self):
    self.diam=(self.x_diam+self.y_diam)/2.

class Cluster:
  def __init__(self,asset_id):
    self.asset_id=asset_id
    self.x=0
    self.x_err=0.
    self.x_diam=0
    self.x_diam_err=0.
    self.y=0
    self.y_err=0.
    self.y_diam=0
    self.y_diam_err=0.
    self.iden=0
    self.nClass=0
    self.nMarked=0.
    self.links=[]
    self.iClass=0
  def circularize(self):
    self.diam=(self.x_diam+self.y_diam)/2.


# Read in data from ascii file (sql header and footers stripped)
def parseFile(fileName):
  f = open(fileName, mode='r')
  llines = []
  nobs=0
  #import crater
  #Crater MZ
  for line in f:
    line.strip();
    if line.startswith('|'):
      line.replace("|"," ")
      line.lstrip("|")
      words = line.split("|")
      if (len(words)>2):
        crater = Crater(words[1],words[2],-999,-999,-999,-999,-999,-999)
        if(words[3].strip()!="No craters"):
          craterList.append(crater)
          nobs = nobs+1
    else:
      words = line.split("\"");
      if (words[0].strip()=="x:"):
        crater.x=float(words[1])
      elif (words[0].strip()=="x_diameter:"):
        crater.x_diam=float(words[1])
      elif (words[0].strip()=="y:"):
        crater.y=float(words[1])
      elif (words[0].strip()=="y_diameter:"):
        crater.y_diam=float(words[1])
      elif (words[0].strip()=="id:"):
        crater.iden=int(words[1])
      elif (words[0].strip()=="angle:"):
        crater.angle=float(words[1])
      elif (words[0].strip()=="boulderyness:"):
        crater.boulderyness=float(words[1])
  print nobs
  
# Read in data from ascii csv file (sql header and footers stripped)
def parseFile2(fileName):
  f = open(fileName, mode='r')
  llines = []
  nobs=0
  #import crater
  #Crater MZ
  for nline,line in enumerate(f):
    if (line[0]!="*" and line[0]!="c"):
      sepline = line.split(',')
      if ((float(sepline[8])+float(sepline[7]))>40.001):
        crater=Crater(int(sepline[0]),int(sepline[1]),int(sepline[2]),float(sepline[5]),float(sepline[7]),float(sepline[6]),float(sepline[8]),float(sepline[9]),int(sepline[3]))
        craterList.append(crater)
        if (nline % 1000 == 0):
          print("linesread = ",nline)
        nobs=nobs+1
  #sort craterList by assetID
  print "sorting..."
  craterList.sort(key=lambda x: x.asset_id)
  print "sort done"

def fieldCraters(imageID):
  global nC
  for n in range(nC,len(craterList)):
    crater=craterList[n]
    if (imageID != crater.asset_id):
      break
    if (int(imageID)==int(crater.asset_id)):
      crater.circularize()
      fieldCraterList.append(crater)
  nC=n
 
def findClusters(tol):
  #creates lists of linked clusters
  nCraters=len(fieldCraterList)
  tminList=zeros(nCraters)
  for i in range(nCraters):
    fieldClusterList.append(Cluster(fieldCraterList[i].asset_id))
    tminList[i]=1000.
  for i in range(nCraters-1):
      craterA=fieldCraterList[i];
      clusterA=Cluster(craterA.asset_id)
      tmin=10.*tol
      for j in range(i,nCraters):
        if (i!=j):
          craterB=fieldCraterList[j]
          delta_x = craterA.x - craterB.x
          delta_y = craterA.y - craterB.y
          delta_d = craterA.diam - craterB.diam
  # dist is the linking length - all the magic is here.
          rad  = (craterA.diam +craterB.diam)/2.
          ws =1. - 0.01*(rad-20.0)
          ws=max(ws,0.25)
          dist = ws*sqrt(delta_x**2 + delta_y**2 + w_diam * delta_d**2)
          if (dist < tol):
            fieldClusterList[i].links.append(j)
            fieldClusterList[j].links.append(i)
          tminList[i]=min(tminList[i],dist)
          tminList[j]=min(tminList[j],dist)
  for i,crater in enumerate(fieldCraterList):
    if (tminList[i]<tol):
      solidCraters.append(crater)
    else:
      strayList.append(crater)

def consolidateClusters():
  #does friends of friends usig the lists generated in findClusters()
  numC=len(fieldClusterList)
  for i in range(numC):
    for index in fieldClusterList[i].links:
      if (i!=index):
        for indexA in fieldClusterList[index].links:
          if (fieldClusterList[i].links.count(indexA)==0):
            fieldClusterList[i].links.append(indexA)
        fieldClusterList[index].links=[]
  for i in range(numC):
    if (len(fieldClusterList[numC-i-1].links)==0):
      fieldClusterList.pop(numC-i-1)
  numC=len(fieldClusterList)
#  for i in range(numC):
#     print fieldClusterList[i].links

def calcClusters(outfile):
# generates output from cluster lists
  classifiers=[]
  for cNum,cluster in enumerate(fieldClusterList):
    cluster.nMarked=len(cluster.links)
    x_array = zeros(cluster.nMarked)
    y_array = zeros(cluster.nMarked)
    xdiam_array = zeros(cluster.nMarked)
    ydiam_array = zeros(cluster.nMarked)
    classifiers[:]=[]
    for i in range(cluster.nMarked):
      x_array[i]=fieldCraterList[cluster.links[i]].x
      y_array[i]=fieldCraterList[cluster.links[i]].y
      xdiam_array[i]=fieldCraterList[cluster.links[i]].x_diam
      ydiam_array[i]=fieldCraterList[cluster.links[i]].y_diam
      if (not(fieldCraterList[cluster.links[i]].class_id in classifiers)):
        classifiers.append(fieldCraterList[cluster.links[i]].class_id)
    #print len(classifiers)
    cluster.x=average(x_array)
    cluster.y=average(y_array)
    cluster.x_diam=average(xdiam_array)
    cluster.y_diam=average(ydiam_array)
    cluster.x_err=std(x_array)
    cluster.y_err=std(y_array)
    cluster.x_diam_err=std(xdiam_array)
    cluster.y_diam_err=std(ydiam_array)
    cluster.circularize()
    cluster.nClass=len(classifiers)
    if (len(classifiers) > 1):
      outfile.write("%d,%d,%d,%d,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n" % (cluster.asset_id,cNum,cluster.nMarked,cluster.nClass,cluster.x,cluster.x_err,cluster.y,cluster.y_err,cluster.x_diam,cluster.x_diam_err,cluster.y_diam,cluster.y_diam_err))

    
def calcTargetNum():
  classifiers=[]
  for crater in solidCraters:
    if (not(crater.class_id in classifiers)):
       classifiers.append(crater.class_id)
  targetNum=int(round(1.0*len(solidCraters)/len(classifiers)))
   
  

  
def plotCraters(asset):
  im=Image.new("RGB",(600,400))
  draw = ImageDraw.Draw(im)
  print len(fieldCraterList), len(fieldClusterList), len(solidCraters), len(strayList)
  for crater in solidCraters:
      draw.ellipse((int(crater.x-crater.diam/2.),int(crater.y-crater.diam/2.),int(crater.x+crater.diam/2.),int(crater.y+crater.diam/2.)),outline=(50,155,50,128))
  for crater in strayList:
      draw.ellipse((int(crater.x-crater.diam/2.),int(crater.y-crater.diam/2.),int(crater.x+crater.diam/2.),int(crater.y+crater.diam/2.)),outline=(255,50,50,128))
  for crater in fieldClusterList:
    if crater.nClass > 1:
      draw.ellipse((int(crater.x-crater.diam/2.),int(crater.y-crater.diam/2.),int(crater.x+crater.diam/2.),int(crater.y+crater.diam/2.)),outline=(50,50,250,128))
    else:
      draw.ellipse((int(crater.x-crater.diam/2.),int(crater.y-crater.diam/2.),int(crater.x+crater.diam/2.),int(crater.y+crater.diam/2.)),outline=(250,50,250,128))   #  if len(fieldCraterList) > 0: 
#  im.save("test_images/" + str(fieldCraterList[0].asset_id) + ".png" , "PNG")
  im.save("test_images/" + str(asset) + ".png" , "PNG")


def makeCatalog(inFilename,outFilename):
  outFile = open(outFilename, 'a')
  outFile.write("asset_id,craterNumber,nMarked,nClass,x_ave,x_err,y_ave,y_err,x_diam_ave,x_diam_err,y_diam_ave,y_diam_err\n")
  craterList[:]=[]
  parseFile2(inFilename)
  print("Number of classifications = ",len(craterList))
  #Get asset_list
  asset_list=[]
  for crater in craterList:
    asset_list.append(crater.asset_id)
  #remove duplicates by setting list to set and back again
  asset_list = list(set(asset_list))
  asset_list.sort()
  print("Number of Assets = ",len(asset_list))
  print asset_list
  for nAsset,asset in enumerate(asset_list):
    fieldCraters(asset)
    findClusters(tolerence)
    consolidateClusters()
    calcClusters(outFile)
    #plotCraters(asset)
    #reset lists
    fieldCraterList[:]=[]
    strayList[:]=[]
    solidCraters[:]=[]
    fieldClusterList[:]=[]
    if (nAsset % 1000 == 0):
      print("Assets done = ",nAsset)
    
def runOvernight():
  for i in range(0,17):
    global nC
    nC=0
    makeCatalog("class%02d.csv"%(i),"MZclusters_2_20p.csv")

def Kprocess(imageID,kC,Rcutoff):
#  file = urllib.urlopen('http://moonzoo.s3.amazonaws.com/v10/slices/000140322.jpg')
#  im = cStringIO.StringIO(file.read()) # constructs a StringIO holding the image
#  im = Image.open("000140322.jpg")
  im=Image.new("RGB",(600,400))
#  im=Image.open("http://moonzoo.s3.amazonaws.com/v10/slices/000116030.jpg")
  draw = ImageDraw.Draw(im)
  index=0
  for crater in craterList:
    crater.circularize()
    if (int(imageID)==int(crater.asset_id) and crater.diam > float(Rcutoff)):
      index=index+1
      draw.ellipse((int(crater.x-crater.diam),int(crater.y-crater.diam),int(crater.x+crater.diam),int(crater.y+crater.diam)),outline=(50,155,50,128))
  print index
  cArray =  zeros((index,3))
  index=0
  for crater in craterList:
    if (int(imageID)==int(crater.asset_id) and float(crater.diam) > Rcutoff):
      cArray[(index,0)]=crater.x
      cArray[(index,1)]=crater.y
      cArray[(index,2)]=crater.diam
      index=index+1
  results = kmeans(cArray,kC)
  codebook = results[0]
  print codebook
  print results[1]
  for i in range(kC):
    draw.ellipse((codebook[i,0]-codebook[i,2],codebook[i,1]-codebook[i,2],codebook[i,0]+codebook[i,2],codebook[i,1]+codebook[i,2]),outline=(255,128,0,255))   
  im.show()
#  im.save("test.jpg", "JPEG")
#parseFile('classificationsoftop100_2.txt')
#makeArray(297235)




