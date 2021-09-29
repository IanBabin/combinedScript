##Pull data from all files on synology regardless of modification date
##This is a backup script in the event of major changes or loss, and should only
##Be run as a last resort. 

import mysql.connector
import os
import tarfile
import re
import time
import glob
import mmap
from fractions import Fraction
import numpy as np
import matplotlib.pyplot as plt

##Create Connection to database

sanddb = mysql.connector.connect(
//
)

keyTableFile = open("/Home/physics/babin/sandDatabaseCsv/KEYTABLE.csv","w")
energyTableFile = open("/Home/physics/babin/sandDatabaseCsv/ENERGYTABLE.csv","w")
waveFunctionFile = open("/Home/physics/babin/sandDatabaseCsv/WAVEFUNCTIONTABLE.csv","w")


##Select the current highest ID 
dbcursor = sanddb.cursor()
dbcursor.execute("SELECT MAX(ID) FROM KEYTABLE")
ID = dbcursor.fetchone()[0]
print(ID)
colors = ['xkcd:blood red','xkcd:red orange','xkcd:goldenrod','xkcd:yellow green','xkcd:bright green','xkcd:seafoam','xkcd:blue green','xkcd:sea blue','xkcd:purplish blue','xkcd:medium purple','xkcd:indigo','xkcd:aubergine']
def numberReturn(LI):
  return re.sub('[^0-9]','', LI)

## Function for determination of ground state parity based on element and mass
def parity(element,mass):
  elements=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']
  Z = elements.index(element)+1
  N = int(mass) - Z
  def nShell(n):
    if 0<n<=2:
      return pow(pow(-1,0),n)
    elif 3<=n<=6:
      return pow(pow(-1,1),(n-2))
    elif 7<=n<=12:
      return pow(pow(-1,2),(n-8))
    elif 13<=n<=20:
      return pow(pow(-1,3),(n-12))
    elif 21<=n<=28:
      return pow(pow(-1,4),(n-20))
  def pShell(n):
    if 0<n<=2:
      return pow(pow(-1,0),n)
    elif 3<=n<=6:
      return pow(pow(-1,1),(n-2))
    elif 7<=n<=12:
      return pow(pow(-1,2),(n-8))
    elif 13<=n<=20:
      return pow(pow(-1,3),(n-12))
    elif 21<=n<=28:
      return pow(pow(-1,4),(n-20))
  parity = pShell(Z)*nShell(N)
  if parity == 1:
    return 0
  else:
    return 1
##Function responsible for extraction, processing, and upload of data
##Defined to facilitate handling both .tar and .tgz format archives
def extraction(file):
  try:
    x0 = 1
    y0 = 1/np.sqrt(3)
    eigenDataFile = None
    global ID
    global colors
  ## Isolate file name and use it to determine mass, element, parity, interaction model, hw, nmax, jj, version, and other misc info
    ID = ID + 1
    fileName = os.path.splitext(file)[0].split("/")[len(os.path.splitext(file)[0].split('/'))-1]
    fileName = fileName.replace("n2lo_opt","N2LOopt")
    fileName = fileName.replace("_Vcoul_","_")
    fileName = fileName.replace("_bare_","_")
    fileName = fileName.replace("_V.cf_","_")
    fileName = fileName.replace("_V_","_")
    bgID = fileName
    mass = re.search('[\d]+', fileName.split("_")[0]).group(0)
    element = re.search('[a-zA-Z]+',fileName.split("_")[0]).group(0)
    fileName = fileName.replace(str(mass)+element,"")
    interactionModel = re.search("N2LOopt|JISP16|N3LO|NNLOopt|NNLOsat",fileName).group(0)
    fileName = fileName.replace(interactionModel,"")
    hw=re.search("hw[\d]+\.?[\d]*MeV",fileName).group(0)
    fileName = fileName.replace(hw,"")
    hw = hw.replace("hw","")
    hw = hw.replace("MeV","")
    if re.search("_v[\d]+_",fileName) is not None:
      version = re.search("_v[\d]+_",fileName).group(0)
      fileName = fileName.replace(version, "_")
    else:
      version = "0"
    version = version.strip("_")
    version = version.replace("v","")
    nMax = re.search("(_Nmax)[\d]{1,3}(_\d+){0,3}",fileName).group(0)
    fileName = fileName.replace(nMax,"")
    nMax = nMax.strip("_")
    nMax = nMax.replace("Nmax","").lstrip("0")
    nMaxArray = [i.lstrip('0') for i in nMax.split("_")]
    nMax = "_".join(nMaxArray)
    if len(nMaxArray) >= 2:
      if (int(nMaxArray[1]) % 2) != 0:
        displayNMax = int(nMaxArray[1]) - 1
      else:
        displayNMax = nMaxArray[1]
    elif len(nMaxArray) == 1:
      if (int(nMaxArray[0]) % 2) != 0:
        displayNMax = int(nMaxArray[0]) - 1
      else:
        displayNMax = nMaxArray[0]
    jj=re.search("_JJ[\d]+",fileName).group(0)
    jj = jj.strip("_")
    fileName = fileName.replace(jj,"")
    jj = float(jj.replace("JJ",""))
    jj = Fraction(jj/2)
    groundStateParity = parity(element,mass)
    fileName = fileName.strip("_")
    fileName = re.sub("__+","_",fileName)
    keyTableFile.write(str(ID)+","+str(mass)+","+element+","+interactionModel+","+str(hw)+","+nMax+","+str(displayNMax)+","+str(groundStateParity)+","+str(version)+","+fileName+",0\n")
    keyInsert = "INSERT INTO KEYTABLE (ID,MASS,ELEMENT,INTERACTIONMODEL,HW,NMAX,DISPNMAX,GSP,VERSION,MISC,PUBLIC) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    keyEntries = (ID,mass,element,interactionModel,hw,nMax,displayNMax,groundStateParity,version,fileName,0)
    dbcursor.execute(keyInsert,keyEntries)
    sanddb.commit()

  ## Excitation Energies - Find the mfdn eigenvalue data file via regular expression. Select the first 10 lines of the last entry
  ## or the last entry if it's less than 10 lines
    tar = tarfile.open(file)
    for member in tar.getmembers():
      if re.search("^.+(mfdn_e).+",member.name):
        eigenDataFile = member.name
    if eigenDataFile is not None:
      f = tar.extractfile(eigenDataFile)
      temp = []
      output = []
      for line in f.readlines():
        temp.append(line.split())
        size = len(temp)
      lastEntrySize = int(temp[size-1][0].decode("utf-8"))
      if size <= 10:
        start = size - lastEntrySize
        end = start + lastEntrySize
        for i in range(start,end):
          energyEntry = (ID,str(jj),float(temp[i][1].decode("utf-8")))
          output.append(energyEntry)
          energyTableFile.write(str(ID)+","+str(jj)+","+str(temp[i][1].decode("utf-8"))+"\n")
      if size > 10:
        start = size - lastEntrySize
        end = start + 10
        for i in range(start,end):
          energyEntry = (ID,str(jj),float(temp[i][1].decode("utf-8")))
          output.append(energyEntry)
          energyTableFile.write(str(ID)+","+str(jj)+","+str(temp[i][1].decode("utf-8"))+"\n")
      energyInsert = "INSERT INTO ENERGYTABLE (ID,SPIN,ENERGY) VALUES (%s,%s,%s)"
      dbcursor.executemany(energyInsert,output)
      sanddb.commit()

    ##Wave Function data 
    for member in tar.getmembers():
        if re.match("^.*eigenvector01.log",member.name):
            eigenVectorLogFile = member.name
            break
        else:
            eigenVectorLogFile = None
    if eigenVectorLogFile != None:
        f = tar.extractfile(eigenVectorLogFile)
        outfile = open("/Home/physics/babin/temp/tempWaveFunc.log","w+")
        for line in f.readlines():
            outfile.write(line.decode("utf-8"))
        f.close()
        outfile.close()
        with open("/Home/physics/babin/temp/tempWaveFunc.log","r+",encoding="utf-8") as inFile:
          fig = plt.figure(figsize=(20,10))
          ax = fig.add_subplot(polar = True)
          ax.set_ylim(-.1,1.2)
          ax.set_thetamin(0)
          ax.set_thetamax(60)
          ax.set_title(str(mass)+""+element+", Spin = "+str(jj), va='bottom')
          ax.set_xlabel("Beta")
          ax.yaxis.set_label_position("right")
          ax.set_ylabel("Gamma")
          lines = []
          infoLine = []
          waveEntry = []
          findNmax = '^Nmax = \d+$'
          for line in inFile:
            if re.search('^Nmax = \d+$',line):
                lines.append(int(re.search(findNmax,line).group(0).replace("Nmax = ","")))
            if re.search('^(Number of protons: [\d]+	 Number of neutrons:[\d]+	Jcut=[\d]+)$',line):
                infoLine = line.split("\t")
          proton = float(re.search(r'\d+', infoLine[0]).group())
          neutron = float(re.search(r'\d+', infoLine[1]).group())
          jcut = int(re.search(r'\d', infoLine[2]).group())
          if 0<proton<=2:
              nProton = (3/2)*proton
          elif 3<=proton<=8:
              nProton = 3+(5/2)*(proton-2)
          elif 9<=proton<=20:
              nProton = 18+(7/2)*(proton-8)
          elif 21<=proton<=28:
              nProton = 60+(9/2)*(proton-20)    
          if 0<neutron<=2:
              nNeutron = (3/2)*neutron
          elif 3<=proton<=8:
              nNeutron = 3+(5/2)*(neutron-2)
          elif 9<=proton<=20:
              nNeutron = 18+(7/2)*(neutron-8)
          elif 21<=proton<=28:
              nNeutron = 60+(9/2)*(neutron-20)
          nSigma = float(nNeutron + nProton - (3/2))
          data = mmap.mmap(inFile.fileno(),0)
          for nSet in range(len(lines)):
            findData = '(?<=Nmax = '+str(lines[nSet])+'\n)(?s)(.*)(?=Saving into Nmax'+str(lines[nSet])+'\.table)'
            findData = findData.encode("utf-8")
            returnMatch = re.search(findData,data,re.MULTILINE).group().decode("utf-8").strip()
            returnMatch = returnMatch.split("\n")
            for i in range(len(returnMatch)):
                returnMatch[i] = returnMatch[i].replace("\t"," ")
            for line in range(len(returnMatch)):
                currentLine = re.sub("S=\d+",'\g<0> ',returnMatch[line])
                deformation = re.search("\(\d+ \d+\)",currentLine).group(0).replace("(","").replace(")","")
                currentLine = re.sub("\(\d+ \d+\)",deformation,currentLine).split()
                for i in range(3):
                  currentLine[i] = numberReturn(currentLine[i])
                if(float(currentLine[6])>1e-4):
                  N = float(currentLine[2])
                  lmbda = float(currentLine[3])
                  mu = float(currentLine[4])
                  k = (np.sqrt(5/np.pi)*(nSigma+N))/3
                  area = 2000*np.pi*float(currentLine[6])
                  x = (2*lmbda+mu+3)/3
                  y = (mu+1)/np.sqrt(3)
                  xRel = x - x0
                  yRel = y - y0
                  beta = np.sqrt((pow(xRel,2)/pow(k,2))+(pow(yRel,2)/pow(k,2)))
                  gamma = np.arctan(yRel/xRel)
                  ax.scatter(gamma,beta,s=area,color=colors[int(currentLine[2])],marker="o",alpha = .5)
                  entry = [ID] + currentLine
                  waveEntry.append(entry)
                  waveFunctionFile.write(str(ID)+","+str(currentLine[0])+","+str(currentLine[1])+","+str(currentLine[2])+","+str(currentLine[3])+","+str(currentLine[4])+","+str(currentLine[5])+","+str(currentLine[6])+"\n")
          keyInsert = "INSERT INTO WAVEFUNCTION (ID,PSPIN,NSPIN,TOTALSPIN,LAMBDA,MU,BASISSTATES,PROBABILITY) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
          dbcursor.executemany(keyInsert,waveEntry)
          sanddb.commit()
          plt.savefig("/Home/physics/babin/betaGammaOutput/"+str(ID)+".jpg")
          plt.close()
  except Exception as e:
    print("Error in processing: "+file)
    print("Exception: "+str(e))
i = 1
## Change to main SA-NCSM directory on synology 
os.chdir("/synology/Common/SA-NCSM")
fileCount = len(glob.glob('**/*.tgz', recursive=True))+len(glob.glob('**/*.tar', recursive=True))
if ID is None:
  ID = 0
now = time.time()
## Recursively select all tgz format files below SA-NCSM
for file in glob.glob('**/*.tgz', recursive=True):
  modified = os.path.getmtime(file)
 # if modified > now - 2629743:
  print(str(i)+"/"+str(fileCount))
  extraction(file)
  i=i+1
for file in glob.glob('**/*.tar', recursive=True):
  modified = os.path.getmtime(file)
 # if modified > now - 2629743:
  print(str(i)+"/"+str(fileCount))
  extraction(file)
  i=i+1
keyTableFile.close()
energyTableFile.close()
waveFunctionFile.close()