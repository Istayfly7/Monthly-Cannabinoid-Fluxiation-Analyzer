import sys
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def Err(code = 0):
    if(code == 0):
        print("Error Misc.")
        return -1
    elif(code == 1):
        print("Error: Incorrect Number of Arguments")
        return -1
    else:
        print("Unknown Error")
        return -1
       
#reads file and turns data into a numpy array
def readFile(file):
    df = pd.read_excel(file)
    df = df.fillna(0)
    #df = df.replace({'%':''})
    data = df.to_numpy()
    return data

def cutoffSample(data):
    flag = True
    fnanrow = 0

    for i in range(data.shape[0]):
        if((data[i, 0] == 0) & (flag == True)):
            fnanrow = i
            flag = False

    #print("First row of nan is: ", fnanrow)
    return data[:fnanrow, :]

#returns a dictionary of the sample profile -> includes the different types like clean, iso, oleo 
def getSampleTypes(data):
    #dictionary to store all the cannabinoids of 1 sample
    sample = dict()
    clean = -1
    oleo = -1
    iso = -1
    wax = -1
    dist = -1
    spent = -1

    #save the cannabinoids in my order in order to do comparisons - clean, oleo, iso, wax, dist, spent
    #wax, oleo - hom and degass
    for i in range(0, len(data)):
        if("CLEAN" in ((data[i])[0]).upper()):
            #save i value. data[i] = clean bio data
            clean = i
        elif("OLEO" in ((data[i])[0]).upper()):
            ##save i values. data[i] = oleoresin
            oleo = i
        elif("ISO" in ((data[i])[0]).upper()):
            ##save i values. data[i] = oleoresin
            iso = i
        elif("WAX" in ((data[i])[0]).upper()):
            ##save i values. data[i] = oleoresin
            wax = i
        elif("DIST" in ((data[i])[0]).upper()):
            ##save i values. data[i] = oleoresin
            dist = i
        elif("SPENT" in ((data[i])[0]).upper()):
            ##save i values. data[i] = oleoresin
            spent = i

    #adding stages of product to sample dictionary
    if(clean != -1):    
        sample[(data[clean])[0]] = (data[clean])[1:]
    if(oleo != -1):
        sample[(data[oleo])[0]] = (data[oleo])[1:]
    if(iso != -1):
        sample[(data[iso])[0]] = (data[iso])[1:]
    if(wax != -1):
        sample[(data[wax])[0]] = (data[wax])[1:]
    if(dist != -1):    
        sample[(data[dist])[0]] = (data[dist])[1:]
    if(spent != -1):    
        sample[(data[spent])[0]] = (data[spent])[1:]

    return sample

#find actual key name using short description-> finds the oleo type using oleo and isolate type using iso
def findType(sample, type1, type2):
    for t in sample.keys():
        if(type1 in t.lower()):
            ty1 = t
        elif(type2 in t.lower()):
            ty2 = t

    return ty1, ty2

#return list of compared values
def compareTypes(list2, list1):
    fin = np.array([])

    #remove % symbols first
    list1 = ((str(list1)[1:-1]).replace('\n', ' ')).split()
    list1 = (str(list1[:-1])[2:-2]).split('\', \"\'%\'\", \'')
    list2 = ((str(list2)[1:-1]).replace('\n', ' ')).split()
    list2 = (str(list2[:-1])[2:-2]).split('\', \"\'%\'\", \'')

    #now do work
    length = len(list1)
    for i in range(length):
        fin = np.append(fin, float(list2[i]) - float(list1[i]))
    return fin

#checking what values i have for sample type per sample
def checkAvailTypes(sample):
    #[clean, oleo, iso, wax, dist, spent] -> 1 means got values for that type*******account for nonhom to hom
    stypelist = [0, 0, 0, 0, 0, 0]
    for sampType in sample.keys():
        if('clean' in sampType.lower()):
            stypelist[0] = 1
        elif('oleo' in sampType.lower()):
            stypelist[1] = 1
        elif('iso' in sampType.lower()):
            stypelist[2] = 1
        elif('wax' in sampType.lower()):
            stypelist[3] = 1
        elif('dist' in sampType.lower()):
            stypelist[4] = 1
        elif('spent' in sampType.lower()):
            stypelist[5] = 1

    return stypelist

#average out the flux of the cannabinoids of each sample of same sample type
def averageSamples(allSampleComps):
    averagedSamples = dict()
    cocount, cwcount, sccount, cicount, cdcount = 0, 0, 0, 0, 0

    for i in range(len(allSampleComps)):
        if(allSampleComps[i]):
            for k in allSampleComps[i].keys():
                if(k == 'clean-oleo'):
                    if(cocount > 0):
                        for j in range(len(averagedSamples[k])):
                            (averagedSamples[k])[j] = (averagedSamples[k])[j] + ((allSampleComps[i])[k])[j]
                    else:
                        averagedSamples[k] = (allSampleComps[i])[k]
                    cocount += 1
                elif(k == 'clean-iso'):
                    if(cicount > 0):
                        for j in range(len(averagedSamples[k])):
                            (averagedSamples[k])[j] = (averagedSamples[k])[j] + ((allSampleComps[i])[k])[j]
                    else:
                        averagedSamples[k] = (allSampleComps[i])[k]
                    cicount += 1
                elif(k == 'clean-wax'):
                    if(cwcount > 0):
                        for j in range(len(averagedSamples[k])):
                            (averagedSamples[k])[j] = (averagedSamples[k])[j] + ((allSampleComps[i])[k])[j]
                    else:
                        averagedSamples[k] = (allSampleComps[i])[k]
                    cwcount += 1
                elif(k == 'clean-dist'):
                    if(cdcount > 0):
                        for j in range(len(averagedSamples[k])):
                            (averagedSamples[k])[j] = (averagedSamples[k])[j] + ((allSampleComps[i])[k])[j]
                    else:
                        averagedSamples[k] = (allSampleComps[i])[k]
                    cdcount += 1
                elif(k == 'spent-clean'):
                    if(sccount > 0):
                        for j in range(len(averagedSamples[k])):
                            (averagedSamples[k])[j] = (averagedSamples[k])[j] + ((allSampleComps[i])[k])[j]
                    else:
                        averagedSamples[k] = (allSampleComps[i])[k]
                    sccount += 1

    for key in averagedSamples.keys():
        for j in range(len(averagedSamples[key])):
            if(key == 'clean-oleo'):
                (averagedSamples[key])[j] = (averagedSamples[key])[j] / cocount
            elif(key == 'clean-iso'):
                (averagedSamples[key])[j] = (averagedSamples[key])[j] / cicount
            elif(key == 'clean-wax'):
                (averagedSamples[key])[j] = (averagedSamples[key])[j] / cwcount
            elif(key == 'clean-dist'):
                (averagedSamples[key])[j] = (averagedSamples[key])[j] / cdcount
            elif(key == 'spent-clean'):
                (averagedSamples[key])[j] = (averagedSamples[key])[j] / sccount

    return averagedSamples

def plot(allAveDict):
    canlist = ['THCa','Delta 9 THC','CBDa','CBD','Delta 8 THC','CBNa','CBN','CBGa','CBG','THCVa','THCV','CBDVa','CBDV','CBCa','CBC','Total Cannabinoids','Total Potential THC','Total Potential CBD']
    print("Averaged Data")
    for key in allAveDict.keys():
        plt.plot(canlist, allAveDict[key], label=key)
        print(key, ":", allAveDict[key])

    plt.ylabel('Average Potency Fluxiation')
    plt.xlabel('Stages')
    plt.xticks(rotation=50)
    plt.yticks(np.arange(0, 110, 10))
    plt.axes().yaxis.set_minor_locator(MultipleLocator(2))
    plt.legend()
    plt.grid(linestyle='--')
    plt.show()

def createReport(allAveDict):
    file = open("MonthlyReport.txt", "w")
    canlist = 'THCa    Delta 9 THC   CBDa   CBD   Delta 8 THC   CBNa   CBN   CBGa   CBG   THCVa   THCV   CBDVa   CBDV   CBCa   CBC   Total Cannabinoids    Total Potential THC   Total Potential CBD'
    file.write(canlist + '\n')

    for key in allAveDict.keys():
        #print(key, ":", allAveDict[key])
        file.write(str(key) + " : " + str(allAveDict[key]) + '\n')

    file.close()

if __name__ == "__main__":
    #get folder name of the files to be analyzed
    if(len(sys.argv) == 2):
        folder = sys.argv[1]
    else:
        Err(1)

    #get file names of files in the folder
    infolder = [f for f in listdir(folder) if isfile(join(folder, f))]

    #list of all the samples/files dictionaries
    allSamples = np.array([])
    
    #read the 1 file at a time in a loop. Read the file ananlyze and save data until no more files
    for file in infolder:
        lfile = "C:/Users/Rico-Porter/workspace/Analyzers/MonthlyCanFlux/" + folder + "/" + file
        data = readFile(lfile)
        data = cutoffSample(data[7:, 2:])
        #print("data for " + file + ":")
        #print(data)
        #print()

        #adding the sample dictionary to list of sample dictionaries*******account for nonhom to hom
        allSamples = np.append(allSamples,getSampleTypes(data))
    #print("All the picked samples:", allSamples)

    #list of all the samples/files dictionaries with comparisons
    allSampleComps = np.array([])

    for sample in allSamples:
        #dictionary to store all the cannabinoid fluxiations from step to step of 1 sample
        fluxComp = dict()
        stypelist = checkAvailTypes(sample)
        #print("stypelist:", stypelist)

        #compare the flux of each cannabinoid from one step to another(compare types)  - clean:oleo, clean:iso, clean:wax, clean:dist,  dist*******account for nonhom to hom and spents
        if(stypelist[0] == 1 and stypelist[1] == 1):
            c, o = findType(sample, 'clean', 'oleo')
            fluxComp['clean-oleo'] = compareTypes(sample[o], sample[c])
        if(stypelist[0] == 1 and stypelist[2] == 1):
            c, i = findType(sample, 'clean', 'iso')
            fluxComp['clean-iso'] = compareTypes(sample[i], sample[c])
        if(stypelist[0] == 1 and stypelist[3] == 1):
            c, w = findType(sample, 'clean', 'wax')
            fluxComp['clean-wax'] = compareTypes(sample[w], sample[c])
        if(stypelist[0] == 1 and stypelist[4] == 1):
            c, d = findType(sample, 'clean', 'dist')
            fluxComp['clean-dist'] = compareTypes(sample[d], sample[c])
        if(stypelist[0] == 1 and stypelist[5] == 1):
            c, s = findType(sample, 'clean', 'spent')
            #print("sample[c]:", sample[c])
            #print("sample[s]:", sample[s])
            fluxComp['spent-clean'] = compareTypes(sample[c], sample[s])

        allSampleComps = np.append(allSampleComps, fluxComp) 
    #print("Fluxiation Data:", allSampleComps)

    #average the fluxiation values over all the samples
    allAveDict = averageSamples(allSampleComps) 

    #plot the values
    plot(allAveDict)

    #create report**********************TO DO
    createReport(allAveDict)
