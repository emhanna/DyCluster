# Modified by E.M. Hanna
#!/usr/bin/python
# c 2011 Jos Berenger, Nazar Zaki and Dmitry Efimov., 2011.JUN.22
# UAE University, FIT, Bioinformatics Lab more on @bioAE Twitter

from optparse import OptionParser
from collections import defaultdict
from math import *
import string
import itertools
import operator
from numpy import *
import numpy as num
import numpy.matlib as M
from numpy.matlib import rand,zeros,ones,empty,eye
import os
import re

def sortedDictValues1(adict):
    items = adict.items()
    items.sort()
    return [value for key, value in items]

#Read proteins from file of interactions
def readFileInteractions(filename):
    relations = defaultdict(list)
    label2id = {}
    id2label = {}
    total = 0
    inp = open (filename,"r");
    for line in inp.readlines():
        p = string.split(line,"\t")
        x = p[0].strip();
        y = p[1].strip();
        if not label2id.has_key(x):
            label2id[x] = len(label2id);
            id2label[len(id2label)] = x
        if not label2id.has_key(y):
            label2id[y] = len(label2id);
            id2label[len(id2label)] = y
        if ((not label2id[x] in relations[label2id[y]]) and x != y):
            total = total + 1
            relations[label2id[y]].append(label2id[x]);
            relations[label2id[x]].append(label2id[y]);
    print "Total number of interactions (duplications are removed):", total
    print "Total number of proteins:", len(label2id)
    
      
    return relations,label2id,id2label

def findHubs(rel,id2label,th):
    hubs = []
    for id1 in id2label:
        neigh1 = set(rel[id1])
        if len(neigh1) > th:
            hubs.append(id1)
    return hubs

##### Pass file name         
def removeLowWeightedEdges(fname,rel,id2p,W,th):
    #count = 0
    for id1 in id2p:
        neigh1 = set(rel[id1])
        k = len(neigh1)
        if k == 0:
            continue
        #calculate average in weights
        sum = 0.0
        for id2 in neigh1:
            sum = sum + W[id1,id2]
        avg = sum/k
        if avg == 0.0:
            avg0 = 1.0
        else:
            avg0 = avg
        for id2 in neigh1:
            if ((W[id1,id2]-avg)/avg0 < -1.0*float(th)):
                rel[id1].remove(id2)
                #count = count + 1
                
    #print count," interactions removed";   
     
    #########write filtered PPI file
    with open("filtered1.txt",'w+') as outfile:
        for idd1 in id2p:
            neighb1 = set(rel[idd1])
            k1 = len(neighb1)
            if k1 ==0:
                continue
                
            for idd2 in neighb1:
                
                string1 = id2p[idd1]+'\t'+id2p[idd2]+'\n'
                string2 = id2p[idd2]+'\t'+id2p[idd1]#+'\n'
                
                outfile.write(string1)
                
    outfile.close()
    count2 = 0
    with open("filtered1.txt",'r') as infile:
        
        lines = infile.readlines()
        #print lines
        for item in lines:
            #print "in for loop"
            
            strng = re.split("\t|\n",item)
            #print strng
            str1 = strng[1]+'\t'+strng[0]+'\n'
            
            #print "str1 = ",str1;
            if str1 in lines:
                count2 +=1
                lines.remove(item)
                
        print count2," lines deleted";
       
        with open("PE_"+fname,'w') as outf:
            for line in lines:
                outf.write(line)
        
        

def calculateEdgeProbability(rel,id2p,hubs,N,itN):
    P = zeros((N,N))
    for id1 in id2p:
        if id1 in hubs:
            continue
        neigh1 = set(rel[id1])
        k = len(neigh1)
        if k == 0:
            continue
        for id2 in neigh1:
            P[id1,id2] = 0.5
    for i in range(0,itN):
        Pnew = zeros((N,N))
        for id1 in id2p:
            if id1 in hubs:
                continue
            neigh1 = set(rel[id1])
            for id2 in neigh1:
                if id2 in hubs:
                    continue
                if id2<=id1:
                    continue
                neigh2 = set(rel[id2])
                prob = 1.0
                for id3 in neigh1&neigh2:
                    prob = prob*(1.0-(P[id1,id3]*P[id2,id3]))
                #neigh21 = neigh2 - neigh1
                #neigh12 = neigh1 - neigh2
                #for id21 in neigh21:
                #    for id12 in neigh12:
                #        if id12 in rel[id21]:
                #            prob = prob*(1.0-(P[id1,id12]*P[id12,id21]*P[id21,id2]))
                prob = 1.0 - prob
                Pnew[id1,id2] = prob
                Pnew[id2,id1] = prob
        P = Pnew.copy()
    return P


def createComplexesMaxClusteringCoeff(pr,id2p,rel,thhub):
    k = len(pr)
    for id1 in id2p:
        neigh1 = set(rel[id1])
        if len(neigh1)<=2:
            continue
        if len(neigh1)>=thhub:
            continue
        k = k+1
        id2dl = {}
        for id2 in neigh1:
            neigh2 = set(rel[id2])
            id2dl[id2] = len(neigh2&neigh1)
        id2dls = sorted(id2dl.iteritems(), key=operator.itemgetter(1))
        l = len(id2dls)
        idseq = []
        for i in range(0,l):
            idseq.append(id2dls[i][0])
        idseq.reverse()
        cl = 0
        for id2 in neigh1:
            neigh2 = set(rel[id2])
            cl = cl + len(neigh1&neigh2)
        ccprev = cl*1.0/(l*l*(l-1))
        idchosen = id2dls[0][0]
        while l >= 3:
            id2 = id2dls[0][0]
            cl = 0
            for id3 in neigh1:
                neigh2 = set(rel[id3])&neigh1
                cl = cl + len(neigh2)
            cc = cl*1.0/(l*l*(l-1))
            if cc > ccprev:
                ccrprev = cc
                idchosen = id2
            neigh2 = set(rel[id2])&neigh1
            for id3 in neigh2:
                if id2dl.has_key(id3):
                    id2dl[id3] = id2dl[id3] - 1
                    if id2dl[id3] == 0:
                        del id2dl[id3]
            del id2dl[id2]
            neigh1 = neigh1 - set([id2])
            id2dls = sorted(id2dl.iteritems(), key=operator.itemgetter(1))
            l = len(id2dls)

        for i in range(0,len(idseq)):
            if idseq[i] == idchosen:
                pr[k].append(idseq[i])
                break
            pr[k].append(idseq[i])
        pr[k].append(id1)

def deleteSmallComplexes(pr,th_small):
    idc2len = {}
    for idc in pr.iterkeys():
        complex = set(pr[idc])
        idc2len[idc] = len(complex)
    for id in idc2len:
        if (idc2len[id] > th_small):
            continue
        try:
            print "deleted"
            del pr[id]
        except:
            continue
       

def mergeDuplicatedComplexes(pr,rel,id2p,th_merge=0.7):
    a = set(pr.iterkeys())
    b = set(pr.iterkeys())
    for c1 in a:
        for c2 in b:
            if c1 == c2:
                continue
            if not pr.has_key(c1):
                continue
            if not pr.has_key(c2):
                continue
            complex1 = set(pr[c1])
            complex2 = set(pr[c2])
            l1 = len(complex1)
            l2 = len(complex2)
            inter = complex1 & complex2
            q = (1.0*len(inter)/l1)*(1.0*len(inter)/l2)
            #q = 1.0*len(inter)/(1.0*l1) #+1.0*len(inter)/(1.0*l2)
            #print (" merging? " + str(q) + " "
            #       + str(inter) + " " + str(complex1) + " " + str(complex2))
            if q >= th_merge:
                #print(str(q))
                for id in complex2-complex1:
                    pr[c1].append(id)
                del pr[c2]

def complexEstimation(pr,rel,id2p,th_addback=0.5):
    idfordel = []
    for key1 in pr.iterkeys():
        a = pr[key1]
        inside = 0
        for id in a:
            neigh = set(rel[id])&set(a)
            inside = inside + len(neigh)
        inside = (1.0*th_addback) * (1.0 *inside)
        #print "len:",len(a),"; inside:",inside
        if inside<len(a)-1:
            idfordel.append(key1)
    for key1 in idfordel:
        del pr[key1]    

def findAttachedProteins(pr,rel):
    for key in pr.iterkeys():
        complex = set(pr[key])
        size = len(complex)
        for id in complex:
            neigh_id = set(rel[id]) - complex
            for id2 in neigh_id:
                neigh_id2 = set(rel[id2])&complex
                nneigh = len(neigh_id2)
                if (nneigh*2>size and (not id2 in pr[key])):
                    pr[key].append(id2)

def createResultFile(pr,id2protein,f = "result.txt"):
    f = open(f, 'w')
    k = 0
    for compl in pr.iterkeys():
        k = k + 1
        s = "C" + str(k) + ": "
        for prot in pr[compl]:
            s = s + " "+ id2protein[prot]
        print >> f, s
        #print s
    f.close()

def writeP(id2protein, protein2id, P, filename):
    f = open('resultP.txt', 'w')
    inp = open (filename,"r");
    for line in inp.readlines():
        p = string.split(line,"\t")
        x = p[0].strip();
        y = p[1].strip();
        prob = P[protein2id[x],protein2id[y]]
        if prob== 0:
            prob = 0.01
        s = x + "\t" + y + "\t" + str(prob)
        print >> f, s
        print s

    f.close()
    inp.close()

def readFileInteractions(filename):
    relations = defaultdict(list)
    label2id = {}
    id2label = {}
    total = 0
    inp = open (filename,"r");
    for line in inp.readlines():
        p = string.split(line,"\t")
        x = p[0].strip();
        y = p[1].strip();
        if not label2id.has_key(x):
            label2id[x] = len(label2id);
            id2label[len(id2label)] = x
        if not label2id.has_key(y):
            label2id[y] = len(label2id);
            id2label[len(id2label)] = y
        if ((not label2id[x] in relations[label2id[y]]) and x != y):
            total = total + 1
            relations[label2id[y]].append(label2id[x]);
            relations[label2id[x]].append(label2id[y]);
    print "Total number of interactions (duplications are removed):", total
    print "Total number of proteins:", len(label2id)

    return relations,label2id,id2label


def main(options):
    
    print (options)
    
    relations,protein2id,id2protein = readFileInteractions(options.fileinteractions)
    hubs = []
    pr = defaultdict(list)
    N = len(protein2id)
    hubs = findHubs(relations,id2protein,610)
    P = calculateEdgeProbability(relations,id2protein,hubs,N,1)
    removeLowWeightedEdges(options.fileinteractions,relations,id2protein,P,float(options.removeLowWeightedEdges))
    
    #createComplexesMaxClusteringCoeff(pr,id2protein,relations,float(options.minhubsize))
    #complexEstimation(pr,relations,id2protein,float(options.addback))
    #mergeDuplicatedComplexes(pr,relations,id2protein,float(options.overlap))
    #findAttachedProteins(pr,relations)
    #deleteSmallComplexes(pr,2)
    #createResultFile(pr,id2protein)
    #writeP(id2protein, protein2id, P, options.fileinteractions)
