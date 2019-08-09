# -*- coding: utf-8 -*-


# Imports
from __future__ import unicode_literals
import os
import glob
import re
import argparse
import numpy as np
#from matplotlib import pyplot
import math
import itertools
import copy
from schrodinger import structure
from schrodinger.application.bioluminate import protein
import datetime
import time
import fnmatch
from schrodinger.structutils import minimize
from operator import itemgetter
parser=argparse.ArgumentParser()
args=parser.parse_args()


# Script information
__author__ = "Marc Domingo"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Marc Domingo"
__email__ = "marc.domingo@e-campus.uab.cat"


#Functions
def filter(input):

    # Pick alpha carbon condition
    for line in input:
        if '  CA ' in line and 'ATOM' in line:
            output.write(line)

def list_of_coordinates(input):
    print ("extracting coordinates")
    geometry=[]
    #Pick the coordinates from pdb
    for line in input:
        geometry.append(line[31:38])   #x_coordinate
        geometry.append(line[39:46])   #y coordinate
        geometry.append(line[47:54])   #z coordinate
        geometry.append(line[23:26])   #residue number

    #Transform to float
    for i in range(len(geometry)):
        #if i==0 or (i+1)%4 != 0:
        geometry[i]=float(geometry[i])


    #Split the list in list of list
    geometry2=[geometry[i*4:(i+1)*4] for i in range((len(geometry)+4-1)//4)]

    return geometry2

def combination(geometry2, n):
    print("Calculating all combinations")

    # Make groups of the number of residues coordinating the metal and that we want to mimic
    combinated_list=list(itertools.combinations(geometry2, n))
    return combinated_list


def zerolistmaker(n):
    listofzeros=[]
    three_zeros = [0.000]*3
    for i in range(n):
        listofzeros.append(three_zeros)

    return listofzeros


def centroid(combinated,n):
    print("Calculating centroid")
    zero_list = []
    three_zeros = [0.000] * 3
    for i in range(len(combinated)):
        zero_list.append(three_zeros)
    centroid = []
    print(combinated[0])
    np.asarray(combinated)
    np.asarray(centroid)

    # Calculate centrioid
    #print(combinated[0])
    for i in range(len(combinated)):
        centroid.append(np.mean(combinated[i],axis=0))
    #print centroid[0]
    uu=[]
    for i in range(len(centroid)):
        uu.append(np.array(centroid[i]).tolist())
        #print (uu[i])
        del uu[i][3] #borrar la mitjana
    #for i in range(len(centroid)):
        for j in range(0,4):
            uu[i].append(combinated[i][j][3])

    #print (uu[0])

    return uu

    #centroid[12230][0][0] =  combinated[1][1][1]
    #for i in range(len(centroid)):
        #del centroid[i][1:]


def module(combination_list,centroid_list,n):
    print("calculating module")

    module_list = [ [] for i in range(len(combination_list)) ]
    #print(centroid_list)
    #Calculate the module from each vertix of the prism to the centroid
    for i in range(len(combination_list)):
        for j in range(0,2*n):
            if j < n:
                module_list[i].append(math.sqrt(pow(combination_list[i][j][0]-centroid_list[i][0],2)+pow(combination_list[i][j][1]-centroid_list[i][1],2)+pow(combination_list[i][j][2]-centroid_list[i][2],2)))
            else:
                module_list[i].append(centroid_list[i][j-1])  #j-1 amb 4
    #print(module_list)
    return module_list
def score(module_list):
    print("Calculating score")

    score_function= [ [] for i in range(len(module_list)) ]
    #Calculate the score of each possible combination of alpha carbons
    
    x=3.112877226136617
    y=4.473735645408656
    z=4.7342913012403445
    t=7.851186033014375
        
    for i in range(len(module_list)):
        score_function[i].append(min(math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- t,2)),math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- z,2)),math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- t,2)), math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- t,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- t,2)),math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-x,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- z,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- t,2)),math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- t,2)), math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - z,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- z,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-y,2) + pow(module_list[i][1] - t,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- z,2)),math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- x,2)),math.sqrt(pow(module_list[i][0]-z,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- t,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- z,2)+pow(module_list[i][3]- y,2)),math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - x,2) + pow(module_list[i][2]- y,2)+pow(module_list[i][3]- z,2)),math.sqrt(pow(module_list[i][0]-t,2) + pow(module_list[i][1] - y,2) + pow(module_list[i][2]- x,2)+pow(module_list[i][3]- z,2))))  #Aquests numeros son les distancies entre els carbonis alfa i centroide del metal binding site que vull mimetitzar
        score_function[i].append(int(module_list[i][4]))
        score_function[i].append(int(module_list[i][5]))
        score_function[i].append(int(module_list[i][6]))
        score_function[i].append(int(module_list[i][7]))


    sorted_score=sorted(score_function, key=itemgetter(0))
    #output.write("%s" %len(sorted_score))

    #prova = [0.022,2,3,4]
    #["A:{}".format(prova[i]) for i in range(1,4)]


    formated_function= [ [] for i in range(len(module_list)) ]

    for i in range(len(sorted_score)):
        formated_function[i].append(score_function[i][0])
        for j in range(1, 5):# 4 changed to 5
            formated_function[i].append('A:{}'.format(score_function[i][j]))



        #["A:{}".format(prova[i][j]) for j in range(1, 4)]

    for item in sorted_score:
        output.write("%s\n" %item)
    print("End")
    #return formated_function, sorted_score


def mutate(infile,outfile):
    uu=[]
    for line in infile:
        uu.append(line)
    sorted_score = [l.split(', ') for l in uu]
    outfile.write(sorted_score[0][0])
    
    sorted_sc = sorted(sorted_score, key=itemgetter(0))
    print(sorted_sc[0])
    i=0
    for i in range(len(sorted_score)):
        list_mutations=[]
        for j in range(1,5):#4 changed to 5

            # Get the input structure
            if j==1:
                reference_st = structure.StructureReader('PELE_PaDaI.mae').next()
            if j>1:
                reference_st = structure.StructureReader('mutated_structures.mae').next()
            # Create the writer for the output file and append the reference

            writer = structure.StructureWriter('mutated_structures.mae')
            # writer.append(reference_st)

            list_mutations.append(str(sorted_sc[i][j]))
            # Define the residues and mutations

            residues = ["'A:"+sorted_sc[i][j]+"'"]
            muts = ['HIS']

            # Get a compatible list of mutations. The above turns into
            # [('A', 22, 'DENQ')]
            mutations = protein.Mutator.convert_residue_list(residues, muts)

            # Construct the mutator
            mutator = protein.Mutator(reference_st, mutations)

            # Loop over each mutation
            for mutation in mutator.generate():
                #
                mutated_structure = mutation.struct
                residue_map = mutation.residue_map

                res_str = ", ".join(str(res) for res in residue_map.values())
                print 'Residues affected by this mutation: %s' % res_str

                # Do something with the mutated structure (refine maybe)

                writer.append(mutated_structure)

        os.system("/gpfs/projects/bsc72/Schrodinger_SMP2/utilities/structconvert -imae mutated_structures.mae -opdb merged_"+str(i)+"_"+ list_mutations[0]+"_"+list_mutations[1]+"_"+list_mutations[2]+"_"+list_mutations[3].rstrip()+".pdb")
        i = i + 1
        if i>5000:
            break


def propka():
    pH=7
    for filename in os.listdir("."):
        if filename.startswith("merged_"):
            os.system("/gpfs/projects/bsc72/Schrodinger_SMP2/utilities/prepwizard -propka_pH %s -f 2005 %s %s" %(pH,filename,"propka"+filename))

def minimization():
    for filename in os.listdir("."):
        if filename.startswith("propka"):
            st = structure.StructureReader(filename).next()
            writer = structure.StructureWriter("min"+filename)
            minimize.minimize_structure(st)
            st.append("min"+filename)



def nitro_filter(infile):
   # print("Extracting nitrogen coordinates")
   # print (datetime.datetime.now())

    time.sleep(240)
    uu=[]
    for line in infile:
        uu.append(line)
    sorted_score = [l.split(', ') for l in uu]

    sorted_sc = sorted(sorted_score, key=itemgetter(0))
    j=0
    pdbs = [file for file in os.listdir(".") if file.startswith("propkamerged_")]
    print(sorted_sc[0])
    for j, pdb in enumerate(pdbs):
        print pdb
        match=re.search(r'(\d+)',pdb)
        print(match.group())
        i=int(match.group())
        #i = int(pdb.split("_A")[0].split("_")[-1])
        #print i
    #for filename in os.listdir('.'):
     #   if filename.startswith("propkamerged_"):
        input = open(pdb, "rt")
        previous_line='hola'
        for line in input:
            if (('ND1 HIS A '+str(sorted_sc[i][1])+' ') in line or ('ND1 HIS A  '+str(sorted_sc[i][1])+' ') in line or ('ND1 HIS A   '+str(sorted_sc[i][1])+' ') in line  or ('ND1 HIS A '+str(sorted_sc[i][2])+' ') in line or ('ND1 HIS A  '+str(sorted_sc[i][2])+' ') in line or ('ND1 HIS A   '+str(sorted_sc[i][2])+' ') in line or ('ND1 HIS A '+str(sorted_sc[i][3])+' ') in line or ('ND1 HIS A  '+str(sorted_sc[i][3])+' ') in line or ('ND1 HIS A   '+str(sorted_sc[i][3])+' ') in line or ('ND1 HIS A '+str(sorted_sc[i][4].rstrip())+' ') in line or ('ND1 HIS A  '+str(sorted_sc[i][4].rstrip())+' ') in line or ('ND1 HIS A   '+str(sorted_sc[i][4].rstrip())+' ') in line) and 'ATOM' in line:
                nd1=line
            elif (('NE2 HIS A '+str(sorted_sc[i][1])+' ') in line or ('NE2 HIS A  '+str(sorted_sc[i][1])+' ') in line or ('NE2 HIS A   '+str(sorted_sc[i][1])+' ') in line  or ('NE2 HIS A '+str(sorted_sc[i][2])+' ') in line or ('NE2 HIS A  '+str(sorted_sc[i][2])+' ') in line or ('NE2 HIS A   '+str(sorted_sc[i][2])+' ') in line  or ('NE2 HIS A '+str(sorted_sc[i][3])+' ') in line or ('NE2 HIS A  '+str(sorted_sc[i][3])+' ') in line or ('NE2 HIS A   '+str(sorted_sc[i][3])+' ') in line  or ('NE2 HIS A '+str(sorted_sc[i][4].rstrip())+' ') in line or ('NE2 HIS A  '+str(sorted_sc[i][4].rstrip())+' ') in line or ('NE2 HIS A   '+str(sorted_sc[i][4].rstrip())+' ') in line )and 'ATOM' in line:
                ne2=line
            elif (('HD1 HIS A '+str(sorted_sc[i][1])+' ') in line or ('HD1 HIS A  '+str(sorted_sc[i][1])+' ') in line or ('HD1 HIS A   '+str(sorted_sc[i][1])+' ') in line  or ('HD1 HIS A '+str(sorted_sc[i][2])+' ') in line or ('HD1 HIS A  '+str(sorted_sc[i][2])+' ') in line or ('HD1 HIS A   '+str(sorted_sc[i][2])+' ') in line  or ('HD1 HIS A '+str(sorted_sc[i][3])+' ') in line or ('HD1 HIS A  '+str(sorted_sc[i][3])+' ') in line or ('HD1 HIS A   '+str(sorted_sc[i][3])+' ') in line  or ('HD1 HIS A '+str(sorted_sc[i][4].rstrip())+' ') in line or ('HD1 HIS A  '+str(sorted_sc[i][4].rstrip())+' ') in line or ('HD1 HIS A   '+str(sorted_sc[i][4].rstrip())+' ') in line )and 'ATOM' in line and ('HE2 HIS A') not in previous_line:
                output.write(ne2)
            elif (('HE2 HIS A '+str(sorted_sc[i][1])+' ') in line or ('HE2 HIS A  '+str(sorted_sc[i][1])+' ') in line or ('HE2 HIS A   '+str(sorted_sc[i][1])+' ') in line  or ('HE2 HIS A '+str(sorted_sc[i][2])+' ') in line or ('HE2 HIS A  '+str(sorted_sc[i][2])+' ') in line or ('HE2 HIS A   '+str(sorted_sc[i][2])+' ') in line  or ('HE2 HIS A '+str(sorted_sc[i][3])+' ') in line or ('HE2 HIS A  '+str(sorted_sc[i][3])+' ') in line or ('HE2 HIS A   '+str(sorted_sc[i][3])+' ') in line  or ('HE2 HIS A '+str(sorted_sc[i][4].rstrip())+' ') in line or ('HE2 HIS A  '+str(sorted_sc[i][4].rstrip())+' ') in line or ('HE2 HIS A   '+str(sorted_sc[i][4].rstrip())+' ') in line )and 'ATOM' in line and ('HD1 HIS A') not in previous_line:
                output.write(nd1)
            previous_line=line


                #if ('HD1 HIS A  47') in line and 'ATOM' in line:
                  #ND1 = True
        j=j+1
        input.close()

def nitro_coordinates(input):
    print("Extracting nitrogen coordinates")
    geometry=[]

    #Pick the coordinates from pdb
    for line in input:
        geometry.append(line[31:38])   #x_coordinate
        geometry.append(line[39:46])   #y coordinate
        geometry.append(line[47:54])   #z coordinate
        geometry.append(line[23:26])   #residue number

    #Transform to float
    for i in range(len(geometry)):
        #if i==0 or (i+1)%4 != 0:
        geometry[i]=float(geometry[i])


    #Split the list in list of list
    geometry2=[geometry[i*16:(i+1)*16] for i in range((len(geometry)+16-1)//16)]
    geometry3=list(geometry2)
    #print geometry2
    for j in range(len(geometry2)):
        geometry3[j]=[geometry2[j][i*4:(i+1)*4] for i in range(4)]#el range4 era un 3
    #print geometry3
    return geometry3

def nitro_centroid(nitro_list):
    print("Calculating centroid")
    zero_list = []
    three_zeros = [0.000] * 3
    for i in range(len(nitro_list)):
        zero_list.append(three_zeros)
        print i
        print nitro_list[i]
    centroid = []
    #print nitro_list[0]
    #print nitro_list[88]
    np.asarray(nitro_list)
    np.asarray(centroid)

    # Calculate centrioid
    #print(combinated[0])
    for i in range(len(nitro_list)):
        centroid.append(np.mean(nitro_list[i],axis=0))
    #print centroid[0]
    uu=[]
    for i in range(len(centroid)):
        uu.append(np.array(centroid[i]).tolist())
        #print (uu[i])
        del uu[i][3] #borrar la mitjana


    for i in range(len(centroid)):
        for j in range(0,4):#el 4 era un 3
            uu[i].append(nitro_list[i][j][3])

    return uu

def nitro_module(combination_list,centroid_list,n):
    print("calculating module")

    module_list = [ [] for i in range(len(combination_list)) ]
    #print(centroid_list)
    #Calculate the module from each vertix of the prism to the centroid
    for i in range(len(combination_list)):
        for j in range(0,2*n):
            if j < n:
                module_list[i].append(math.sqrt(pow(combination_list[i][j][0]-centroid_list[i][0],2)+pow(combination_list[i][j][1]-centroid_list[i][1],2)+pow(combination_list[i][j][2]-centroid_list[i][2],2)))
            else:
                module_list[i].append(centroid_list[i][j-1])  #j-1 amb 4
    print(module_list)
    return module_list


def nitro_score(module_list):
    print("Calculating score")

    score_function= [ [] for i in range(len(module_list)) ]
    #Calculate the score of each possible combination of alpha carbons

    for i in range(len(module_list)):
        score_function[i].append(math.sqrt(pow(module_list[i][0]-2.0,2) + pow(module_list[i][1] - 2.0,2) + pow(module_list[i][2]- 2.0,2)+pow(module_list[i][3]- 2.0,2)))  #Aquests numeros son les distancies entre els carbonis alfa i centroide del metal binding site que vull mimetitzar
        score_function[i].append(module_list[i][4])
        score_function[i].append(module_list[i][5])
        score_function[i].append(module_list[i][6])
        score_function[i].append(module_list[i][7])


    sorted_score=sorted(score_function, key=itemgetter(0))
    #output.write("%s" %len(sorted_score))


        #["A:{}".format(prova[i][j]) for j in range(1, 4)]

    for item in sorted_score:
        output.write("%s\n" %item)
    print ("End")
    return sorted_score




input = open ("PELE_PaDaI.pdb", "rt")
output = open("filtred", "wt")
filter(input)
input.close()
output.close()
input = open ("filtred", 'rt')
combination_list=combination(list_of_coordinates(input),4)
input.close()
centroid_list=centroid(combination_list,4)
module_list=module(combination_list,centroid_list,4)
output = open("final_score", "wt")
#formated_function, sorted_score=score(module_list)


#infile =open("final_score_changed","rt")
#outfile =open("logfile","wt")
#mutate(infile,outfile)
#outfile.close()
#infile.close()

#output.close()
propka()
#minimization()

infile =open("final_score_changed","rt")
output = open("nitrogen_filtred", "wt")
nitro_filter(infile)
infile.close()
output.close()

input = open("nitrogen_filtred", "rt")
nitro_list=nitro_coordinates(input)
input.close()
centroid_list=nitro_centroid(nitro_list)
nitro_module_list=nitro_module(nitro_list, centroid_list, 4)
output= open("final_nitro_score", "wt")
nitro_score(nitro_module_list)
output.close()

