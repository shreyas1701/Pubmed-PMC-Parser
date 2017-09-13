# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 02:57:48 2017

@author: shrey
"""

import csv
import test_pmc as mp
from Bio import Entrez

####pmc database
pmc_dict_list=[]
Entrez.email = "shreyas1701@gmail.com"
data_pmc=Entrez.esearch(db="pmc", retmax='100',sort='relevance',term = "breast+cancer")
res_pmc=Entrez.read(data_pmc)
PMCID = res_pmc["IdList"]

for each in PMCID:
    path_xml = mp.load_pmc_xml(each) # loads xml and returns tree
    pmc_dict=  mp.parse_pmc_web_tree(path_xml) # dictionary output
    pmc_dict_list.append(pmc_dict)   

with open('C:\Research\dict_pmc.csv', 'w',encoding='UTF8',newline='') as f:  
    w = csv.DictWriter(f, pmc_dict_list[0].keys())
    w.writerow(dict((fn,fn) for fn in pmc_dict_list[0].keys()))
    w.writerows(pmc_dict_list)