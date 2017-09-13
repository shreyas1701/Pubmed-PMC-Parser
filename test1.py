# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 12:01:23 2017

@author: shrey
"""
import csv
import test_pubmed as pp
from Bio import Entrez
###pubmed database
pubmed_dict_list=[]
Entrez.email = "shreyas1701@gmail.com"
data = Entrez.esearch(db="pubmed", retmax='100',sort="relevance",term = "breast+cancer")#can use sort=best match
res=Entrez.read(data)
PMID = res["IdList"]

for each in PMID:
    path_xml = pp.load_xml(each) # loads xml and returns tree
    pubmed_dict=pp.parse_pubmed_web_tree(path_xml) # dictionary output
    pubmed_dict_list.append(pubmed_dict)   

with open('C:\Research\dict_pub_bs.csv', 'w',encoding='UTF8',newline='') as f:  
    w = csv.DictWriter(f, pubmed_dict_list[0].keys())
    w.writerow(dict((fn,fn) for fn in pubmed_dict_list[0].keys()))
    w.writerows(pubmed_dict_list)


