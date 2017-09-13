# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:24:45 2017

@author: shrey
"""

import time
from Bio import Entrez
import requests
from lxml import etree
from lxml import html
from utils import stringify_children
from nltk import tokenize



def get_links_id(pmid):
	 link_list = []
	 links = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")	
	 record = Entrez.read(links)
	 try:
	 	 records = record[0][u'LinkSetDb'][0][u'Link']
	 	 for link in records:
		     link_list.append(link[u'Id'])
	 	 citations = len(link_list)
	 except:
	 	 citations = 0
	 return citations

def load_xml(pmid, sleep=None):
    """
    Load XML file from given pmid from eutils site
    return a dictionary for given pmid and xml string from the site
    sleep: how much time we want to wait until requesting new xml
    """
    link = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=%s" % str(pmid) 
    page = requests.get(link)  
    tree = html.fromstring(page.content)
    if sleep is not None:
        time.sleep(sleep)
    return tree


def parse_pubmed_web_tree(tree):
    """
    Giving tree, return simple parsed information from the tree
    """

    if tree.xpath('//articletitle') is not None:
        title = ' '.join([title.text for title in tree.xpath('//articletitle')])
    else:
        title = ''

    abstract_tree = tree.xpath('//abstract/abstracttext')
    abstract = ' '.join([stringify_children(a).strip() for a in abstract_tree])


    con=tokenize.sent_tokenize(abstract)
    conclusion =''.join(con[len(con)-1])

    if tree.xpath('//article//title') is not None:
        journal = ';'.join([t.text.strip() for t in tree.xpath('//article//title')])
    else:
        journal = ''

    pubdate = tree.xpath('//pubmeddata//history//pubmedpubdate[@pubstatus="medline"]')
    if len(pubdate) >= 1 and pubdate[0].find('year') is not None:
        year = pubdate[0].find('year').text
    else:
        year = ''

    affiliations = list()
    if tree.xpath('//affiliationinfo/affiliation') is not None:
        for affil in tree.xpath('//affiliationinfo/affiliation'):
            affiliations.append(affil.text)
        affiliations_text = '; '.join(affiliations)
    else:
         affiliations_text= '' 

    keywords = list()
    if tree.xpath('//keywordlist/keyword')is not None:
         for keyword in tree.xpath('//keywordlist/keyword'):
            keywords.append(keyword.text)
         keywords_text = '; '.join(keywords)
    else:
         keywords_text= ''    
         
         
    authors_tree = tree.xpath('//authorlist/author')
    authors = list()
    if authors_tree is not None:
        for a in authors_tree:
            firstname = a.find('forename').text if a.find('forename') is not None else ''
            lastname = a.find('lastname').text if a.find('forename') is not None else ''
            fullname = (firstname + ' ' + lastname).strip()
            if fullname == '':
                fullname = a.find('collectivename').text if a.find('collectivename') is not None else ''
            authors.append(fullname)
        authors_text = '; '.join(authors)
    else:
        authors_text = ''
        
    if tree.xpath('//abstract/copyrightinformation') is not None:
        copyright_text = ';'.join([c.text.strip() for c in tree.xpath('//abstract/copyrightinformation')])
    else:
        copyright_text = '' 
        
    if tree.xpath('//articleidlist/articleid[@idtype="pubmed"]') is not None:
        pmid = ';'.join([p.text.strip() for p in tree.xpath('//articleidlist/articleid[@idtype="pubmed"]')])
    else:
        pmid = ''

    if tree.xpath('//articleidlist/articleid[@idtype="doi"]') is not None:
        doi = ';'.join([p.text.strip() for p in tree.xpath('//articleidlist/articleid[@idtype="doi"]')])
    else:
        doi = ''
    
    n_citations=get_links_id(pmid)

# title,abstract,journal,affiliation,authors,year,keywords,conclusion,copyright,pmid,doi,citationindex
    dict_out = {'title': title,
                'abstract': abstract,
                'journal': journal,
                'affiliation': affiliations_text,
                'authors': authors_text,
                'year': year,
                'keywords': keywords_text,
                'conclusion':conclusion,
                'copyright':copyright_text,
                'pmid':pmid,
                'doi':doi,
                'citation_index':n_citations}
    return dict_out

def parse_xml_web(pmid, sleep=None, save_xml=False):
    """
    Give pmid, load and parse xml from Pubmed eutils
    if save_xml is True, save xml output in dictionary
    """
    tree = load_xml(pmid, sleep=sleep)
    dict_out = parse_pubmed_web_tree(tree)
    dict_out['pmid'] = str(pmid)
    if save_xml:
        dict_out['xml'] = etree.tostring(tree)
    return dict_out



