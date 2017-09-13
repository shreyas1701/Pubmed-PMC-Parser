# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:16:28 2017

@author: shrey
"""
from Bio import Entrez
import time
import requests
from lxml import html
from utils import stringify_children,stringify_affiliation_rec
from nltk import tokenize

def get_links_id(pmid):
	link_list = []
	links = Entrez.elink(dbfrom="pmc", id=pmid, linkname="pmc_pmc_citedby")	
	record = Entrez.read(links)
	try:
		records = record[0][u'LinkSetDb'][0][u'Link']
		for link in records:
			link_list.append(link[u'Id'])
		citations= len(link_list)
	except:
		citations = 0
	return citations

def load_pmc_xml(pmcid, sleep=None):
    """
    Load XML file from given pmid from eutils site
    return a dictionary for given pmid and xml string from the site
    sleep: how much time we want to wait until requesting new xml
    """
    link= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&retmode=xml&sort=relevance&id=%s" % str(pmcid)
    page = requests.get(link)    
    tree = html.fromstring(page.content)
    if sleep is not None:
        time.sleep(sleep)
    return tree

def parse_pmc_web_tree(tree):
     """
    Giving tree, return simple parsed information from the tree
    """
     
     tree_title = tree.xpath('//title-group/article-title')
     if tree_title is not None:
        title = ' '.join([stringify_children(a).strip() for a in tree_title])
     else:
        title = ''

     try:
        abstracts = list()
        abstract_tree = tree.xpath('//abstract')
        for a in abstract_tree:
            for t in a.itertext():
                text = t.replace('\n', ' ').replace('\t', ' ').strip()
                abstracts.append(text)
        abstract = ' '.join(abstracts)
     except:
        abstract = ''
  

     if tree.xpath('//journal-title-group/journal-title') is not None:
        journal = ';'.join([t.text.strip() for t in tree.xpath('//journal-title-group/journal-title')])
     else:
        journal = ''
     
     pubdate = tree.xpath('//pub-date[@pub-type="ppub"]')
     if len(pubdate) >= 1 and pubdate[0].find('year') is not None:
        year = pubdate[0].find('year').text
     else:
        year = ''
          
     # create affiliation dictionary
     affil_id = tree.xpath('//aff[@id]/@id')
     if len(affil_id) > 0:
        affil_id = list(map(str, affil_id))
     else:
        affil_id = ['']  # replace id with empty list

     affil_name = tree.xpath('//aff[@id]')
     affil_name_list = list()
     for e in affil_name:
        name = stringify_affiliation_rec(e)
        name = name.strip().replace('\n', ' ')
        affil_name_list.append(name)
     affiliation_list = [[idx, name] for idx, name in zip(affil_id, affil_name_list)]
    
     if tree.xpath('//article-meta/article-id[@pub-id-type="pmc"]') is not None:
        pmcid = ';'.join([p.text.strip() for p in tree.xpath('//article-meta/article-id[@pub-id-type="pmc"]')])
     else:
        pmcid = ''
    # print(pmcid)
     
     try:
         keywords = list()
         if tree.xpath('//kwd-group/kwd')is not None:
             for keyword in tree.xpath('//kwd-group/kwd'):
                    keywords.append(keyword.text)
             keywords_text ='; '.join(keywords)
         else:
             keywords_text= ' '   
     except:
         kwds =[]
         for c_node in tree.findall('.//kwd-group'):
             test_node = c_node.find('title')
             if test_node is not None and  'Keywords' in test_node.text:
                 for id_node in c_node.findall('kwd'):
                     for t in id_node.itertext():
                         kwds.append(t)   
         keywords_text= ';'.join(kwds)
         
     authors_tree = tree.xpath('//contrib-group/contrib[@contrib-type="author"]/name')
     authors = list()
     if authors_tree is not None:
        for a in authors_tree:
            firstname = a.find('given-names').text if a.find('given-names') is not None else ''
            lastname = a.find('surname').text if a.find('surname') is not None else ''
            fullname = (firstname + ' ' + lastname).strip()
            if fullname == '':
                fullname = a.find('collectivename').text if a.find('collectivename') is not None else ''
            authors.append(fullname)
        authors_text = '; '.join(authors)
     else:
        authors_text = ''
     #print(authors_text)
     
     if tree.xpath('//permissions/copyright-statement') is not None:
        copyright_text = ';'.join([c.text.strip() for c in tree.xpath('//permissions/copyright-statement')])
     else:
        copyright_text = '' 

     
     if tree.xpath('//article-meta/article-id[@pub-id-type="pmid"]') is not None:
        pmid = ';'.join([p.text.strip() for p in tree.xpath('//article-meta/article-id[@pub-id-type="pmid"]')])
     else:
        pmid = ''

     n_citations=get_links_id(pmcid)


     if tree.xpath('//article-meta/article-id[@pub-id-type="doi"]') is not None:
        doi = ';'.join([p.text.strip() for p in tree.xpath('//article-meta/article-id[@pub-id-type="doi"]')])
     else:
        doi = ''
     
     try:
         cons =[]
         for c_node in tree.findall('.//sec'):
             test_node = c_node.find('title')
             if test_node is not None and  'Conclusion' in test_node.text:
                id_node = c_node.find('p')
                for t in id_node.itertext():
                    text =(t.replace('\n', ' ').replace('\t', ' ').strip())
                    cons.append(text)
         conclusion=''.join(cons)
         if conclusion=="":
             con=tokenize.sent_tokenize(abstract)
             conclusion =''.join(con[len(con)-1])
     except:
         con=tokenize.sent_tokenize(abstract)
         conclusion =''.join(con[len(con)-1].strip())
#     print(conclusion)
         
# title,abstract,journal,affiliation,authors,year,keywords,conclusion,copyright,pmid,doi,citationindex
     dict_out = {'title': title,
                'abstract': abstract,
                'journal': journal,
                'affiliation': affiliation_list,
                'authors': authors_text,
                'year': year,
                'keywords': keywords_text,
                'conclusion':conclusion,
                'copyright':copyright_text,
                'pmid':pmid,
                'doi':doi,
                'citation_index':n_citations}
     return dict_out
        
