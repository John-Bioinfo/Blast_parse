#!/usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = "Zongyun Qiao"
__copyright__ = "Copyright 2017, A Biotech"
__credits__ = [
    "Zongyun Qiao"]  # remember to add yourself
__license__ = "GPL"
__version__ = "0.1-dev, 20171204"
__maintainer__ = "Zongyun Qiao"
__email__ = "gulile@yeah.net"

import os, re
import json, io
import urllib.request
from urllib.request import Request
from urllib.request import urlopen
from urllib.parse import urlencode, quote_plus

ctURL = "https://clinicaltrials.gov/ct2/results"

headers = [('User-Agent', 'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-US; rv:1.9.1.6) Gecko/20091201 Firefox/3.5.6')]

### fullURL = "https://clinicaltrials.gov/ct2/results?term=EGFR&type=&rslt=With&recrs=e&cond=NSCLC"

#values = { 'term' : 'EGFR',
#           'type' : '',
#           'rslt' : 'With',
#           'recrs': 'e',
#           'cond' : 'NSCLC'  }

values1 = { 'term' : 'ALK',
           'type' : '',
           'rslt' : 'With',
           'recrs': 'e',
           'cond' : 'NSCLC'  }

## encode search values
d_paras = urlencode(values1, quote_via = quote_plus)

## ref : https://stackoverflow.com/questions/40557606/how-to-url-encode-in-python-3

data = d_paras.encode('ascii')

if not os.path.exists("database.json"):

    ## request url and get data, then store data into json file
    opener = urllib.request.build_opener()
    opener.addheaders = headers
    resources = opener.open(ctURL, data).read()
    decData = resources.decode("utf-8")

## ----------------------------------------------------------------------  
    jsonData = json.dumps({"test_data": decData} , ensure_ascii=False)
    
    # use dump , not dumps          
    with open("database.json", "w") as dump_file:
        json.dump(jsonData, dump_file)

# load data from the stored json file
with io.open("database.json", "r") as source_d:
    a = json.load(source_d)
    newDict = json.loads( a )
    ## json.loads(a) will return string object, and it need string input, as "a" above

newData = newDict["test_data"]

## res_data = re.findall("table \= \$\('#theDataTable'\)\.DataTable.*?\(\{\s+(.*?),\s+\"columns\"\:", decData, re.S)

## ----------------------------------------------------------------------

res_data = re.findall("table \= \$\('#theDataTable'\)\.DataTable.*?\(\{\s+(.*?),\s+\"columns\"\:", newData, re.S)

for res in res_data:
    print(res)
    if res.find("\"data\"") != -1:
        searchT = re.search( "var tableData1 = \[(.*)\]\s+];", newData, re.S)
        if searchT:
            tableData = searchT.group(1)
            # print(tableData)
            tables = re.findall("\"<a title=\\\"(.*?)<\\/a>", tableData, re.S)
            
            for t in tables:
                print(t)
