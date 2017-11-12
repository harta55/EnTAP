from __future__ import print_function

import os
import argparse
import sys
import re
from xml.etree import cElementTree as ET

# ******************* Defines ***********************
OUTPUT_FILE          = "ncbi_tax.entp"
TAX_ROOT             = "root[Subtree]"
TAX_DATABASE         = "taxonomy"
TAX_BASE_URL         = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
PACK_URLLIB2         = "urllib2"
PACK_URLLIB          = "urllib.request"
RET_MAX              = 10000
PERCENT              = 100.0
URL_TIMEOUT          = 5
VERSION_NUM          = '1'
VERSION_NUM2         = '0'
__version_info__     = (VERSION_NUM, VERSION_NUM2)
__version__          = '.'.join(__version_info__)

EXIT_SUCCESS = 0
EXIT_ERROR   = 1
EXIT_EXISTS  = 2
EXIT_VERSION = 3     # Python version not compatbile
# ***************************************************

# ******************* Globals ***********************
_outpath = ""
_version = 0.0
url_lib = None
# ***************************************************


def user_parse():
    global _outpath

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=('''\
        This script is used to download the taxonomic database for use by EnTAP from the NCBI
        Taxonomy database. It then has to be configured by EnTAP in order to be used.

        '''))

    parser.add_argument('-o', type=str, action='store',
                        default=os.path.join(os.getcwd(), OUTPUT_FILE),
                        help='Path to save Taxonomic database to')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    inputs = parser.parse_args()

    _outpath = inputs.o
    if os.path.exists(_outpath):
        print('File already exists at: ' + _outpath)
        exit(EXIT_EXISTS)


def verify_package(package):
    is_2 = 2 <= _version < 3
    has_package = False
    if is_2:
        import imp
        try:
            imp.find_module(package)
            has_package= True
        except ImportError:
            has_package = False
    elif _version >= 3.4:
        import importlib.util
        check = importlib.util.find_spec(package)
        has_package = check is not None
    elif not is_2 and _version < 3.4:
        import importlib
        spam_loader = importlib.find_loader(package)
        has_package = spam_loader is not None
    else:
        print("Python version not compatible.")
        exit(EXIT_VERSION)
    return has_package


def init_version():
    global _version
    _version = sys.version_info[0] + 0.1*sys.version_info[1]


def download_tax():
    global url_lib
    if verify_package(PACK_URLLIB):
        import urllib.request as url_lib
    elif verify_package(PACK_URLLIB2):
        import urllib2 as url_lib
    else:
        print("Unable to find a correct urllib package")
        exit(1)

    print('Downloading taxonomic database...')
    tax_results = get_ncbi_query_results(TAX_DATABASE, TAX_ROOT)
    fetch_url(tax_results, TAX_DATABASE)


def get_ncbi_query_results(database, query):
    web_env = ""
    key     = 0
    count   = 0

    query_url = TAX_BASE_URL + "esearch.fcgi?db={}&term={}&usehistory=y".format(database, query)

    try:
        response = url_lib.urlopen(query_url, timeout=URL_TIMEOUT)
        result = response.read().decode('utf-8')

        exp = re.compile('<WebEnv>(\S+)<\/WebEnv>')
        web_env = exp.search(result).group()[8:-9]

        exp = re.compile('<QueryKey>(\S+)<\/QueryKey>')
        key = int(exp.search(result).group()[10:-11])

        exp = re.compile('<Count>(\S+)<\/Count>')
        count = int(exp.search(result).group()[7:-8])
    except ...:
        print("Error in searching tax database: " + query_url)
        exit(1)
    return [web_env, key, count]


def fetch_url(search_results, database):
    response = None

    out = open(_outpath, 'w')
    for x in range(0, int(search_results[2]/RET_MAX + 2)):
        try:
            url = TAX_BASE_URL + "efetch.fcgi?db={}&query_key={}&WebEnv={}"    \
                "&retstart={}&retmax={}&rettype=null&retmode=xml". \
                format(database, str(search_results[1]), search_results[0], str(x*RET_MAX-RET_MAX), RET_MAX)
            response = url_lib.urlopen(url)
            data = response.read().decode('utf-8')
            out.writelines(process_fetch(data))
            print("Status: {:0.2f}%".format(PERCENT*float((x*RET_MAX)/(search_results[2]+RET_MAX))))
        except Exception as e:
            print("Error in fetching data")
            print(e)
            exit(1)
    print('Download complete!')
    out.close()


def process_fetch(data):
    lines = []
    lineage = ""
    sci_name = ""
    othernames = None
    synonym = ""

    processed = ET.fromstring(data)

    for tax in processed.findall('Taxon'):
        lineage = tax.find('Lineage').text
        if lineage is not None: lineage = lineage.lower()
        sci_name = tax.find('ScientificName').text.lower()
        tax_id = tax.find('TaxId').text
        othernames = tax.find('OtherNames')
        if othernames is not None:
            syn = othernames.findall('Synonym')
            for el in syn:
                synonym = el.text.replace("'", "").lower()
                if not synonym: continue
                lines.append("{}\t{}\t{}\n".format(synonym, tax_id, lineage))
        lines.append("{}\t{}\t{}\n".format(sci_name, tax_id, lineage))
    return lines


# **************************** End Declarations ********************************#


# ************* Entry ****************
def main():
    user_parse()
    init_version()
    download_tax()


if __name__ == "__main__":
    main()