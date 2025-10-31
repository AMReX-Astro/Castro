import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import *
import re
import urllib.request

"""
This is a general parser that will take ADS bibtex listings for
papers and output a list of Paper objects that contain the
bibliographic information in convenient form.  This can then be used
to write a webpage listing of papers.

It optionally supports a list of ADS URLs, and will fetch the bibtex
for each paper.
"""


replace_str = {
    r"$^{4}$": r"<sup>4</sup>"
}

mdash1 = "-{2,3}"

class Paper:
    def __init__(self, authors, title, year, journal,
                 month=None, booktitle=None, editors=None,
                 volume=None, pages=None, link=None, note=None,
                 subject=None):

        self.authors = list(authors)
        self.title = title
        self.year = int(year)
        self.journal = journal
        self.month = month
        self.booktitle = booktitle
        self.editors = editors
        self.volume = volume
        self.pages = pages
        self.link = link
        self.note = note
        self.subject = subject

    def __lt__(self, other):
        if not self.year == other.year:
            return self.year < other.year
        else:
            if not (self.month == None or other.month == None):
                return self.month < other.month
            else:
                return self.year < other.year

    def jstring(self):
        t_str = self.title
        for k, v in replace_str.items():
            t_str = t_str.replace(k,v)
        t_str = re.sub(mdash1, "&mdash;", t_str)

        out_str = name_string(self.authors) + " "
        out_str += f"{self.year}, "

        if not self.journal == None:
            out_str += f"{self.journal}, "

        if not self.booktitle == None:
            out_str += f"in {self.booktitle}, "

        if not self.editors == None:
            out_str += f"ed. {name_string(self.editors)}, "

        if not self.volume == None:
            out_str += f"{self.volume}, "

        if not self.pages == None:
            out_str += f"p. {self.pages}, "

        if not self.note == None:
            out_str += f"{self.note}, "

        out_str = out_str.strip()

        if len(out_str) > 0:
            if out_str[len(out_str)-1] == ",":
                out_str = out_str[:len(out_str)-1]

        if not self.link == None:
            l_str = f"{self.link}"
        else:
            l_str = ""

        return t_str, out_str, l_str


def name_string(names):
    nm_str = ""
    if len(names) == 1:
        nm_str = f"{names[0]}"
    else:
        for n, a in enumerate(names):
            if n < len(names)-1:
                astr = "{}, "
            else:
                astr = "&amp; {}"
            nm_str += astr.format(a)
    return nm_str

def get_item(dict, name):
    if name in dict.keys():
        return dict[name]
    else:
        return None

def translate_journal(j):
    if j == None:
        return None
    else:
        jn = j["name"].strip()
        if jn.lower() == r"\apj":
            return "ApJ"
        elif jn.lower() == r"\apjs":
            return "ApJS"
        elif jn.lower() == r"\mnras":
            return "MNRAS"
        else:
            return jn

def fix_pages(p):
    if p == None:
        return None
    else:
        return p.replace("--","&ndash;")

def clean_names(a):
    if a == None:
        return None
    else:
        a_new = []
        for name in a:
            a_new.append(name.replace("{","").replace("}","").replace("~"," "))
        return a_new

def clean_ednames(a):
    if a == None:
        return None
    else:
        a_new = []
        for ed_dict in a:
            a_new.append(ed_dict["name"].replace("{","").replace("}","").replace("~"," "))
        return a_new

def customizations(record):
    """Use some functions delivered by the library

    :param record: a record
    :returns: -- customized record

    """
    record = convert_to_unicode(record)
    record = type(record)    # lowercase
    record = author(record)
    record = editor(record)
    record = journal(record)
    record = keyword(record)
    record = link(record)
    record = page_double_hyphen(record)
    record = doi(record)
    return record


def extract_paper_info(e):
    """ take a BibDatabase entry and make a Paper object from it """

    if not "title" in e.keys():
        print( "no title: ", e)
        return None
    else:
        title = e["title"]

    if not "author" in e.keys():
        print( "no author: ", e)
        return None
    else:
        authors = e["author"]

    authors = clean_names(authors)

    volume = get_item(e, "volume")
    journal = translate_journal(get_item(e, "journal"))
    year = get_item(e, "year")
    month = get_item(e, "month")
    editors = clean_ednames(get_item(e, "editor"))
    booktitle = get_item(e, "booktitle")
    pages = fix_pages(get_item(e, "pages"))
    note = get_item(e, "note")
    subject = get_item(e, "subject")

    if "adsurl" in e.keys():
        link = get_item(e, "adsurl")
    else:
        l = get_item(e, "link")
        if not l == None:
            link = l[0]["url"]
        else:
            link = None

    return Paper(authors, title, year, journal,
                 month=month, editors=editors,
                 booktitle=booktitle,
                 volume=volume, pages=pages,
                 link=link, note=note, subject=subject)


def parse_urlfile(url_file):
    """
    take a file of the form

    category: ads url

    and get the bibtex from the URL and return a list of Paper objects
    with the category stored as the subject

    """

    papers = []

    with open(url_file) as f:

        parser = BibTexParser(common_strings=True)
        parser.customization = customizations

        for line in f:
            if line.startswith("#") or line.strip() == "": continue

            subject, url = line.split(": ")

            # for the ADS bibtex URL, lop off the paper_id
            paper_id = url.strip().split("/")[-1]
            bibtex_url = f"http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode={paper_id}&data_type=BIBTEX"

            # get the bibtex in html -- this is a little tricky, since
            # urlopen gives us a byte object that we need to decode
            # into unicode before we can play with it.
            print(bibtex_url)
            with urllib.request.urlopen(bibtex_url) as response:
                bibtex_html = response.read()

            raw_bibtex_html = bibtex_html.splitlines()

            bibtex_string = ""
            for line in raw_bibtex_html:
                bibtex_string += "{}\n".format(line.decode("utf8"))

            print(bibtex_string)

            continue
            # strip off any header and just leave the bibtex
            found_start = False
            bibtex = ""
            for line in bibtex_string:
                if line.startswith("@"):
                    found_start = True
                if found_start:
                    bibtex += line

            print(bibtex)
            # parse the bibtex string
            # we need to ensure that the month string is quoted
            new_str = ""
            for line in bibtex.splitlines():
                if line.strip().startswith("month"):
                    fields = line.strip().split("=")
                    new_str += "month = '{}',\n".format(fields[-1].strip().replace(",",""))
                else:
                    new_str += line + "\n"

            print(new_str)

            bib_database = bibtexparser.loads(new_str, parser=parser)

            for e in bib_database.entries:
                p = extract_paper_info(e)
                if not e is None:
                    p.subject = subject
                    papers.append(p)

    papers.sort(reverse=True)
    return papers


def parse_bibfile(bibfile):

    with open(bibfile) as bibtex_file:
        parser = BibTexParser(common_strings=True)
        parser.customization = customizations
        bib_database = bibtexparser.load(bibtex_file, parser=parser)

        papers = []

        for e in bib_database.entries:
            p = extract_paper_info(e)
            if not e is None:
                papers.append(p)

    papers.sort(reverse=True)

    return papers
