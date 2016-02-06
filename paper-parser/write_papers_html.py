# due to Unicode, this requires python 3+

from __future__ import print_function

import parser

papers = parser.parse_urlfile("castro-papers.txt")


# sorted by topic
tf = open("papers.template", "r")
dh = open("papers.html", "w")

subs = list(set([p.subject for p in papers]))
subs.sort(key=str.lower)

papers_by_subj = {}

for p in papers:
    subj = p.subject
    if not subj in papers_by_subj.keys():
        papers_by_subj[subj] = [p]
    else:
        papers_by_subj[subj].append(p)


# now loop over subject
ostr = ""
for s in sorted(papers_by_subj, key=str.lower):
    ps = papers_by_subj[s]
    ps.sort(reverse=True)

    ostr += "<header class='major'>\n<h3>{}</h3>\n</header>\n".format(s)

    ostr += "<div class='table-wrapper'>\n"
    ostr += "  <table>\n"

    for p in ps:

        t, o, l = p.jstring()
        ostr += "<td><td>"
        if not l == "":
            ostr += "<a href='{}'>{}</a><br>\n".format(l, t)
        else:
            ostr += "{}<br>\n".format(t)

        ostr += "{}</td></tr>".format(o)

    ostr += "  </table>\n"
    ostr += "</div>\n"

for line in tf:
    dh.write(line.replace("@@pub-list@@", ostr))

dh.close()
tf.close()

