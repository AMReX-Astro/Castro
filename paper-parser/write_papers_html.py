# due to Unicode, this requires python 3+


import parser

papers = parser.parse_bibfile("papers.bib")


# sorted by topic
tf = open("papers.template")
dh = open("papers.html", "w")

subs = list({p.subject for p in papers})
print(subs)
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

    ostr += f"<header class='major'>\n<h3>{s}</h3>\n</header>\n"

    ostr += "<div class='table-wrapper'>\n"
    ostr += "  <table>\n"

    for p in ps:

        t, o, l = p.jstring()
        ostr += "<tr><td>"
        if not l == "":
            ostr += f"<a href='{l}'>{t}</a><br>\n"
        else:
            ostr += f"{t}<br>\n"

        ostr += f"{o}</td></tr>\n"

    ostr += "  </table>\n"
    ostr += "</div>\n"

for line in tf:
    dh.write(line.replace("@@pub-list@@", ostr))

dh.close()
tf.close()

