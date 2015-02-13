EPStoPDF = epstopdf

ALL: CastroUserGuide.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard ConvertCheckpoint/*.eps) \
             $(wildcard Scaling/*.eps) \
             $(wildcard CastroVerification/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/CastroIntroduction.tex \
                  GettingStarted/CastroGettingStarted.tex \
	          Inputs/CastroInputs.tex \
                  UnitsAndConstants/CastroUnits.tex \
                  Equations/Equations.tex \
	          Gravity/CastroGravity.tex \
	          FlowChart/FlowChart.tex \
                  LevelSet/LevelSet.tex \
                  AMR/AMR.tex \
                  EOSNetwork/EOSNetwork.tex \
                  ConvertCheckpoint/ConvertCheckpoint.tex \
                  Visualization/CastroVisualization.tex \
                  Software/CastroSoftware.tex \
                  CastroVerification/Verification.tex \
                  Parallel/Parallel.tex \
                  Scaling/Scaling.tex \
                  Suggestions/Suggestions.tex \
		  MAESTRO_restart/MAESTRO_restart.tex \
	          Rotation/CastroRotation.tex

CastroUserGuide.dvi: CastroUserGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex CastroUserGuide.tex < /dev/null
	bibtex CastroUserGuide
	latex CastroUserGuide.tex < /dev/null
	latex CastroUserGuide.tex < /dev/null

CastroUserGuide.pdf: CastroUserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex CastroUserGuide.tex < /dev/null
	bibtex CastroUserGuide
	pdflatex CastroUserGuide.tex < /dev/null
	pdflatex CastroUserGuide.tex < /dev/null

pdf:	CastroUserGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) CastroUserGuide.dvi
	$(RM) CastroUserGuide.blg CastroUserGuide.log CastroUserGuide.toc
	$(RM) CastroUserGuide.aux CastroUserGuide.ps CastroUserGuide.bbl CastroUserGuide.pdf
	$(RM) *~

.PHONY: clean
