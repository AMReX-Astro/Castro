EPStoPDF = epstopdf

ALL: CastroUserGuide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/CastroIntroduction.tex \
                  GettingStarted/CastroGettingStarted.tex \
	          Inputs/CastroInputs.tex \
                  UnitsAndConstants/CastroUnits.tex \
                  Software/CastroSoftware.tex \
                  Parallel/Parallel.tex \
	          FlowChart/FlowChart.tex \
                  runtime_parameters/runtime_parameters.tex \
                  runtime_parameters/rp_intro.tex \
                  ManagingJobs/managingjobs.tex \
                  FAQ/faq.tex \
                  Equations/Equations.tex \
	          Gravity/CastroGravity.tex \
	          Diffusion/CastroDiffusion.tex \
	          Rotation/CastroRotation.tex \
                  EOSNetwork/EOSNetwork.tex \
                  AMR/AMR.tex \
                  LevelSet/LevelSet.tex \
                  ConvertCheckpoint/ConvertCheckpoint.tex \
		  MAESTRO_restart/MAESTRO_restart.tex \
                  Visualization/CastroVisualization.tex \
                  CastroVerification/Verification.tex \
                  castrosymbols.tex




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
	$(RM) CastroUserGuide.lot CastroUserGuide.lof CastroUsetGuide.out
	$(RM) *~

.PHONY: clean
