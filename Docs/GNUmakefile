EPStoPDF = epstopdf

ALL: CastroUserGuide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Preface/Preface.tex \
                  Introduction/Introduction.tex \
                  GettingStarted/CastroGettingStarted.tex \
	          Inputs/CastroInputs.tex \
                  Software/CastroSoftware.tex \
	          FlowChart/FlowChart.tex \
                  runtime_parameters/runtime_parameters.tex \
                  runtime_parameters/rp_intro.tex \
                  ManagingJobs/managingjobs.tex \
                  FAQ/faq.tex \
                  Equations/Equations.tex \
	          Gravity/CastroGravity.tex \
	          Diffusion/CastroDiffusion.tex \
	          Rotation/CastroRotation.tex \
                  Radiation/CastroRadiation.tex \
                  EOSNetwork/EOSNetwork.tex \
                  AMR/AMR.tex \
                  LevelSet/LevelSet.tex \
                  ConvertCheckpoint/ConvertCheckpoint.tex \
		  MAESTRO_restart/MAESTRO_restart.tex \
                  Visualization/CastroVisualization.tex \
                  Verification/Verification.tex \
                  castrosymbols.tex




CastroUserGuide.pdf: CastroUserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex CastroUserGuide.tex < /dev/null
	bibtex CastroUserGuide
	pdflatex CastroUserGuide.tex < /dev/null
	makeindex CastroUserGuide < /dev/null
	pdflatex CastroUserGuide.tex < /dev/null
	pdflatex CastroUserGuide.tex < /dev/null

pdf:	CastroUserGuide.pdf 

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source)
	$(RM) CastroUserGuide.blg CastroUserGuide.log CastroUserGuide.toc
	$(RM) CastroUserGuide.aux CastroUserGuide.ps CastroUserGuide.bbl CastroUserGuide.pdf
	$(RM) CastroUserGuide.lot CastroUserGuide.lof CastroUserGuide.out
	$(RM) CastroUserGuide.idx CastroUserGuide.ind CastroUserGuide.ilg
	$(RM) *~

.PHONY: clean
