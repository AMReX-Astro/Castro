# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = -E -a
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = castro
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

clean:
	rm -rf build doxy_files source/class source/file source/namespace source/preprocessed_files
	rm -f source/classlist.rst source/filelist.rst source/namespacelist.rst
	rm -f source/yt_example.rst
	rm -rf source/preprocessed_files
	rm -rf source/*_files.rst

html:
	./rp.py > source/runtime_parameters.rst
ifneq ($(NO_DOXYGEN),TRUE)
	python3 preprocess_files.py
	doxygen Doxyfile
	breathe-apidoc --o source doxy_files/xml -g class,file
	python3 make_api.py
endif
	python parse_changelog.py ../CHANGES.md > source/changelog.md
	jupyter nbconvert --to rst source/yt_example.ipynb
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	python3 preprocess_files.py
ifneq ($(NO_DOXYGEN),TRUE)
	doxygen Doxyfile
	breathe-apidoc --o source doxy_files/xml -g class,file
endif
	python3 make_api.py
	jupyter nbconvert --to rst source/yt_example.ipynb
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
