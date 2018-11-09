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
	rm -rf build doxy_files source/class source/file source/namespace
	rm -f source/classlist.rst source/filelist.rst source/namespacelist.rst
	rm -rf ../docs

html:
	./rp.py > source/runtime_parameters.rst
	doxygen Doxyfile
	breathe-apidoc --o source doxy_files/xml
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	cp -r build/html/* ../docs/

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	doxygen Doxyfile
	breathe-apidoc --o source doxy_files/xml
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
