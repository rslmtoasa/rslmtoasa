# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
PAPER         ?=
BUILDDIR      ?= _build

# Internal variables.
PAPEROPT_letter = -D latex_paper_size=letter
PAPEROPT_a4     = -D latex_paper_size=a4
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source
# the i18n builder cannot share the output directory and must create a temporary
# one, hence we create a tmp directory.
I18NDIR         = $(BUILDDIR)/locale

# the i18n builder cannot share the output directory and must create a temporary
# one, hence we create a temporary directory for the documents in this directory.
.PHONY: help Makefile

# Put it first so that "make" without argument is still useful for
help:
	@$(SPHINXBUILD) -M help $(ALLSPHINXOPTS) $(BUILDDIR)

.DEFAULT_TARGET := help

%: Makefile
	@$(SPHINXBUILD) -M $@ $(ALLSPHINXOPTS) $(BUILDDIR)
