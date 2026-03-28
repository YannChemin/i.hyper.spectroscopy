MODULE_TOPDIR = ../..

PGM = i.hyper.spectroscopy

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make

default: script html $(TEST_DST)

# Ensure HTML manual is installed
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html
