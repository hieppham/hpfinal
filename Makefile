WGET?=wget
TAR?=tar
ZIP?=zip

PVRPTW_GALIB_BUNDLE=hgafinal
PVRPTW_GALIB_BUNDLE_ZIP=$(PVRPTW_GALIB_BUNDLE).zip
PVRPTW_GALIB_BUNDLE_CONTENT=$(GALIB_DIR) Makefile src

GALIB_DIR=galib247
GALIB_TGZ=$(GALIB_DIR).tgz
GALIB_URL="http://lancet.mit.edu/ga/dist/$(GALIB_TGZ)"
GALIB_LIB=$(GALIB_DIR)/ga/libga.a

.PHONY: all build_dir clean distclean galib $(PVRPTW_GALIB_BUNDLE_ZIP)

all: build_dir
	$(MAKE) all -C build

build_dir: $(GALIB_LIB)
	test -d build || mkdir build
	cd build && cmake -D GALIB_DIR=../$(GALIB_DIR) ../src

$(GALIB_LIB): galib
galib: $(GALIB_DIR)
	$(MAKE) lib -C $(GALIB_DIR)

$(GALIB_DIR):
	test -f $(GALIB_TGZ) || $(WGET) $(GALIB_URL)
	$(TAR) xzvf $(GALIB_TGZ)

clean:
	rm -rfv build

distclean: clean
	rm -rfv $(GALIB_DIR) $(GALIB_TGZ)
	rm -rfv $(PVRPTW_GALIB_BUNDLE) $(PVRPTW_GALIB_BUNDLE_ZIP)

$(PVRPTW_GALIB_BUNDLE_ZIP): distclean
	$(MAKE) $(GALIB_DIR)
	mkdir $(PVRPTW_GALIB_BUNDLE)
	cp -Rv $(PVRPTW_GALIB_BUNDLE_CONTENT) $(PVRPTW_GALIB_BUNDLE)
	$(ZIP) -r $@ $(PVRPTW_GALIB_BUNDLE)
