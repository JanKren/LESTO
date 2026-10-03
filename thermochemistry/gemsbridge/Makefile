# libgemsbridge.so: plain-C interface to GEMS3K for OpenFOAM (built with an
# older gcc).  GEMS3K, libstdc++ and libgcc are linked in statically and
# hidden; only the gemsb_* symbols are exported (see gemsbridge.map).
# Prerequisite: ./build-gems3k-static.sh
GEMS3K_SRC ?= $(HOME)/src/GEMS3K
GEMS_ENV   ?= $(HOME)/opt/gems-env
CXX        := $(GEMS_ENV)/bin/x86_64-conda-linux-gnu-g++
GEMS3K_LIB := $(GEMS3K_SRC)/build-static-bridge/lib/libGEMS3K-static.a
CXXFLAGS   := -std=c++20 -O2 -fPIC -pthread -DNDEBUG -DNO_THERMOFUN -DFMT_HEADER_ONLY \
              -fvisibility=hidden -fvisibility-inlines-hidden \
              -I$(GEMS3K_SRC) -I$(GEMS3K_SRC)/GEMS3K -I$(GEMS_ENV)/include

libgemsbridge.so: gemsbridge.cpp gemsbridge.h gemsbridge.map $(GEMS3K_LIB)
	$(CXX) $(CXXFLAGS) -shared -o $@ gemsbridge.cpp $(GEMS3K_LIB) \
	  -static-libstdc++ -static-libgcc -Wl,--exclude-libs,ALL \
	  -Wl,--version-script=gemsbridge.map -Wl,--no-undefined

check: libgemsbridge.so
	@echo "exported: $$(nm -D --defined-only $< | grep -c ' gemsb_') gemsb_*, $$(nm -D --defined-only $< | grep -vc ' gemsb_') other"
	@readelf -d $< | grep NEEDED

clean:
	rm -f libgemsbridge.so
