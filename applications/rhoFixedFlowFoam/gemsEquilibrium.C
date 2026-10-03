/*------------------------------------------------------------------------------

gemsEquilibrium.C

PURPOSE

Implementation of the optional GEMS3K backend.  Everything that touches the
bridge is inside #ifdef LESTO_HAVE_GEMS, which the generated gemsConfig.H
defines only for a build with LESTO_GEMSBRIDGE set.  Without it the class is
a stub whose requireBridge() stops the run with a rebuild hint.

gemsConfig.H must be included here, before any test of LESTO_HAVE_GEMS: the
include is what makes wmake recompile this file when the bridge is toggled.

The bridge header is included with angle brackets on purpose.  wmkdepend
follows only quoted includes and ignores #ifdef, so a quoted include would
make it report a missing gemsbridge.h in every build without the bridge.
wmkdepend therefore does not record gemsbridge.h itself.  Instead, a bridge
build's gemsConfig.H carries the checksum of gemsbridge.h, so an edit of the
bridge header rewrites gemsConfig.H and recompiles this file.

------------------------------------------------------------------------------*/

#include "gemsConfig.H"
#include "gemsEquilibrium.H"
#include "error.H"

#ifdef LESTO_HAVE_GEMS
#include <gemsbridge.h>
#endif

bool LESTO::gemsEquilibrium::available() {
#ifdef LESTO_HAVE_GEMS
  return true;
#else
  return false;
#endif
}


Foam::string LESTO::gemsEquilibrium::bridgeDirectory() {
  return Foam::string(LESTO_GEMSBRIDGE_DIR);
}


Foam::string LESTO::gemsEquilibrium::buildState() {
  if (available()) {
    return "GEMS3K bridge: compiled in from " + bridgeDirectory();
  }
  return "GEMS3K bridge: not compiled in (set LESTO_GEMSBRIDGE and rerun "
         "Allwmake to enable it)";
}


void LESTO::gemsEquilibrium::requireBridge(const Foam::dictionary& dict) {

  using Foam::nl;

  if (!available()) {
    FatalIOErrorInFunction(dict)
      << "equilibrium GEMS needs rhoFixedFlowFoam built with the GEMS3K "
      << "bridge, but this executable was built without it." << nl
      << "Rebuild with" << nl
      << "    export LESTO_GEMSBRIDGE=<directory with gemsbridge.h and "
      << "libgemsbridge.so>" << nl
      << "    ./Allwmake" << nl
      << "or select equilibrium table." << Foam::exit(Foam::FatalIOError);
  }

  FatalIOErrorInFunction(dict)
    << "The GEMS3K backend is not implemented yet; this executable only "
    << "links the bridge from " << bridgeDirectory() << "." << nl
    << "Select equilibrium table." << Foam::exit(Foam::FatalIOError);
}
