/*------------------------------------------------------------------------------

gemsEquilibrium.H

PURPOSE

Optional GEMS3K backend of the phase-change model: equilibrium vapour
pressures from a Gibbs-energy minimisation instead of a p_v table.

The GEMS3K bridge (thermochemistry/gemsbridge, a plain-C ABI exported by
libgemsbridge.so) is an OPTIONAL build dependency:

  - it is compiled in only when LESTO_GEMSBRIDGE names a directory holding
    gemsbridge.h and libgemsbridge.so at build time (see Make/options);
  - Make/gemsConfig.sh records that decision in the generated header
    gemsConfig.H, which defines LESTO_HAVE_GEMS or not.  gemsEquilibrium.C
    includes it, so a toggle of LESTO_GEMSBRIDGE recompiles it and relinks
    the executable even with a plain wmake;
  - this header never includes gemsbridge.h.  Only gemsEquilibrium.C does,
    and only under #ifdef LESTO_HAVE_GEMS.

At this development step only the build plumbing and a stub exist:

  - available() and bridgeDirectory() report how the executable was built;
  - buildState() is the one-line summary shown by rhoFixedFlowFoam -help;
  - requireBridge() stops a run that selects the GEMS backend with a
    FatalIOError: a rebuild hint without the bridge, and a "not implemented
    yet" message with it.

The engine itself (one per MPI rank, frozen and local modes, guards, table
fallback) belongs to a later milestone.

------------------------------------------------------------------------------*/

#ifndef gemsEquilibrium_H
#define gemsEquilibrium_H

#include "dictionary.H"
#include "string.H"

namespace LESTO {

class gemsEquilibrium {

public:

  /* true when the executable was built with the GEMS3K bridge */
  static bool available();

  /* bridge directory used at build time; empty without the bridge */
  static Foam::string bridgeDirectory();

  /* one-line summary of the build state */
  static Foam::string buildState();

  /*--------------------------------------------------------------------------
  Stops the run with a FatalIOError reported against dict, the dictionary
  that selected the GEMS backend.  Without the bridge the message explains
  how to rebuild with it; with the bridge it states that the backend is not
  implemented yet.
  --------------------------------------------------------------------------*/
  static void requireBridge(const Foam::dictionary& dict);
};

}

#endif
