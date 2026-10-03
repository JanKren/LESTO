/*------------------------------------------------------------------------------

binaryIOList.H

PURPOSE

Lists that OpenFOAM writes whenever it writes the fields, and always in
BINARY, whatever writeFormat says.  rhoFixedFlowFoam keeps the restart state
of its phase change in such lists: the mass ledger in <time>/uniform/, and
the regimes, deposits and reservoirs of the interface elements in
<time>/uniform/phaseChange/ (phaseChangeLedger.H, interfaceExchange.H).
Two kinds:

  binaryIOList<T>        a PER-RANK list (an IOList): every rank holds and
                         writes its own values (the element state of its
                         part of the mesh);
  binaryGlobalIOList<T>  a GLOBAL list (a GlobalIOList): the same values on
                         every rank, read by the master and broadcast (the
                         ledger).

Why lists of their own

  - Exact.  A binary list of a contiguous type (scalar, label) is written as
    a raw block of bytes, so it reads back bit for bit.  Scalars written as
    text are rounded to writePrecision; that happens in an ascii file, and
    also in a binary field file for a field or patch whose values are all
    equal: Field::writeEntry writes it as 'uniform <value>' in text (e.g. a
    rank with a single WALL face).
  - Written with the fields.  The object is registered AUTO_WRITE, so every
    write of Time stores it: the write times of controlDict, and the writes
    that function objects start between two steps with Time::writeNow() or
    Time::writeAndEnd() (e.g. runTimeControl with satisfiedAction end, or
    its nWriteStep).  A list written explicitly at the solver's write times
    only would be missing from those directories, and a restart from them
    would be refused.  Nothing is written unless Time writes (review item
    B4 of doc/phase-change-plan.md).

Per-rank and global objects with the collated file handler.  With
-fileHandler collated a per-rank object is written into ONE file per time
and IO rank, processors<N>/<time>/<local>/<name> (or processors<N>_<lo>-<hi>/
with several IO ranks, -ioRanks), a decomposedBlockData container of one
block per rank of the group; only a run on N ranks with the same IO ranks
can read it.  A global object is written by the master as a plain file
(collatedFileOperation::writeObject: io.global()), which any run can read,
like OpenFOAM's own uniform/time.  reconstructPar and decomposePar copy
uniform/ verbatim, so the ledger, which is the same on every rank, must be
global: as a per-rank list it stopped every serial restart after
reconstructPar ('could not detect processor number') and every restart on
another number of ranks (review of M5, round 3).  The per-rank lists are
read only by a run they can belong to (interfaceExchange::perRankFile()).
Without the collated handler both kinds are written by every rank into its
own processor<i>/<time>/.

A per-rank list must never have the local directory "uniform" itself: the
master-reading file handlers (masterUncollated, and collated, which
derives from it) take every object of uniform/ for global data and
resolve it to the master's file on all ranks
(masterUncollatedFileOperation::filePath(), io.local() == "uniform").
Every rank then read rank 0's list, or the container of rank 0's IO group,
and every restart on the same decomposition lost the exact state (review
of M5, round 4).  They lie in uniform/phaseChange/, which decomposePar and
reconstructPar copy (or link) with uniform/.

The format is forced in writeObject(): Time passes the case's writeFormat
and writeCompression, and the list is written in BINARY with that
compression.  Reading needs nothing special: the header says binary.  The
class name in the header is that of the list (e.g. scalarList, labelList;
GlobalIOList<scalar> has the same one), so the files of both kinds can be
read as ordinary IOLists.

------------------------------------------------------------------------------*/

#ifndef binaryIOList_H
#define binaryIOList_H

#include "IOList.H"
#include "GlobalIOList.H"

namespace LESTO {

template<class T>
class binaryIOList : public Foam::IOList<T> {

public:

  using Foam::IOList<T>::IOList;
  using Foam::IOList<T>::operator=;

  /* the other (deprecated) overload stays visible */
  using Foam::IOList<T>::writeObject;

  /*--------------------------------------------------------------------------
  Called by Time (objectRegistry::writeObject) for every write; streamOpt
  holds the case's writeFormat, which is replaced by BINARY.
  --------------------------------------------------------------------------*/
  virtual bool writeObject (
    Foam::IOstreamOption streamOpt,
    const bool           writeOnProc
  ) const override {
    return Foam::IOList<T>::writeObject (
      Foam::IOstreamOption (
        Foam::IOstreamOption::BINARY,
        streamOpt.compression()
      ),
      writeOnProc
    );
  }
};


/*------------------------------------------------------------------------------
The global kind: identical values on every rank.  GlobalIOList reads in its
own constructor (master only, then broadcast), where its global() and
filePath() already apply; the inherited constructors keep that.
------------------------------------------------------------------------------*/

template<class T>
class binaryGlobalIOList : public Foam::GlobalIOList<T> {

public:

  using Foam::GlobalIOList<T>::GlobalIOList;
  using Foam::GlobalIOList<T>::operator=;

  /* the other (deprecated) overload stays visible */
  using Foam::GlobalIOList<T>::writeObject;

  /* as binaryIOList::writeObject() */
  virtual bool writeObject (
    Foam::IOstreamOption streamOpt,
    const bool           writeOnProc
  ) const override {
    return Foam::GlobalIOList<T>::writeObject (
      Foam::IOstreamOption (
        Foam::IOstreamOption::BINARY,
        streamOpt.compression()
      ),
      writeOnProc
    );
  }
};

}

namespace Foam {

/* a global file type, as GlobalIOList (IOobject::typeHeaderOk) */
template<class T>
struct is_globalIOobject<LESTO::binaryGlobalIOList<T>> : std::true_type {};

}

#endif
