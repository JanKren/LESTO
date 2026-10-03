// End-to-end check: open-data PbI2-He GEMS3K system -> C bridge -> OpenFOAM.
//   gemsPbI2Check <PbI2He-dat.lst>
#include "fvCFD.H"
#include "gemsbridge.h"
#include <chrono>
#include <vector>

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::validArgs.append("system-dat.lst");
    argList args(argc, argv, false, false);

    gemsb_set_log_level(4);
    char err[256] = "";
    gemsb_engine* e = gemsb_create_from_lst(args.get<fileName>(1).c_str(), err, 256);
    if (!e) FatalErrorInFunction << "create failed: " << err << exit(FatalError);

    const int nIC = gemsb_num_elements(e), nDC = gemsb_num_species(e), nPH = gemsb_num_phases(e);
    const int iHe = gemsb_element_index(e, "He"), iI = gemsb_element_index(e, "I"),
              iPb = gemsb_element_index(e, "Pb");
    const int jPbI2 = gemsb_species_index(e, "PbI2(g)");
    double Tmin, Tmax, Pmin, Pmax;
    gemsb_TP_range(e, &Tmin, &Tmax, &Pmin, &Pmax);
    Info<< "system: " << nIC << " elements, " << nDC << " species, " << nPH
        << " phases; T grid " << Tmin << "-" << Tmax << " K" << nl;
    for (int k = 0; k < nPH; ++k) Info<< "  phase " << k << ": " << word(gemsb_phase_name(e, k)) << nl;

    const double P = 101325.0;
    std::vector<double> b(nIC, 0.0), p(nDC), nph(nPH);

    // 1) excess PbI2 (10 mol per mol He): condensate present, p(PbI2) = p_v(T)
    Info<< nl << "T[K]   log10 p(PbI2)[Pa]   ref(Gurvich)   diff   condensed phases" << nl;
    const double Tref[]  = {500, 600, 700, 800, 900, 1000, 1100};
    const double lgRef[] = {-3.1532, -0.3623, 1.5553, 2.7683, 3.6747, 4.3701, 4.9150};
    for (int t = 0; t < 7; ++t)
    {
        b[iHe] = 1.0; b[iPb] = 10.0; b[iI] = 20.0*(1.0 + 1e-6);
        int rc = gemsb_equilibrate(e, Tref[t], P, b.data(), GEMSB_COLD, nullptr, nullptr);
        gemsb_gas_partial_pressures(e, p.data());
        gemsb_phase_amounts(e, nph.data());
        const double lg = std::log10(p[jPbI2]);
        Info<< Tref[t] << "   " << lg << "   " << lgRef[t] << "   " << lg - lgRef[t] << "   rc=" << rc;
        for (int k = 0; k < nPH; ++k) if (nph[k] > 0 && word(gemsb_phase_name(e, k)).find("PbI2") != std::string::npos)
            Info<< "  " << word(gemsb_phase_name(e, k)) << "=" << nph[k];
        Info<< nl;
    }

    // 2) CFD-like pattern: nCell wall cells on a 1100 -> 400 K profile,
    //    dilute PbI2 slowly varying in time, per-cell warm start
    const int nCell = 1000, nStep = 50, sz = gemsb_state_size(e);
    std::vector<double> state(size_t(nCell)*sz);
    std::vector<int> valid(nCell, 0);
    int nFail = 0, nRetry = 0; long iters = 0; long calls = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < nStep; ++s)
    {
        for (int c = 0; c < nCell; ++c)
        {
            const double T = 1100.0 - 700.0*c/(nCell - 1);
            const double xPb = 1e-3*(1.0 + 0.5*std::sin(0.1*s + 0.01*c));
            b[iHe] = 1.0; b[iPb] = xPb; b[iI] = 2.0*xPb*(1.0 + 1e-6);
            int rc = gemsb_equilibrate(e, T, P, b.data(), GEMSB_WARM_CELL,
                                       &state[size_t(c)*sz], &valid[c]);
            if (rc == GEMSB_OK_RETRIED) ++nRetry;
            else if (rc != GEMSB_OK) ++nFail;
            iters += gemsb_last_iterations(e); ++calls;
        }
    }
    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Info<< nl << "CFD-like loop: " << calls << " calls, " << 1e6*dt/calls << " us/call, "
        << double(iters)/calls << " iterations/call, cold retries " << nRetry
        << ", failures " << nFail << nl;

    gemsb_destroy(e);
    Info<< "End" << nl;
    return 0;
}
