// SPDX-License-Identifier: GPL-3.0-or-later
#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    argList::addBoolOption("audit", "Check a converged rhoSimpleFoam carrier");
    argList::addBoolOption("project", "Remove small pressure roundoff defects from frozen flux");
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    const instantList times = timeSelector::select0(runTime, args);
    if ((args.found("audit") || args.found("project")) && times.size()) runTime.setTime(times.last(), times.size()-1);
    #include "createMesh.H"
    IOdictionary settings(IOobject("liuCarrierDict", runTime.system(), mesh,
        IOobject::MUST_READ, IOobject::NO_WRITE));
    const scalar target = settings.get<scalar>("massFlowRate");
    if (!(target > 0)) FatalErrorInFunction << "Positive mass flow required" << exit(FatalError);
    if (args.found("audit") || args.found("project"))
    {
        surfaceScalarField phi(IOobject("phi", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
        const label inlet = mesh.boundaryMesh().findPatchID("INLET");
        const label outlet = mesh.boundaryMesh().findPatchID("OUTLET");
        if (inlet < 0 || outlet < 0) FatalErrorInFunction << "Missing flow patches" << exit(FatalError);
        if (args.found("project"))
        {
            const scalarField before(fvc::div(phi)().primitiveField()*mesh.V());
            const scalar originalCellDefect = gMax(mag(before))/target;
            const scalar originalImbalance = mag(gSum(phi.boundaryField()[inlet])
                + gSum(phi.boundaryField()[outlet]))/target;
            // This corrects finite pressure precision, never an unconverged
            // flow. Solve about zero, avoiding subtraction of ~1e5 Pa values.
            if (originalCellDefect > 1e-6 || originalImbalance > 1e-6)
                FatalErrorInFunction << "Carrier not converged enough for projection" << exit(FatalError);
            wordList types(mesh.boundary().size(), "zeroGradient");
            forAll(types, patchi)
                if (mesh.boundaryMesh()[patchi].type() == "wedge"
                    || mesh.boundaryMesh()[patchi].type() == "empty")
                    types[patchi] = mesh.boundaryMesh()[patchi].type();
            types[outlet] = "fixedValue";
            volScalarField correction(IOobject("fluxPotential", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE), mesh,
                dimensionedScalar(dimMass/dimLength/dimTime, 0), types);
            mesh.setFluxRequired(correction.name());
            fvScalarMatrix equation(fvm::laplacian(correction) == fvc::div(phi));
            dictionary controls;
            controls.add("solver", "PCG"); controls.add("preconditioner", "DIC");
            controls.add("tolerance", scalar(1e-12)); controls.add("relTol", scalar(0));
            controls.add("maxIter", label(10000));
            equation.solve(controls);
            surfaceScalarField fluxCorrection(equation.flux());
            const scalar maximumCorrection = gMax(mag(fluxCorrection.primitiveField()))/target;
            if (maximumCorrection > 1e-5)
                FatalErrorInFunction << "Flux projection is too large" << exit(FatalError);
            volVectorField U(IOobject("U", runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
            volScalarField rho(IOobject("rho", runTime.timeName(), mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
            const volVectorField deltaU(fvc::reconstruct(fluxCorrection)/rho);
            const scalar relativeVelocityCorrection = gMax(mag(deltaU.primitiveField()))
                /max(gMax(mag(U.primitiveField())), SMALL);
            if (relativeVelocityCorrection > 1e-5)
                FatalErrorInFunction << "Velocity projection is too large" << exit(FatalError);
            phi -= fluxCorrection;
            U -= deltaU;
            U.correctBoundaryConditions();
            phi.write(); U.write();
            OFstream projection(runTime.path()/"carrier-projection.json");
            projection.precision(17);
            projection << "{\n  \"original_max_cell_flux_defect_over_inflow\": " << originalCellDefect
                << ",\n  \"original_relative_flow_imbalance\": " << originalImbalance
                << ",\n  \"max_internal_face_flux_correction_over_inflow\": " << maximumCorrection
                << ",\n  \"max_velocity_correction_over_max_velocity\": " << relativeVelocityCorrection << "\n}\n";
        }
        const scalar incoming = -gSum(phi.boundaryField()[inlet]);
        const scalar outgoing = gSum(phi.boundaryField()[outlet]);
        const scalarField defect(fvc::div(phi)().primitiveField()*mesh.V());
        const scalar cellDefect = gMax(mag(defect))/target;
        scalar other = 0;
        scalar wedgeNet = 0;
        forAll(phi.boundaryField(), patchi)
            if (patchi != inlet && patchi != outlet)
            {
                if (mesh.boundaryMesh()[patchi].type() == "wedge")
                    wedgeNet += gSum(phi.boundaryField()[patchi]);
                else other += gSum(mag(phi.boundaryField()[patchi]));
            }
        const scalar inletError = mag(incoming/target - 1);
        const scalar flowError = mag(outgoing - incoming)/target;
        // Low-Mach pressure is ~1e5 Pa with a ~Pa drop. The carrier check is
        // separate from the M10e trace-element ledger and linear-solve gates.
        const scalar tolerance = 1e-12;
        const bool passed = inletError <= tolerance && flowError <= tolerance
            && cellDefect <= tolerance && other/target <= tolerance && mag(wedgeNet)/target <= tolerance;
        OFstream report(runTime.path()/"carrier-audit.json");
        report.precision(17);
        report << "{\n  \"status\": \"" << (passed ? "PASS" : "FAIL")
            << "\",\n  \"time\": " << runTime.value()
            << ",\n  \"target_mass_flow_kg_s\": " << target
            << ",\n  \"relative_tolerance\": " << tolerance
            << ",\n  \"inlet_mass_flow_kg_s\": " << incoming
            << ",\n  \"outlet_mass_flow_kg_s\": " << outgoing
            << ",\n  \"relative_inlet_error\": " << inletError
            << ",\n  \"relative_flow_imbalance\": " << flowError
            << ",\n  \"max_cell_flux_defect_over_inflow\": " << cellDefect
            << ",\n  \"wall_and_axis_absolute_flux_over_inflow\": " << other/target
            << ",\n  \"wedge_net_flux_over_inflow\": " << wedgeNet/target << "\n}\n";
        Info << "Carrier audit: " << (passed ? "PASS" : "FAIL")
             << "; max cell flux defect/inflow = " << cellDefect << endl;
        return passed ? 0 : 1;
    }
    IFstream input(runTime.path()/runTime.constant()/"wallTemperature");
    List<Tuple2<scalar, scalar>> profile(input);
    if (profile.size() < 2) FatalErrorInFunction << "Need at least two temperature points" << exit(FatalError);
    forAll(profile, i)
        if (!(profile[i].second() > 0) || (i && profile[i].first() <= profile[i-1].first()))
            FatalErrorInFunction << "Invalid temperature profile" << exit(FatalError);
    auto temperature = [&](scalar x)
    {
        if (x <= profile.first().first()) return profile.first().second();
        for (label i=1; i<profile.size(); ++i)
            if (x <= profile[i].first())
            {
                scalar f = (x-profile[i-1].first())/(profile[i].first()-profile[i-1].first());
                return (1-f)*profile[i-1].second()+f*profile[i].second();
            }
        return profile.last().second();
    };
    const scalar radius = settings.get<scalar>("radius");
    const scalar pressure = settings.get<scalar>("pressure");
    volScalarField T(IOobject("T", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField p(IOobject("p", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volVectorField U(IOobject("U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    const label inlet = mesh.boundaryMesh().findPatchID("INLET");
    const scalar area = gSum(mesh.magSf().boundaryField()[inlet]);
    auto velocity = [&](const vector& c)
    {
        const scalar rho = pressure*.0040026/(8.31446261815324*temperature(c.x()));
        return vector(2*target/(rho*area)*max(scalar(0), 1-sqr(c.y()/radius)), 0, 0);
    };
    forAll(T, celli)
    {
        T[celli] = temperature(mesh.C()[celli].x());
        p[celli] = pressure;
        U[celli] = velocity(mesh.C()[celli]);
    }
    forAll(T.boundaryField(), patchi)
        forAll(T.boundaryField()[patchi], facei)
        {
            const vector& c = mesh.Cf().boundaryField()[patchi][facei];
            T.boundaryFieldRef()[patchi][facei] = temperature(c.x());
            p.boundaryFieldRef()[patchi][facei] = pressure;
            U.boundaryFieldRef()[patchi][facei] = mesh.boundary()[patchi].name() == "WALL"
                ? vector::zero : velocity(c);
        }
    T.write(); p.write(); U.write();
    Info << "Initialised measured-profile boundary values; solve rhoSimpleFoam next" << endl;
    return 0;
}
