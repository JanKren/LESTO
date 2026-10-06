const dictionary& coefficients = speciesDict.subDict("ChapmanEnskogCoeffs");
const scalar carrierM = coefficients.get<scalar>("carrierMolarMass");
const scalar sigma = coefficients.get<scalar>("sigma");
const scalar carrierSigma = coefficients.get<scalar>("carrierSigma");
const scalar epsilon = coefficients.get<scalar>("epsilonOverK");
const scalar carrierEpsilon = coefficients.get<scalar>("carrierEpsilonOverK");
for (const scalar value : {molarMass[speciesi], carrierM, sigma, carrierSigma,
                          epsilon, carrierEpsilon}) {
  if (!std::isfinite(value) || value <= 0) {
    FatalIOErrorInFunction(coefficients)
      << "ChapmanEnskog requires finite positive masses and LJ parameters."
      << exit(FatalIOError);
  }
}
if (speciesDict.found("D")) {
  FatalIOErrorInFunction(speciesDict)
    << "Remove constant D when using ChapmanEnskog." << exit(FatalIOError);
}
const volScalarField& T = thermo.T();
const volScalarField& p = thermo.p();
if (T.dimensions() != dimTemperature || p.dimensions() != dimPressure) {
  FatalErrorInFunction << "Require T [K] and absolute p [Pa]." << exit(FatalError);
}
molecularD.set(speciesi, new volScalarField (
  IOobject("D_" + speciesName, runTime.timeName(), mesh,
           IOobject::NO_READ, frozenWrite),
  T*dimensionedScalar("zeroDPerK", dimViscosity/dimTemperature, 0)
));
auto evaluateD = [&](scalar temperature, scalar pressure, scalar density) {
  if (!std::isfinite(temperature) || !std::isfinite(pressure)
    || !std::isfinite(density) || temperature <= 0 || pressure <= 0
    || density <= 0) {
    FatalErrorInFunction << "ChapmanEnskog requires finite positive T, p, rho."
      << exit(FatalError);
  }
  const scalar D = LESTO::chapmanEnskogDiffusivity (
    temperature, pressure, molarMass[speciesi], carrierM, sigma, carrierSigma,
    epsilon, carrierEpsilon
  );
  if (!std::isfinite(D) || D <= 0) {
    FatalErrorInFunction << "Invalid ChapmanEnskog diffusivity."
      << exit(FatalError);
  }
  return D;
};
forAll(T, celli) {
  molecularD[speciesi][celli] = evaluateD(T[celli], p[celli], rho[celli]);
}
forAll(T.boundaryField(), patchi) {
  forAll(T.boundaryField()[patchi], facei) {
    molecularD[speciesi].boundaryFieldRef()[patchi][facei] = evaluateD (
      T.boundaryField()[patchi][facei], p.boundaryField()[patchi][facei],
      rho.boundaryField()[patchi][facei]
    );
  }
}
rhoD.set(speciesi, new volScalarField (
  IOobject(rhoDName, runTime.timeName(), mesh, IOobject::NO_READ, frozenWrite),
  rho*molecularD[speciesi]
));
Info<< "Species " << speciesName
  << ": diffusivityModel = ChapmanEnskog; dilute binary gas, LJ parameters "
  << "sigma [angstrom], epsilonOverK [K]; D min/max [m2/s] "
  << gMin(molecularD[speciesi].primitiveField()) << " "
  << gMax(molecularD[speciesi].primitiveField()) << nl;
