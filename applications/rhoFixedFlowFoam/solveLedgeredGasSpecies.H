// Zero-exchange counterpart of solvePairedGasSpecies.H. The same conservative
// unrelaxed transport equation, with a final solve and flux-form defect.
volScalarField& Y=species[speciesi];
respeciation.checkTransport(speciesi);
const word divScheme("div(phi,"+Y.name()+")");
const word finalName(Y.name()+"Final");
const bool haveFinal=mesh.solution().solversDict().found(finalName);
scalar residual=GREAT;
for(label corr=0;corr<=nCorr;++corr) {
  fvScalarMatrix eqn(fvm::ddt(rho,Y)+fvm::div(phi,Y,divScheme)-fvm::laplacian(rhoD[speciesi],Y));
  residual=eqn.solve(mesh.solution().solverDict(Y.name())).initialResidual();
  if(residual<tolerance || corr==nCorr) {
    if(!haveFinal)respeciation.recordSolve(speciesi,eqn);
    break;
  }
}
if(haveFinal) {
  fvScalarMatrix eqn(fvm::ddt(rho,Y)+fvm::div(phi,Y,divScheme)-fvm::laplacian(rhoD[speciesi],Y));
  eqn.solve(mesh.solution().solverDict(finalName));respeciation.recordSolve(speciesi,eqn);
}
Info<<Y.name()<<" min/max = "<<gMin(Y.primitiveField())<<' '<<gMax(Y.primitiveField())<<nl;
