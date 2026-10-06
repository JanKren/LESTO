// M10d: the same bounded predictor/guard and fixed transport matrix as pairs.
volScalarField& Y=species[speciesi];
respeciation.checkTransport(speciesi);
LESTO::wallChannels& wall=respeciation.walls();
wall.prepare(speciesi);
const word divScheme("div(phi,"+Y.name()+")"),finalName(Y.name()+"Final");
const bool haveFinal=mesh.solution().solversDict().found(finalName);
autoPtr<fvScalarMatrix> transport;
scalar residual=GREAT;
for(label corr=0;corr<=nCorr+1;++corr){
  fvScalarMatrix eqn(fvm::ddt(rho,Y)+fvm::div(phi,Y,divScheme)-fvm::laplacian(rhoD[speciesi],Y));
  wall.linearise(speciesi,eqn);transport.reset(new fvScalarMatrix(eqn));
  eqn+=fvm::Sp(wall.Sp(speciesi),Y);eqn-=wall.Su(speciesi);
  const bool final=corr==nCorr+1 || residual<tolerance;
  residual=eqn.solve(mesh.solution().solverDict(final && haveFinal?finalName:Y.name())).initialResidual();
  if(final)break;
}
for(label guard=0;guard<8 && wall.guard(speciesi);++guard){
  fvScalarMatrix eqn(transport());eqn+=fvm::Sp(wall.Sp(speciesi),Y);eqn-=wall.Su(speciesi);
  eqn.solve(mesh.solution().solverDict(haveFinal?finalName:Y.name()));
}
wall.realise(speciesi);
respeciation.recordSolve(speciesi,transport());
Info<<Y.name()<<" min/max = "<<gMin(Y.primitiveField())<<' '<<gMax(Y.primitiveField())<<nl;
