import std.stdio;
import std.random;
import model.msmc_model;
import model.data;

void main(string[] args) {
  
  auto mu = 0.0003578;
  auto r = 0.0001;
  auto l = [0.265,0.963603,1.94999,3.13356,4.06754,3.37182,2.47337,1.5496,0.903973,0.655254,0.495602,0.749641];
  auto T = 12UL;
  auto M = 4UL;
  auto L = 10_000_000UL;
  auto subpopLabels = new size_t[M];
  
  auto model = new MSMCmodel(mu, r, subpopLabels, l, T, 1, true);
  auto alleles = canonicalAlleleOrder(M);
  auto nrObs = alleles.length;
  auto nrS = model.nrStates;
  
  auto transitionMatrix = new double[][](nrS, nrS);
  auto emissionMatrix = new double[][](nrS, nrObs);
  auto equilibriumProb = new double[nrS];
  
  foreach(aij; 0 .. nrS) {
    foreach(bkl; 0 .. nrS) {
      transitionMatrix[bkl][aij] = model.transitionRate.transitionProbability(aij, bkl);
    }
    foreach(o; 0 .. nrObs)
      emissionMatrix[aij][o] = model.emissionRate.emissionProb(alleles[o], aij);
    equilibriumProb[aij] = model.transitionRate.equilibriumProbability(aij);
  }
  
  auto state = dice(equilibriumProb);
  auto called_sites = 0UL;
  foreach(pos; 0UL .. L) {
    auto obs = dice(emissionMatrix[state]);
    called_sites += 1;
    if(obs > 0) {
      writefln("chr_sim\t%s\t%s\t%s\t%s", pos, called_sites, alleles[obs], state);
      called_sites = 0;
    }
    
    auto newState = dice(transitionMatrix[state]);
    state = newState;
  }
}