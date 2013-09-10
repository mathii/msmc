/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of msmc.
 * msmc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
module model.emission_rate;
import std.stdio;
import std.math;
import std.exception;
import std.conv;
import std.range;
import std.algorithm;
import std.mathspecial;
import model.triple_index;
import model.msmc_model;
import model.time_intervals;
import model.coalescence_rate;

class EmissionRate {
  double[][] upperTreeEmissions;
  const TripleIndex tripleIndex;
  const TimeIntervals timeIntervals;
  const TimeIntervals tTotIntervals;
  const PiecewiseConstantCoalescenceRate coal;
  double mu;
  size_t nrHaplotypes;
  size_t T;
  bool directedEmissions;
  
  this(in TripleIndex tripleIndex, in TimeIntervals timeIntervals, in TimeIntervals tTotIntervals,
       in PiecewiseConstantCoalescenceRate coal, double mu, bool directedEmissions=false)
  {
    this.tripleIndex = tripleIndex;
    this.mu = mu;
    this.timeIntervals = timeIntervals;
    this.tTotIntervals = tTotIntervals;
    this.coal = coal;
    this.directedEmissions = directedEmissions;
    nrHaplotypes = tripleIndex.nrIndividuals;
    T = timeIntervals.nrIntervals;
    if(nrHaplotypes > 2) {
      computeUpperTreeLengths();
    }
  }
  
  private void computeUpperTreeLengths() {
    upperTreeEmissions = new double[][](T, nrHaplotypes - 1);
    foreach(a; 0 .. T) {
      upperTreeEmissions[a][] = 0.0;
      auto t = timeIntervals.meanTimeWithLambda(a, coal.getTotalMarginalLambda(a));
      foreach(k; 1 .. nrHaplotypes - 1) {
        auto val = mutationTreeLength(t, nrHaplotypes - 1, k);
        // auto val = mutationTreeLength(nrHaplotypes - 1, k);
        upperTreeEmissions[a][k] = val;
        upperTreeEmissions[a][0] += binomial(nrHaplotypes - 1, k) * val;
      }
    }
  }
  
  double mutationTreeLength(double start_time, size_t M, size_t k) const {
    auto sum = 0.0;
    foreach(j; 2 .. M - k + 1 + 1) {
      sum += binomial(j, 2) * binomial(M - j, k - 1) * getExpectedCoalescenceTime(start_time, M, j);
    }
    auto result_wakeley_book_eq_4_22 = sum / (binomial(M - 1, k) * k);
    return 2.0 * result_wakeley_book_eq_4_22 / binomial(M, k);
  }

  double getExpectedCoalescenceTime(double u, size_t m, size_t k) const {
    auto beta = timeIntervals.findIntervalForTime(u);
  
    auto firstTerm = 0.0;
    foreach(j; k .. m + 1) {
      auto val = c(m, k, j);
      val *= integralHelper(binomial(j, 2) * coal.getAvgLambda(beta), u, timeIntervals.rightBoundary(beta));
      firstTerm += val;
    }
  
    auto secondTerm = 0.0;
    foreach(gamma; beta + 1 .. T) {
      foreach(j; k .. m + 1) {
        auto inner_sum = 0.0;
        if(gamma - 1 >= beta + 1)
          inner_sum = iota(beta + 1, gamma).map!(nu => coal.getAvgLambda(nu) * timeIntervals.delta(nu)).reduce!"a+b"();
        auto val = exp(-binomial(j, 2) * (timeIntervals.rightBoundary(beta) - u) * coal.getAvgLambda(beta));
        val *= exp(-binomial(j, 2) * inner_sum);
        val *= c(m, k, j);
        val *= integralHelper(binomial(j, 2) * coal.getAvgLambda(gamma), timeIntervals.leftBoundary(gamma), timeIntervals.rightBoundary(gamma));
        secondTerm += val;
      }
    }
    return firstTerm + secondTerm;
  }
  
  double emissionProb(string alleles, size_t aij, size_t i_tTot) const {
    auto triple = tripleIndex.getTripleFromIndex(aij);
    if(nrHaplotypes == 2) {
      auto t = timeIntervals.meanTimeWithLambda(triple.time, coal.getTotalMarginalLambda(triple.time));
      if(alleles[0] == alleles[1])
        return exp(-2.0 * mu * t);
      else
        return (1.0 - exp(-2.0 * mu * t));
    }
    
    auto emissionId = getEmissionId(alleles, triple.ind1, triple.ind2);
    return emissionProb(emissionId, triple.time, i_tTot);
  }
  
  int getEmissionId(string alleles, size_t ind1, size_t ind2) const
  out(res) {
    assert(res == -1 || (res >= 0 && res < getNrEmissionIds));
  }
  body {
    if(directedEmissions)
      return getDirectedEmissionId(alleles, ind1, ind2);
    else
      return getSymmetricEmissionId(alleles, ind1, ind2);
    
  }
  
  private int getDirectedEmissionId(string alleles, size_t ind1, size_t ind2) const {
    auto count_derived = count(alleles, '1');
    auto count_ancestral = nrHaplotypes - count_derived;
    if(count_derived == 0)
      return 0;
    if(count_ancestral == 0)
      return -1;
    if(count_derived == 1) {
      if(alleles[ind1] != alleles[ind2])
        return 1;
      else {
        return cast(int)nrHaplotypes;
      }
    }
    else {
      if(alleles[ind1] != alleles[ind2])
        return -1;
      else {
        if(alleles[ind1] == '1')
          return cast(int)count_derived;
        else
          return cast(int)nrHaplotypes + cast(int)count_derived - 1;
      }
    }
  }
  
  private int getSymmetricEmissionId(string alleles, size_t ind1, size_t ind2) const {
    auto count_0 = count(alleles, '0');
    auto count_1 = nrHaplotypes - count_0;
    if(count_0 == 0 || count_1 == 0)
      return 0;
    if(count_0 == 1 || count_1 == 1) {
      if(alleles[ind1] != alleles[ind2])
        return 1;
      else
        return cast(int)nrHaplotypes - 1;
    }
    else {
      if(alleles[ind1] != alleles[ind2])
        return -1;
      else
        return cast(int)count(alleles, alleles[ind1]);
    }
  }
  
  double emissionProb(int emissionId, size_t timeIndex, size_t i_tTot) const
  in {
    assert(emissionId == -1 || (emissionId >= 0 && emissionId < getNrEmissionIds()), text(emissionId));
    assert(timeIndex < T);
  }
  body {
    if(directedEmissions)
      return directedEmissionProb(emissionId, timeIndex, i_tTot);
    else
      return symmetricEmissionProb(emissionId, timeIndex, i_tTot);
  }
  
  // double symmetricEmissionProb(int emissionId, size_t timeIndex, size_t i_tTot) const
  // {
  //   auto t = timeIntervals.meanTimeWithLambda(timeIndex, coal.getTotalMarginalLambda(timeIndex));
  //   auto tTot = nrHaplotypes * t + upperTreeEmissions[timeIndex][0];
  //   if(emissionId < 0)
  //     return 0.0;
  //   if(emissionId == 0) {
  //     return exp(-mu * tTot);
  //   }
  //   if(emissionId == 1)
  //     return (1.0 - exp(-mu * tTot)) * t / tTot;
  //   else {
  //     auto term = (1.0 - exp(-mu * tTot)) * upperTreeEmissions[timeIndex][emissionId - 1] / tTot;
  //     term += (1.0 - exp(-mu * tTot)) * upperTreeEmissions[timeIndex][nrHaplotypes - emissionId] / tTot;
  //     if(emissionId == nrHaplotypes - 1)
  //       term += (1.0 - exp(-mu * tTot)) * t / tTot;
  //     return term;
  //   }
  // }

  double symmetricEmissionProb(int emissionId, size_t timeIndex, size_t i_tTot) const
  {
    auto t = timeIntervals.meanTimeWithLambda(timeIndex, nrHaplotypes);
    auto tLeaf = tTotIntervals.meanTime(i_tTot, 2);
    if(tLeaf < t * nrHaplotypes)
      tLeaf = t * nrHaplotypes;
    auto tTot = tLeaf + 2.0 * iota(2, nrHaplotypes).map!"1.0 / a"().reduce!"a+b"();
    if(emissionId < 0)
      return 0.0;
    if(emissionId == 0) {
      return exp(-mu * tTot);
    }
    if(emissionId == 1)
      return (1.0 - exp(-mu * tTot)) * t / tTot;
    if(emissionId == nrHaplotypes - 1) {
      auto val = (1.0 - exp(-mu * tTot)) * (tLeaf - 2.0 * t) / (nrHaplotypes - 2.0) / tTot;
      // val += (1.0 - exp(-mu * tTot)) * 2.0 / (nrHaplotypes - 2.0) / (nrHaplotypes - 1.0) / tTot;
      return val;
    }
    else {
      auto val = 2.0 / (nrHaplotypes - 1.0); // doubleton in pair
      auto cnt = 1.0;
      foreach(k; 2 .. nrHaplotypes - 1) {
        val += 2.0 / k;
        cnt += binomial(nrHaplotypes - 1, k);
      }
      return (1.0 - exp(-mu * tTot)) * val / tTot / cnt;
    }
  }
  
  double directedEmissionProb(int emissionId, size_t timeIndex, size_t i_tTot) const {
    auto t = timeIntervals.meanTimeWithLambda(timeIndex, coal.getTotalMarginalLambda(timeIndex));
    auto tTot = nrHaplotypes * t + upperTreeEmissions[timeIndex][0];
    if(emissionId < 0 )
      return 0.0;
    if(emissionId == 0) {
      return exp(-mu * tTot);
    }
    if(emissionId == 1)
      return (1.0 - exp(-mu * tTot)) * t / tTot;
    if(emissionId < nrHaplotypes)
      return (1.0 - exp(-mu * tTot)) * upperTreeEmissions[timeIndex][emissionId - 1] / tTot;
    else {
      auto freq = emissionId - nrHaplotypes + 1;
      if(freq == 1)
        return (1.0 - exp(-mu * tTot)) * (t + upperTreeEmissions[timeIndex][1]) / tTot;
      else
        return (1.0 - exp(-mu * tTot)) * upperTreeEmissions[timeIndex][freq] / tTot;
    }
  }
  
  size_t getNrEmissionIds() const {
    if(directedEmissions)
      return 2 * nrHaplotypes - 2;
    else 
      return nrHaplotypes;
  }
  
}

unittest {
  auto mu = 0.001;
  auto e = MSMCmodel.withTrivialLambda(mu, mu, [0UL, 0, 0, 0, 0, 0], 10, 1, true);
  foreach(k; 1 .. 4) {
    assert(approxEqual(e.emissionRate.mutationTreeLength(1.24, 5, k), mutationTreeLength(5, k), 1e-8, 0.0));
  }
  assert(e.emissionRate.getEmissionId("001100", 1, 2) == -1);
  assert(e.emissionRate.getEmissionId("001100", 2, 3) == 2);
  assert(e.emissionRate.getEmissionId("001100", 0, 1) == 7);
  assert(e.emissionRate.getEmissionId("111111", 0, 1) == -1);
  assert(e.emissionRate.getEmissionId("111110", 0, 1) == 5);

  e = MSMCmodel.withTrivialLambda(mu, mu, [0UL, 0, 0, 0, 0, 0], 10, 1, false);
  assert(e.emissionRate.getEmissionId("001100", 1, 2) == -1);
  assert(e.emissionRate.getEmissionId("001100", 2, 3) == 2);
  assert(e.emissionRate.getEmissionId("001100", 0, 1) == 4);
  assert(e.emissionRate.getEmissionId("111111", 0, 1) == 0);
  assert(e.emissionRate.getEmissionId("111110", 0, 1) == 5);
  
}

double mutationTreeLength(size_t m, size_t freq)
in {
  assert(m >= freq);
  assert(m > 0);
  assert(freq > 0);
}
out(res) {
  assert(res > 0.0);
  assert(res < 1.0);
}
body {
  return 2.0 / freq * (1.0 / binomial(m, freq));
}


double binomial(size_t m, size_t k) {
  auto prod = 1.0;
  foreach(i; 1 .. k + 1) {
    prod *= cast(double)(m - (k - i)) / i;
  }
  return prod;
}

unittest {
  assert(binomial(6,2) == 15);
  assert(binomial(4,2) == 6);
  assert(binomial(6,5) == 6);
  assert(binomial(6,6) == 1);
  assert(binomial(4,4) == 1);
}

double integralHelper(double lambda, double lower, double upper) {
  return (1.0 - exp(-lambda * (upper - lower))) / lambda;
}

double c(size_t n, size_t m, size_t j) {
  auto first = (-1.0) ^^ (j - m);
  auto second = (2.0 * j - 1.0) / cast(double)(fact(m) * fact(j - m));
  auto third = gamma(m + j - 1) / gamma(m);
  auto fourth = gamma(n) / gamma(n + j);
  auto fifth = gamma(n + 1) / gamma(n - j + 1);
  return first * second * third * fourth * fifth;
}

double fact(size_t k) {
  return cast(double).reduce!"a*b"(1UL, iota(1, k + 1));
}

unittest {
  assert(fact(1) == 1);
  assert(fact(2) == 2);
  assert(fact(3) == 6);
}
