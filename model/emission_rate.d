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
import model.triple_index;
import model.msmc_model;

double emissionProb(string alleles, in MSMCmodel model, double t, size_t ind1, size_t ind2, double tTot) {
  if(alleles.length == 2) {
    auto mu = model.mutationRate;
    if(alleles[0] == alleles[1])
      return exp(-2.0 * mu * t);
    else
      return (1.0 - exp(-2.0 * mu * t));
  }
  
  auto p_0 = emissionProb(alleles, '0', model, t, ind1, ind2);
  auto p_1 = emissionProb(alleles, '1', model, t, ind1, ind2);
  return p_0 + p_1;
}

double emissionProb(string alleles, char ancestralAllele, in MSMCmodel model, double t, size_t ind1, size_t ind2)
in {
  assert(t > 0.0);
  assert(t < double.infinity);
}
out(res) {
  assert(res >= 0.0);
}
body {
  auto M = alleles.length;
  auto mu = model.mutationRate;
  auto count_ancestral = count(alleles, ancestralAllele);
  auto count_derived = M - count_ancestral;
  
  auto tTotUpper = 0.0;
  foreach(k; 1 .. M - 1) {
    tTotUpper += binomial(M - 1, k) * mutationTreeLength(model, t, M - 1, k);
    // stderr.writefln("%s %s %s", k, mutationTreeLength(model, t, M - 1, k), mutationTreeLength(M - 1, k));
  }
  
  auto tTot = M * t + tTotUpper;

  if(count_derived == 0) {
    return exp(-mu * tTot);
  }
  if(count_ancestral == 0) {
    return 0.0;
  }
  if(count_derived == 1) {
    if(alleles[ind1] != alleles[ind2])
      return (1.0 - exp(-mu * tTot)) * t / tTot;
    else
      return (1.0 - exp(-mu * tTot)) * (t + mutationTreeLength(M - 1, 1)) / tTot;
  }
  else {
    if(alleles[ind1] != ancestralAllele && alleles[ind2] != ancestralAllele) {
      // return (1.0 - exp(-mu * tTot)) * mutationTreeLength(model, t, M - 1, count_derived - 1) / tTot;
      auto length = mutationTreeLength(M - 1, 1); // this is the doubleton above the pair of first coalescence
      auto norm = 1.0;
      
      foreach(k; 2 .. M - 1) {
        length += mutationTreeLength( M - 1, k) * binomial(M - 2, k - 1);
        norm += binomial(M - 2, k - 1);
      }
      length /= norm;
          
      return (1.0 - exp(-mu * tTot)) * length / tTot;
    }
    else if(alleles[ind1] == ancestralAllele && alleles[ind2] == ancestralAllele) {
      // return (1.0 - exp(-mu * tTot)) * mutationTreeLength(model, t, M - 1, count_derived) / tTot;
      auto length = 0.0;
      auto norm = 0.0;
      
      foreach(k; 2 .. M - 1) {
        length += mutationTreeLength(M - 1, k) * binomial(M - 2, k);
        norm += binomial(M - 2, k);
      }
      length /= norm;
          
      return (1.0 - exp(-mu * tTot)) * length / tTot;
    }
    else
      return 0.0;
  }
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

double mutationTreeLength(in MSMCmodel model, double start_time, size_t m, size_t k) {
  auto coalescenceTimes = getExpectedCoalescenceTimeIntervals(model, start_time, m);
  auto sum = 0.0;
  foreach(j; 2 .. m - k + 1 + 1) {
    sum += binomial(j, 2) * binomial(m - j, k - 1) * coalescenceTimes[m - j];
  }
  auto result_wakeley_book_eq_4_22 = sum / (binomial(m - 1, k) * k);
  return 2.0 * result_wakeley_book_eq_4_22 / binomial(m, k);
}

unittest {
  auto mu = 0.001;
  auto subpopLabels = new size_t[6];
  auto model = MSMCmodel.withTrivialLambda(mu, mu, subpopLabels, 10, 1);
  foreach(k; 1 .. 4) {
    assert(approxEqual(mutationTreeLength(model, 1.24, 5, k), mutationTreeLength(5, k), 1e-8, 0.0));
  }
}

double[] getExpectedCoalescenceTimeIntervals(in MSMCmodel model, double start_time, size_t M) {
  double[] ret;
  auto last_t = start_time;
  foreach_reverse(m; 2 .. M + 1) {
    auto t = getExpectedCoalescenceTime(last_t, model, m) - last_t;
    ret ~= t;
    last_t += t;
  }
  return ret;
}

double getExpectedCoalescenceTime(double u, in MSMCmodel model, size_t M) {
  auto beta = model.timeIntervals.findIntervalForTime(u);
  auto lambda_b = model.lambda(beta) * binomial(M, 2);
  auto firstTerm = integralHelper(lambda_b, u, model.timeIntervals.rightBoundary(beta));
  
  auto sum = 0.0;
  foreach(gamma; beta + 1 .. model.nrTimeIntervals) {
    auto term = integralHelper(model.lambda(gamma) * binomial(M, 2),
                               model.timeIntervals.leftBoundary(gamma), model.timeIntervals.rightBoundary(gamma));
    auto inner_sum = 0.0;
    if(gamma > beta + 1) {
      foreach(nu; beta + 1 .. gamma)
        inner_sum += model.lambda(nu) * binomial(M, 2) * model.timeIntervals.delta(nu);
    }
    term *= exp(-inner_sum);
    sum += term;
  }
  sum *= exp(-(model.timeIntervals.rightBoundary(beta) - u) * lambda_b);
  return firstTerm + sum;
}

double integralHelper(double lambda, double lower, double upper) {
  if(upper == double.infinity) {
    return 1.0 / lambda + lower;
  }
  return 1.0 / lambda + lower - (1.0 / lambda + upper) * exp(-(upper - lower) * lambda);
}