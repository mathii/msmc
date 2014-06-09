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

module model.transition_rate;
import std.math;
import std.conv;
import std.stdio;
import std.exception;
import std.parallelism;
import std.range;
import model.time_intervals;

class TransitionRate {
  double rho;
  const TimeIntervals timeIntervals;
  const double[] lambdaVec;
  size_t nrTimeIntervals;

  private double[][] transitionProbabilities;
  
  this(in TimeIntervals timeIntervals, double rho, in double[] lambdaVec) {
    enforce(rho > 0.0, "need positive recombination rate");
    this.timeIntervals = timeIntervals;
    this.nrTimeIntervals = timeIntervals.nrIntervals;
    this.rho = rho;
    this.lambdaVec = lambdaVec;
    fillTransitionProbabilitiesParallel();
    // fillTransitionProbabilitiesSingleThread();
  }
  
  private double integrateLambda(double from, double to, size_t fromIndex, size_t toIndex) const {
    double sum = 0.0;
    if(fromIndex == toIndex) {
      return exp(-(to - from) * lambdaVec[toIndex]);
    }
    foreach(kappa; fromIndex + 1 .. toIndex) {
      sum += lambdaVec[kappa] * timeIntervals.delta(kappa);
    }
    
    double ret = exp(-(timeIntervals.rightBoundary(fromIndex) - from) *
               lambdaVec[fromIndex] - sum - (to - timeIntervals.leftBoundary(toIndex)) * lambdaVec[toIndex]);
    return ret;
  }
  
  private void fillTransitionProbabilitiesSingleThread() {
    transitionProbabilities = new double[][](nrTimeIntervals, nrTimeIntervals);
    foreach(b; 0 .. nrTimeIntervals) {
      auto sum = 0.0;
      foreach(a; 0 .. nrTimeIntervals) {
        if(a != b) {
          transitionProbabilities[a][b] = transitionProbabilityOffDiagonal(a, b);
          sum += transitionProbabilities[a][b];
        }
      }
      transitionProbabilities[b][b] = 1.0 - sum;
    }
  }

  private void fillTransitionProbabilitiesParallel() {
    transitionProbabilities = new double[][](nrTimeIntervals, nrTimeIntervals);
    foreach(b; taskPool.parallel(iota(nrTimeIntervals))) {
      auto sum = 0.0;
      foreach(a; 0 .. nrTimeIntervals) {
        if(a != b) {
          transitionProbabilities[a][b] = transitionProbabilityOffDiagonal(a, b);
          sum += transitionProbabilities[a][b];
        }
      }
      transitionProbabilities[b][b] = 1.0 - sum;
    }
  }
  
  private double transitionProbabilityOffDiagonal(size_t a, size_t b) const {
    if(a < b) {
      return q2IntegralSmaller(a, b);
    }
    if(a > b) {
      return q2IntegralGreater(a, b);
    }
    assert(false);
  }
  
  private double q2IntegralSmaller(size_t a, size_t b) const
    in {
      assert(a < b);
    }
  body {
    auto meanTime = timeIntervals.meanTimeWithLambda(b, lambdaVec[b]);
    double integ = (1.0 - exp(-(timeIntervals.delta(a) * lambdaVec[a]))) / lambdaVec[a];
    double sum = 0.0;
    foreach(g; 0 .. a) {
      sum += 2.0 * (1.0 - exp(-timeIntervals.delta(g) * lambdaVec[g])) *
        integrateLambda(timeIntervals.rightBoundary(g), timeIntervals.leftBoundary(a), g + 1, a) / lambdaVec[g] * integ;
    }

    sum += 2.0 * (timeIntervals.delta(a) - integ) / lambdaVec[a];
    
    double ret = (1.0 - exp(-rho * 2 * meanTime)) / (meanTime * 2) * lambdaVec[a] * sum;
    return ret;
  }

  private double q2IntegralGreater(size_t a, size_t b) const
    in {
      assert(a > b);
    }
  body {
    auto meanTime = timeIntervals.meanTimeWithLambda(b, lambdaVec[b]);
    double integ = integrateLambda(meanTime, timeIntervals.leftBoundary(a), b, a) / 
                   lambdaVec[a] * (1.0 - exp(-(timeIntervals.delta(a)) * lambdaVec[a]));
    double sum = 0.0;
    foreach(g; 0 .. b) {
      sum += 2.0 * (1.0 - exp(-lambdaVec[g] * timeIntervals.delta(g))) / lambdaVec[g] * 
             integrateLambda(timeIntervals.rightBoundary(g), meanTime, g + 1, b);
    }
    sum += 2.0 * (1.0 - exp(-lambdaVec[b] * (meanTime - timeIntervals.leftBoundary(b)))) / lambdaVec[b];
      
    return integ * (1.0 - exp(-rho * 2.0 * meanTime)) / (meanTime * 2.0) * lambdaVec[a] * sum;
  }
  
  double transitionProbability(size_t a, size_t b) const {
    return transitionProbabilities[a][b];
  }
  
  double equilibriumProbability(size_t a) const {
    return integrateLambda(0.0, timeIntervals.leftBoundary(a), 0, a) *
        (1.0 - exp(-timeIntervals.delta(a) * lambdaVec[a]));
  }
}

unittest {
  writeln("test transitionRate.equilibriumProbability");
  auto T = 4UL;
  auto r = 0.01;
  auto lambdaVec = [1.0, 0.1, 1, 2];
  auto lvl = 1.0e-8;
  
  auto timeIntervals = TimeIntervals.getQuantileBoundaries(T, 1.0);
  auto transitionRate = new TransitionRate(new TimeIntervals(timeIntervals), r, lambdaVec);
  
  auto sum = 0.0;
  foreach(a; 0 .. T)
    sum += transitionRate.equilibriumProbability(a);
  assert(approxEqual(sum, 1.0, lvl, 0.0), text(sum, " should be 1.0"));
}

unittest {  
  writeln("test transitionRate.transitionProbability");
  auto T = 4UL;
  auto r = 0.01;
  auto lambdaVec = [1, 0.1, 1, 2];
  auto lvl = 1.0e-8;
  
  auto timeIntervals = TimeIntervals.getQuantileBoundaries(T, 1.0);
  auto transitionRate = new TransitionRate(new TimeIntervals(timeIntervals), r, lambdaVec);
  
  foreach(b; 0 .. T) {
    auto sum = 0.0;
    foreach(a; 0 .. T)
      sum += transitionRate.transitionProbability(a, b);
    assert(approxEqual(sum, 1.0, lvl, 0.0), text("transition prob from state ", b, ": ", sum, " should be 1.0"));
  }
}
