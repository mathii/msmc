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
 
module model.psmc_model;
import std.exception;
import std.json;
import std.conv;
import std.file;
import std.stdio;
import std.string;
import std.math;
import std.algorithm;
import model.time_intervals;
import model.transition_rate;

class PSMCmodel {
  const TransitionRate transitionRate;
  const TimeIntervals timeIntervals;
  const double[] lambdaVec;
  double mutationRate;

  this(double mutationRate, double recombinationRate, in double[] lambdaVec, in double[] timeBoundaries) {
    timeIntervals = new TimeIntervals(timeBoundaries ~ [double.infinity]);
    this.mutationRate = mutationRate;
    this.lambdaVec = lambdaVec;
    transitionRate = new TransitionRate(timeIntervals, recombinationRate, lambdaVec);
  }

  override string toString() const {
    return format("<MSMCmodel: mutationRate=%s, recombinationRate=%s, lambdaVec=%s, nrTimeIntervals=%s", mutationRate, recombinationRate, lambdaVec, timeIntervals.nrIntervals);
  }
  
  static PSMCmodel withTrivialLambda(double mutationRate, double recombinationRate, size_t nrTimeIntervals) {
    auto lambdaVec = new double[nrTimeIntervals];
    lambdaVec[] = 1.0;
    auto boundaries = TimeIntervals.getQuantileBoundaries(nrTimeIntervals, 1.0);
    return new PSMCmodel(mutationRate, recombinationRate, lambdaVec, boundaries[0 .. $ - 1]);
  }
  
  double emissionProb(size_t id, size_t timeIndex) const {
    auto t = timeIntervals.meanTimeWithLambda(timeIndex, lambdaVec[timeIndex]);
    switch(id) {
      case 0:
      return 1.0;
      case 1:
      return exp(-2.0 * mutationRate * t);
      case 2:
      return 1.0 - exp(-2.0 * mutationRate * t);
      default:
      assert(false);
    }
  }
  
  double emissionProb(char[2] alleles, size_t timeIndex) const {
    if("ACTG01".canFind(alleles[0]) && "ACTG01".canFind(alleles[1])) {
      if(alleles[0] == alleles[1])
        return emissionProb(1UL, timeIndex);
      else
        return emissionProb(2UL, timeIndex);
    }
    else
      return emissionProb(0UL, timeIndex);
  }
  
  double transitionProb(size_t a, size_t b) const {
    return transitionRate.transitionProbability(a, b);
  }
  
  double equilibriumProb(size_t a) const {
    return transitionRate.equilibriumProbability(a);
  }

  @property size_t nrStates() const {
    return timeIntervals.nrIntervals;
  }
  
  @property double recombinationRate() const {
    return transitionRate.rho;
  }
    
}
