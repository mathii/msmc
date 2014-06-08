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
import model.triple_index_marginal;
import model.time_intervals;
import model.emission_rate;
import model.transition_rate;
import model.coalescence_rate;

class PSMCmodel {
  const TransitionRate transitionRate;
  const TimeIntervals timeIntervals;
  const double[] lambdaVec;

  this(double mutationRate, double recombinationRate, in double[] lambdaVec, in double[] timeBoundaries) {
    timeIntervals = new TimeIntervals(timeBoundaries ~ [double.infinity]);
    this.lambdaVec = lambdaVec;
    transitionRate = new TransitionRate(marginalIndex, coal, timeIntervals, recombinationRate);
  }

  override string toString() const {
    return format("<MSMCmodel: mutationRate=%s, recombinationRate=%s, lambdaVec=%s, nrTimeIntervals=%s", mutationRate, recombinationRate, lambdaVec, nrTimeIntervals);
  }
  
  static MSMCmodel withTrivialLambda(double mutationRate, double recombinationRate, size_t nrTimeIntervals) {
    auto lambdaVec = double[nrTimeIntervals];
    lambdaVec[] = 1.0;
    return new MSMCmodel(mutationRate, recombinationRate, lambdaVec, nrTimeIntervals);
  }
  
  double emissionProb(char[2] alleles, size_t timeIndex) {
    if(alleles[0] == alleles[1])
      return exp(-2.0 * mu * t);
    else
      return 1.0 - exp(-2.0 * mu * t);
  }
  
  double transitionProb(size_t a, size_t b) {
    return transitionRate.transitionProbability(a, b);
  }
  
  double equilibriumProb(size_t a) {
    return transitionRate.equilibriumProbability(a);
  }

  @property size_t nrStates() const {
    return marginalIndex.nrStates;
  }
  
  @property double mutationRate() const {
    return emissionRate.mu;
  }
  
  @property double recombinationRate() const {
    return transitionRate.rho;
  }
  
  @property double[] lambdaVec() const {
    return lambdaVec.dup;
  }
  
}
