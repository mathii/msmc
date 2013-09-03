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
 
import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.getopt;
import std.parallelism;
import std.algorithm;
import std.array;
import std.json;
import std.file;
import std.typecons;
import std.exception;
import std.c.stdlib;
import core.memory;
import model.msmc_hmm;
import model.msmc_model;
import model.triple_index_marginal;
import model.emission_rate;
import model.transition_rate;
import model.time_intervals;
import model.triple_index_marginal;
import model.coalescence_rate;
import model.rate_integrator;
import model.propagation_core_fastImpl;
import model.data;

void estimateTotalBranchlengths(SegSite_t[] inputData, MSMCmodel params, size_t internalNrSegments) {

  auto propagationCore = buildPropagationCore(params, internalNrSegments);
  auto msmc_hmm = buildHMM(inputData, propagationCore);
    
  msmc_hmm.runForward();
  
  auto forwardState = propagationCore.newForwardState();
  auto backwardState = propagationCore.newBackwardState();

  foreach_reverse(dataIndex; 0 .. inputData.length) {
    msmc_hmm.getForwardState(forwardState, inputData[dataIndex].pos);
    msmc_hmm.getBackwardState(backwardState, inputData[dataIndex].pos);
    double ttot = 2.0 * propagationCore.msmc.timeIntervals.meanTimeWithLambda(0, 1.0);
    auto max = forwardState.vec[0] * backwardState.vec[0];
    foreach(i; 0 .. propagationCore.msmc.nrTimeIntervals) {
      auto p = forwardState.vec[i] * backwardState.vec[i];
      if(p > max) {
        max = p;
        // we need the total branch length, so twice the tMRCA with two haplotypes
        ttot = 2.0 * propagationCore.msmc.timeIntervals.meanTimeWithLambda(i, 1.0);
      }
    }
    inputData[dataIndex].i_Ttot = params.tTotIntervals.findIntervalForTime(ttot);
  }
}

void readTotalBranchlengths(SegSite_t[] inputData, MSMCmodel params, size_t internalNrSegments, string treeFileName) {
  auto simTreeParser = new SimTreeParser(treeFileName, internalNrSegments);
  foreach(ref segsite; inputData) {
    auto t = simTreeParser.getTtot(segsite.pos);
    segsite.i_Ttot = params.tTotIntervals.findIntervalForTime(t);
  }
}

class SimTreeParser {
  
  Tuple!(size_t, double)[] data;
  size_t lastIndex;
  
  this(string treeFileName, size_t nrTtotSegments) {
    
    auto treeFile = File(treeFileName, "r");
    auto pos = 0UL;
    foreach(line; treeFile.byLine) {
      auto fields = line.strip().split();
      auto l = fields[0].to!size_t;
      auto str = fields[1];
      auto tTot = getTotLeafLength(str);
      pos += l;
      data ~= tuple(pos, tTot);
    }
  }
  
  double getTtot(size_t pos) {
    auto index = getIndex(pos);
    return data[index][1];
  }
  
  private size_t getIndex(size_t pos) {
    while(data[lastIndex][0] < pos)
      lastIndex += 1;
    while(lastIndex > 0 && data[lastIndex - 1][0] >= pos)
      lastIndex -= 1;
    return lastIndex;
  }
}

double getTotLeafLength(in char[] str) {
  static auto tTotRegex = regex(r"\d+:([\d+\.e-]+)", "g");

  auto matches = match(str, tTotRegex);
  auto times = matches.map!(m => m.captures[1].to!double());
  auto sum = 2.0 * times.reduce!"a+b"();
  return sum;
}

unittest {
  auto tree = "(5:0.1,((2:8.1,1:8.1):0.1,((4:0.1,0:0.1):0.004,3:0.1):0.004):0.01);";
  assert(approxEqual(getTotLeafLength(tree), 2.0 * 16.6, 0.0001, 0.0));
  tree = "((((2:8.3,1:8.3):0.122683,(0:0.11,3:0.11):0.00415405):0.00462688,4:0.12):1.06837,5:1.19);";
  assert(approxEqual(getTotLeafLength(tree), 2.0 * 18.13, 0.0001, 0.0));
}




  
// private PropagationCoreFast buildPropagationCore(MSMCmodel params, size_t internalNrSegments) {
//   auto lambdaVec = new double[internalNrSegments];
//   lambdaVec[] = 1.0;
//   // the factor 2 is just part of the formula for the mean total branch length.
//   // auto expectedTtot = 2.0 * TimeIntervals.computeWattersonFactor(params.nrHaplotypes);
//   auto expectedTtot = 2.0 / (params.nrHaplotypes - 1.0);
//   // the next factor 2 fakes a two haplotype system with the same total branch length (every branch gets half)
//   auto boundaries = TimeIntervals.getQuantileBoundaries(internalNrSegments, expectedTtot / 2.0);
//   auto model = new MSMCmodel(params.mutationRate, params.recombinationRate, [0UL, 0], lambdaVec, boundaries[0 .. $ - 1], 1);
// 
//   auto propagationCore = new PropagationCoreFast(model, 1000);
//   return propagationCore;
// }
//   
// private MSMC_hmm buildHMM(SegSite_t[] inputData, PropagationCoreFast propagationCore) {
//   SegSite_t[] dummyInputData;
//   auto allele_order = canonicalAlleleOrder(propagationCore.msmc.nrHaplotypes);
//   foreach(s; inputData) {
//     auto alleles = allele_order[s.obs[0] - 1];
//     auto counts_0 = count(alleles, '0');
//     if(counts_0 == 1 || alleles.length - counts_0 == 1) {
//       auto dummySite = s.dup;
//       dummySite.obs = [2];
//       dummyInputData ~= dummySite;
//     }
//   }
//     
//   return new MSMC_hmm(propagationCore, dummyInputData);
// }
// 
