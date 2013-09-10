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
import std.getopt;
import std.exception;
import std.c.stdlib;
import std.algorithm;
import std.parallelism;
import std.string;
import std.typecons;
import std.regex;
import std.conv;
import std.range;
import std.math;
import std.functional;
import model.msmc_model;
import model.emission_rate;
import model.data;
import model.triple_index;

auto bottleneck_string = "-eN 0. 5. -eN 0.00210982 3.0303 -eN 0.00427444 1.57203 -eN 0.0064968 1.09639 -eN 0.00878005 0.960283 -eN 0.0111276 0.790608 -eN 0.0135433 0.633806 -eN 0.016031 0.513663 -eN 0.0185953 0.430035 -eN 0.021241 0.374135 -eN 0.0239735 0.297272 -eN 0.0267987 0.297272 -eN 0.0297229 0.247011 -eN 0.0327536 0.247011 -eN 0.0358986 0.244697 -eN 0.039167 0.244697 -eN 0.0425688 0.277391 -eN 0.0461155 0.277391 -eN 0.0498198 0.344184 -eN 0.0536965 0.344184 -eN 0.0577625 0.442999 -eN 0.0620365 0.442999 -eN 0.0665425 0.565825 -eN 0.0713055 0.565825 -eN 0.0763575 0.75083 -eN 0.081736 0.75083 -eN 0.087485 1.01928 -eN 0.093661 1.01928 -eN 0.100331 1.33378 -eN 0.107582 1.33378 -eN 0.115525 1.64472 -eN 0.124305 1.64472 -eN 0.13412 2.01775 -eN 0.178338 2.18741 -eN 0.196522 2.18741 -eN 0.215392 2.11921 -eN 0.235002 2.11921 -eN 0.255413 1.92818 -eN 0.276693 1.92818 -eN 0.298919 1.70894 -eN 0.322179 1.70894 -eN 0.346574 1.51119 -eN 0.37222 1.51119 -eN 0.399254 1.35155 -eN 0.427833 1.35155 -eN 0.458146 1.23164 -eN 0.490415 1.23164 -eN 0.52491 1.14611 -eN 0.561965 1.14611 -eN 0.601985 1.08586 -eN 0.64549 1.08586 -eN 0.693145 1.04152 -eN 0.745825 1.04152 -eN 0.80472 1.00781 -eN 0.871485 1.00781 -eN 0.94856 0.992497 -eN 1.03972 0.992497 -eN 1.1513 1.08308 -eN 1.29514 1.08308 -eN 1.49787 1.7485 -eN 1.84444 1.7485";


double mutationRate, recombinationRate;
size_t nrTimeSegments=20;
string[] inputFileNames;
size_t nrHaplotypes;
string[] treeFileNames;

void main(string[] args) {
  try {
    parseCommandlineArgs(args);
  }
  catch (Exception e) {
    stderr.writeln(e.msg);
    displayHelpMessage();
    exit(0);
  }
  run();
}

void parseCommandlineArgs(string[] args) {
  getopt(args,
         std.getopt.config.caseSensitive,
         "mutationRate|m", &mutationRate,
         "nrTimeSegments|T", &nrTimeSegments,
         "sites_file|s", &inputFileNames,
         "tree_file|t", &treeFileNames);

  enforce(inputFileNames.length == treeFileNames.length);
  nrHaplotypes = getNrHaplotypesFromFile(inputFileNames[0]);
  enforce(mutationRate > 0, "need positive mutationRate");
}

void displayHelpMessage() {
  stderr.writeln("Usage: build/emissionHistogram [options]
Options:
-m, --mutationRate <double>
-T, --nrTimeSegments <int>
-s, --sites_file
-t, --tree_file");
}

void run() {
  
  auto subpopLabels = new size_t[nrHaplotypes];
  auto fields = bottleneck_string.split();

  auto lambdaVec = [0.265,0.963603,1.94999,3.13356,4.06754,3.37182,2.47337,1.5496,0.903973,0.655254,0.495602,0.749641];
  
  // auto times = fields.drop(1).stride(3).map!"2.0 * a.to!double()"().array();
  // auto lambda = fields.drop(2).stride(3).map!"1.0 / a.to!double()"().array();
  // auto model_th = new MSMCmodel(mutationRate, mutationRate / 2.0, subpopLabels, lambda, times, 1, true);
  auto model_th = new MSMCmodel(mutationRate, mutationRate / 2.0, subpopLabels, lambdaVec, nrTimeSegments, 1, true);
  auto model = MSMCmodel.withTrivialLambda(mutationRate, mutationRate / 2.0, subpopLabels, nrTimeSegments, 1, true);
  
  auto nullString = new char[nrHaplotypes];
  nullString[] = '0';

  auto counts = new size_t[string][model.nrStates];
  auto norm = new size_t[model.nrStates];

  writeln("alleles\ttype\tt\ttIndex\tstate\tprob\tprob_th");

  foreach(i; 0 .. inputFileNames.length) {
    auto simStateParser = new SimStateParser(treeFileNames[i], nrTimeSegments, nrTimeSegments);
  
    auto inputFile = File(inputFileNames[i], "r");
    foreach(line; inputFile.byLine) {
      auto line_fields = line.strip().split();
      auto pos = line_fields[1].to!size_t();
      auto called_sites = line_fields[2].to!size_t();
      auto alleles = line_fields[3].idup;

      auto tPair = simStateParser.getFirstPair(pos);
      auto tIndex = model.timeIntervals.findIntervalForTime(tPair[0]);
      auto tState = model.marginalIndex.getIndexFromTriple(Triple(tIndex, tPair[1], tPair[2]));
    
      counts[tState][alleles] += 1;
      counts[tState][nullString.idup] += called_sites - 1;
      norm[tState] += called_sites;
    }
  }
    
  foreach(aij; 0 .. model.nrStates) {
    foreach(al, cnt; counts[aij]) {
      auto tState = model.marginalIndex.getTripleFromIndex(aij);
      auto tIndex = tState.time;
      auto t = model.timeIntervals.meanTime(tIndex, nrHaplotypes);
      auto prob = cast(double)cnt / norm[aij];
      auto eType = model_th.emissionRate.getEmissionId(al, tState.ind1, tState.ind2);
      auto tIndex_th = model_th.timeIntervals.findIntervalForTime(t);
      auto prob_th = model_th.emissionRate.emissionProb(eType, tIndex_th);
      writefln("%s\t%s\t%s\t%s\t%s\t%s\t%s", al, eType, t, tIndex, aij, prob, prob_th);
    }
  }
  
}

class SimStateParser {
  
  Tuple!(size_t, double, size_t, size_t)[] data;
  size_t lastIndex;
  
  this(string treeFileName, size_t nrTimeSegments, size_t nrTtotSegments) {
    
    stderr.writeln("reading tree file ", treeFileName);
    auto treeFile = File(treeFileName, "r");
    auto pos = 0UL;
    foreach(line; treeFile.byLine) {
      auto fields = line.strip().split();
      auto l = fields[0].to!size_t;
      auto str = fields[1];
      auto tFirst = getFirst(str);
      pos += l;
      data ~= tuple(pos, tFirst[0], tFirst[1], tFirst[2]);
    }
  }
  
  auto getFirstPair(size_t pos) {
    auto index = getIndex(pos);
    auto a = data[index][1];
    auto i = data[index][2];
    auto j = data[index][3];
    return tuple(a, i, j);
  }
  
  private size_t getIndex(size_t pos) {
    while(data[lastIndex][0] < pos)
      lastIndex += 1;
    while(lastIndex > 0 && data[lastIndex - 1][0] >= pos)
      lastIndex -= 1;
    return lastIndex;
  }
  
}

private auto getFirst(in char[] str) {
  static auto tfirstRegex = regex(r"\((\d+):([\d\.e-]+),(\d+):[\d\.e-]+\)", "g");
  
  auto matches = match(str, tfirstRegex);
  auto triples = matches.map!(m => tuple(2.0 * m.captures[2].to!double(),
                                         m.captures[1].to!size_t(),
                                         m.captures[3].to!size_t()))();
  auto min_triple = triples.minCount!"a[0] < b[0]"()[0];
  if(min_triple[1] > min_triple[2])
    swap(min_triple[1], min_triple[2]);
  return min_triple;
}

version(unittest) {
  import std.math;
} 

unittest {

  auto tree = "(5:0.1,((2:8.1,1:8.1):0.1,((4:0.1,0:0.1):0.004,3:0.1):0.004):0.01);";
  auto f = getFirst(tree);
  assert(approxEqual(f[0], 2.0 * 0.1, 0.00001, 0.0));
  assert(f[1] == 0);
  assert(f[2] == 4);
  tree = "((((2:8.3,1:8.3):0.12,(0:0.11,3:0.11):0.004):0.0046,4:0.12):1.06,5:1.19);";
  f = getFirst(tree);
  assert(approxEqual(f[0], 2.0 * 0.11, 0.00001, 0.0));
  assert(f[1] == 0);
  assert(f[2] == 3);
}

string getEmissionType(string alleles, size_t ind1, size_t ind2) {
  auto count_derived = count(alleles, '1');
  if(count_derived == 1) {
    if(alleles[ind1] == alleles[ind2])
      return "out_1";
    else
      return "in_1";
  }
  if(count_derived > 1 && alleles[ind1] == alleles[ind2]) {
    if(alleles[ind1] == '1')
      return format("in_%s", count_derived);
    else
      return format("out_%s", count_derived);
  }
  return "weird";
}