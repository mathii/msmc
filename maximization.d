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
import std.concurrency;
import std.algorithm;
import std.array;
import std.json;
import std.file;
import std.typecons;
import std.exception;
import std.c.stdlib;
import core.memory;
import msmc_hmm;
import utils;
import msmc_utils;
import msmc_model;
import triple_index_marginal;
import maximization_step;

class MaximizationApplication {
  
  double mutationRate;
  size_t[] subpopLabels;
  size_t[] timeSegmentPattern;
  
  double[][] expectationResult;
  
  string subpopLabelsString;
  string timeSegmentPatternString;
  string inputFileName;
  
  static void runWithArgs(string[] args) {
    auto app = new MaximizationApplication(args);
    app.run();
  }
  
  this(string[] args) {
    timeSegmentPatternString = "1*4+13*2+1*10";
    parseCommandLine(args);
  }

  void parseCommandLine(string[] args) {
    if(args.length == 1)
      displayHelpMessageAndExit();
    try {
      readArguments(args);
      subpopLabels = parseCommaSeparatedArray(subpopLabelsString);
      timeSegmentPattern = parseTimeSegmentPattern(timeSegmentPatternString);
    }
    catch(Exception e) {
      displayHelpMessageAndExit(e);
    }
  }
  
  void readArguments(string[] args) {
    getopt(args,
        std.getopt.config.caseSensitive,
        "mutationRate|m", &mutationRate,
        "subpopLabels|P", &subpopLabelsString,
        "timeSegmentPattern|p", &timeSegmentPatternString
    );
    enforce(args.length == 2, "need exactly one input file");
    inputFileName = args[1];
  }
  
  static void displayHelpMessageAndExit() {
    stderr.writeln("Usage: msmc maximization [options] datafiles
-p, --timeSegmentPattern=<string> : pattern of fixed segments [=1*4+13*2+1*10]
-m, --mutationRate=<double> : scaled mutation rate to use
-P, --subpopLabels=<string> comma-separated subpopulation labels");
    exit(0);
  }
  
  static void displayHelpMessageAndExit(Exception e) {
    writeln(e.msg);
    displayHelpMessageAndExit();
  }
  
  void run() {
    auto nrTimeSegments = reduce!"a+b"(timeSegmentPattern);
    auto modelParams = MSMCmodel.withTrivialLambda(mutationRate, mutationRate / 2.0, subpopLabels, nrTimeSegments, 10);
    
    readMatrixFromStdIn();
    
    MaximizationStep maximizationStep;
    try {
      maximizationStep = new MaximizationStep(expectationResult, modelParams, timeSegmentPattern);
    }
    catch(Exception e) {
      displayHelpMessageAndExit(e);
    }
    
    maximizationStep.run(true);
  }
  
  void readMatrixFromStdIn() {
    auto inputFile = File(inputFileName, "r");
    expectationResult.length = 0;
    foreach(line; inputFile.byLine) {
      auto row = map!"to!double(a)"(line.strip().split()).array;
      expectationResult ~= row;
    }
  }
}  
