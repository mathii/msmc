#!/usr/bin/env rdmd
import std.stdio;
import std.string;
import std.range;
import std.algorithm;
import std.array;
import std.conv;
import model.time_intervals;

auto bottleneck_string = "-eN 0. 5. -eN 0.00210982 3.0303 -eN 0.00427444 1.57203 -eN 0.0064968 1.09639 -eN 0.00878005 0.960283 -eN 0.0111276 0.790608 -eN 0.0135433 0.633806 -eN 0.016031 0.513663 -eN 0.0185953 0.430035 -eN 0.021241 0.374135 -eN 0.0239735 0.297272 -eN 0.0267987 0.297272 -eN 0.0297229 0.247011 -eN 0.0327536 0.247011 -eN 0.0358986 0.244697 -eN 0.039167 0.244697 -eN 0.0425688 0.277391 -eN 0.0461155 0.277391 -eN 0.0498198 0.344184 -eN 0.0536965 0.344184 -eN 0.0577625 0.442999 -eN 0.0620365 0.442999 -eN 0.0665425 0.565825 -eN 0.0713055 0.565825 -eN 0.0763575 0.75083 -eN 0.081736 0.75083 -eN 0.087485 1.01928 -eN 0.093661 1.01928 -eN 0.100331 1.33378 -eN 0.107582 1.33378 -eN 0.115525 1.64472 -eN 0.124305 1.64472 -eN 0.13412 2.01775 -eN 0.178338 2.18741 -eN 0.196522 2.18741 -eN 0.215392 2.11921 -eN 0.235002 2.11921 -eN 0.255413 1.92818 -eN 0.276693 1.92818 -eN 0.298919 1.70894 -eN 0.322179 1.70894 -eN 0.346574 1.51119 -eN 0.37222 1.51119 -eN 0.399254 1.35155 -eN 0.427833 1.35155 -eN 0.458146 1.23164 -eN 0.490415 1.23164 -eN 0.52491 1.14611 -eN 0.561965 1.14611 -eN 0.601985 1.08586 -eN 0.64549 1.08586 -eN 0.693145 1.04152 -eN 0.745825 1.04152 -eN 0.80472 1.00781 -eN 0.871485 1.00781 -eN 0.94856 0.992497 -eN 1.03972 0.992497 -eN 1.1513 1.08308 -eN 1.29514 1.08308 -eN 1.49787 1.7485 -eN 1.84444 1.7485";

void main(string[] args) {
  auto M = args[1].to!size_t();
  auto T = args[2].to!size_t();
  auto fields = bottleneck_string.split();
  auto times = fields.drop(1).stride(3).map!"2.0 * a.to!double()"().array();
  auto lambda = fields.drop(2).stride(3).map!"1.0/a.to!double()"().array();

  auto timeIntervals = TimeIntervals.standardIntervals(T, M);
  
  auto lambdaVec = new double[T];
  
  foreach(i; 0 .. T) {
    auto mean = 0.0;
    auto cnt = 0.0;
    foreach(j; 0 .. times.length) {
      auto a = timeIntervals.roundToFullInterval(times[j]);
      if(a == i) {
        mean += lambda[j];
        cnt += 1.0;
      }
    }
    lambdaVec[i] = mean / cnt;
  }

  lambdaVec.map!(a => text(a)).joiner(",").writeln;
  
}