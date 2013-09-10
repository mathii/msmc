#!/usr/bin/env rdmd

import std.stdio;
import std.regex;
import std.algorithm;
import std.conv;
import std.typecons;

void main() {
  auto sum = 0.0;
  auto cnt = 0.0;
  auto sum_first = 0.0;
  foreach(line; stdin.byLine) {
    sum += getTotalBranchLength(line);
    sum_first += getFirst(line);
    cnt += 1.0;
  }
  writeln(sum_first / cnt);
  writeln(sum / cnt);
}

double getTotalBranchLength(in char[] str) {
  static auto tTotRegex = regex(r":([\d+\.e-]+)", "g");

  auto matches = match(str, tTotRegex);
  auto times = matches.map!(m => m.captures[1].to!double());
  auto sum = 2.0 * times.reduce!"a+b"();
  return sum;
}

double getFirst(in char[] str) {
  static auto tfirstRegex = regex(r"\((\d+):([\d\.e-]+),(\d+):[\d\.e-]+\)", "g");
  
  auto matches = match(str, tfirstRegex);
  auto triples = matches.map!(m => tuple(2.0 * m.captures[2].to!double(),
                                         m.captures[1].to!size_t(),
                                         m.captures[3].to!size_t()))();
  auto min_triple = triples.minCount!"a[0] < b[0]"()[0];
  if(min_triple[1] > min_triple[2])
    swap(min_triple[1], min_triple[2]);
  return min_triple[0];
}
