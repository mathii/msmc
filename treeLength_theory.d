#!/usr/bin/env rdmd

import std.stdio;
import std.conv;
import std.math;
import std.range;
import std.algorithm;
import std.array;
import std.mathspecial;
import model.time_intervals;

TimeIntervals timeIntervals;
double[] lambda;

void main(string[] args) {
  auto t1 = args[1].to!double();
  auto t2 = args[2].to!double();
  auto r = args[3].to!double();
  auto M = args[4].to!size_t();
  auto bound = [0.0, t1, t2];
  lambda = [1.0, 1.0 / r, 1.0];
  timeIntervals = new TimeIntervals(bound ~ [double.infinity]);
  
  auto coalTimes = iota(2, M + 1).map!(k => getExpectedCoalescenceTime(0.0, M, k)).retro().array();
  writeln(coalTimes);
  auto totalBranchlength = 0.0;
  foreach(i; 0 .. coalTimes.length) 
    totalBranchlength += coalTimes[i] * (M - i);
  writeln(totalBranchlength);
}

double getExpectedCoalescenceTime(double u, size_t m, size_t k) {
  auto beta = timeIntervals.findIntervalForTime(u);
  
  auto firstTerm = 0.0;
  foreach(j; k .. m + 1) {
    auto val = c(m, k, j);
    val *= integralHelper(binomial(j, 2) * lambda[beta], u, timeIntervals.rightBoundary(beta));
    firstTerm += val;
  }
  
  auto secondTerm = 0.0;
  foreach(gamma; beta + 1 .. timeIntervals.nrIntervals) {
    foreach(j; k .. m + 1) {
      auto inner_sum = 0.0;
      if(gamma - 1 >= beta + 1)
        inner_sum = iota(beta + 1, gamma).map!(nu => lambda[nu] * timeIntervals.delta(nu)).reduce!"a+b"();
      auto val = exp(-binomial(j, 2) * (timeIntervals.rightBoundary(beta) - u) * lambda[beta]);
      val *= exp(-binomial(j, 2) * inner_sum);
      val *= c(m, k, j);
      val *= integralHelper(binomial(j, 2) * lambda[gamma], timeIntervals.leftBoundary(gamma), timeIntervals.rightBoundary(gamma));
      secondTerm += val;
    }
  }
  return firstTerm + secondTerm;
}

double integralHelper(double lambda, double lower, double upper) {
  return (1.0 - exp(-lambda * (upper - lower))) / lambda;
}

double binomial(size_t m, size_t k) {
  auto prod = 1.0;
  foreach(i; 1 .. k + 1) {
    prod *= cast(double)(m - (k - i)) / i;
  }
  return prod;
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
