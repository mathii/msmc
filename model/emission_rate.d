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

double emissionProb(string alleles, double mu, double t, size_t ind1, size_t ind2, double tTot) {
  auto M = alleles.length;
  if(M == 2) {
    if(alleles[0] == alleles[1])
      return exp(-2.0 mu * t);
    else
      return (1.0 - exp(-2.0 * mu * t));
  }
  
  if(tTot < M * t) {
    tTot = M * t;
  }
  
  auto p_0 = emissionProb(alleles, '0', mu, t, ind1, ind2);
  auto p_1 = emissionProb(alleles, '1', mu, t, ind1, ind2);
  return p_0 + p_1;
}

double emissionProb(string alleles, char ancestralAllele, double mu, double t, size_t ind1, size_t ind2)
in {
  assert(t > 0.0);
  assert(t < double.infinity);
}
out(res) {
  assert(res >= 0.0);
}
body {
  auto M = alleles.length;
  auto count_ancestral = count(alleles, ancestralAllele);
  auto count_derived = alleles.length - count_ancestral;
  auto wattersonMm1 = iota(1, M - 1).map!"1.0/a"().reduce!"a+b"();

  auto tTotUpper = 2.0 * wattersonMm1;
  auto tTot = tTotUpper + M * t;
  
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
      return (1.0 - exp(-mu * tTot)) * (t + 2.0 / (M - 1.0)) / tTot;
  }
  else {
    if(alleles[ind1] != alleles[ind2])
      return 0.0;
    else
      return (1.0 - exp(-mu * tTot)) * (tTotUpper - 2.0 / (M - 1.0)) * 1.0 / (2.0 ^^ (M - 1) - (M - 2.0)) * 1.0 / tTot;
  }
}

