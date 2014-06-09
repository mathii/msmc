/* Copyright (c) 2012,2013 Genome Research Ltd.
 *
 * Author: Stephan Schiffels <stephan.schiffels@sanger.ac.uk>
 *
 * This file is part of psmc.
 * psmc is free software: you can redistribute it and/or modify it under
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
 
module model.propagation_core_naiveImpl;
import std.stdio;
import std.algorithm;
import std.conv;
import std.string;
import std.exception;
import model.data;
import model.gsl_matrix_vector;
import model.propagation_core;
import model.psmc_model;
import model.stateVec;
import model.stateVecAllocator;

class PropagationCoreNaive {
  
  gsl_matrix*[][] forwardPropagators, backwardPropagators;
  gsl_matrix*[] forwardPropagatorsMissing, backwardPropagatorsMissing;
  
  gsl_vector*[][] emissionProbs; // first index: i_ttot, second index: obs
  gsl_matrix* transitionMatrix;
  
  string[] allele_order;
  
  const PSMCmodel psmc;
  
  this(in PSMCmodel psmc, size_t maxDistance) {
    enforce(maxDistance > 0);
    this.psmc = psmc;
    this.allele_order = ["00", "01"];

    forwardPropagators = new gsl_matrix*[maxDistance];
    backwardPropagators = new gsl_matrix*[maxDistance];
    emissionProbs = new gsl_vector*[3];
      
    foreach(i; 0 .. 3) {
      emissionProbs[i] = gsl_vector_alloc(psmc.nrStates);
      foreach(a; 0 .. psmc.nrStates) {
        if(i == 0)
          gsl_vector_set(emissionProbs[i], a, 1.0); // missing data
        else
          gsl_vector_set(emissionProbs[i], a, psmc.emissionProb(allele_order[i - 1], a));
      }
    }
    
    transitionMatrix = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
    foreach(a; 0 .. psmc.nrStates) {
      foreach(b; 0 .. psmc.nrStates) {
        gsl_matrix_set(transitionMatrix, a, b, psmc.transitionProb(a, b));
      }
    }
      
    foreach(dist; 0 .. maxDistance) {
      forwardPropagators[dist] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
      backwardPropagators[dist] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
    }
    computeForwardPropagators(forwardPropagators, false, maxDistance);
    computeBackwardPropagators(backwardPropagators, false, maxDistance);

    forwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    backwardPropagatorsMissing = new gsl_matrix*[maxDistance];
    foreach(dist; 0 .. maxDistance) {
      forwardPropagatorsMissing[dist] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
      backwardPropagatorsMissing[dist] = gsl_matrix_alloc(psmc.nrStates, psmc.nrStates);
    }
    computeForwardPropagators(forwardPropagatorsMissing, true, maxDistance);
    computeBackwardPropagators(backwardPropagatorsMissing, true, maxDistance);
    
  }
  
  ~this() {
    foreach(dist; 0 .. forwardPropagators.length) {
      gsl_matrix_free(forwardPropagators[dist]);
      gsl_matrix_free(backwardPropagators[dist]);
    }
    foreach(dist; 0 .. forwardPropagatorsMissing.length) {
      gsl_matrix_free(forwardPropagatorsMissing[dist]);
      gsl_matrix_free(backwardPropagatorsMissing[dist]);
    }
    foreach(i; 0 .. emissionProbs.length)
      gsl_vector_free(emissionProbs[i]);
    gsl_matrix_free(transitionMatrix);
  }
  
  private void computeForwardPropagators(gsl_matrix*[] ret, bool missing_data, size_t maxDistance) const
  {
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : gsl_vector_get(emissionProbs[1], a);
      foreach(b; 0 .. psmc.nrStates) {
        auto val = gsl_matrix_get(transitionMatrix, a, b) * e;
        gsl_matrix_set(ret[0], a, b, val);
      }
    }

    foreach(distance; 1 .. maxDistance) {
      gsl_matrix_set_zero(ret[distance]);
      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
                     1.0, ret[0], ret[distance - 1], 0.0, ret[distance]);
      
    }
  }
  
  private void computeBackwardPropagators(gsl_matrix*[] ret, bool missing_data, size_t maxDistance) const
  {
    foreach(a; 0 .. psmc.nrStates) {
      double e = missing_data ? 1.0 : gsl_vector_get(emissionProbs[1], a);
      foreach(b; 0 .. psmc.nrStates) {
        auto val = psmc.transitionRate.transitionProbability(a, b) * e;
        gsl_matrix_set(ret[0], a, b, val);
      }
    }
    
    foreach(distance; 1 .. maxDistance) {
      gsl_matrix_set_zero(ret[distance]);

      gsl_blas_dgemm_checked(CBLAS_TRANSPOSE_t.CblasNoTrans, CBLAS_TRANSPOSE_t.CblasNoTrans,
        1.0, ret[distance - 1], ret[0], 0.0, ret[distance]);
    }
  }
  
  private double fullE(in SegSite_t segsite, size_t a) const {
    double ret = 0.0;
    foreach(o; segsite.obs) {
      ret += gsl_vector_get(emissionProbs[o], a);
    }
    ret /= cast(double)segsite.obs.length;
    return ret;
  }
  
  void propagateSingleForward(in State_t from, State_t to, 
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    to.setZero();
    
    foreach(a; 0 .. psmc.nrStates) {
      auto sum = 0.0;
      foreach(b; 0 .. psmc.nrStates) {
        sum += from.vec[b] * gsl_matrix_get(transitionMatrix, a, b);
      }
      to.vec[a] = fullE(to_segsite, a) * sum;
    }
  }
  
  void propagateSingleBackward(in State_t to, State_t from,
            in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos == from_segsite.pos + 1);
  }
  body {
    foreach(b; 0 .. psmc.nrStates) {
      auto sum = 0.0;
      foreach(a; 0 .. psmc.nrStates) {
        sum += to.vec[a] * fullE(to_segsite, a) * gsl_matrix_get(transitionMatrix, a, b);
      }
      from.vec[b] = sum;
    }
  }
  
  void propagateMultiForward(in State_t from, State_t to,
        in SegSite_t from_segsite, in SegSite_t to_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(a; 0 .. psmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = forwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(b; 0 .. psmc.nrStates) {
          sum += from.vec[b] * gsl_matrix_get(prop, a, b);
        }
        to.vec[a] = sum;
      }
      else {
        auto prop = forwardPropagators[dist - 1];
        auto sum = 0.0;
        foreach(b; 0 .. psmc.nrStates) {
          sum += from.vec[b] * gsl_matrix_get(prop, a, b);
        }
        to.vec[a] = sum;
      }
    }
  }
  
  void propagateMultiBackward(in State_t to, State_t from,
        in SegSite_t to_segsite, in SegSite_t from_segsite) const
  in {
    assert(to_segsite.pos > from_segsite.pos);
    assert(to_segsite.obs[0] < 2);
  }
  body {  
    auto dist = to_segsite.pos - from_segsite.pos;
    foreach(b; 0 .. psmc.nrStates) {
      if(to_segsite.obs[0] == 0) {
        auto prop = backwardPropagatorsMissing[dist - 1];
        auto sum = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          sum += to.vec[a] * gsl_matrix_get(prop, a, b);
        }
        from.vec[b] = sum;
      }
      else {
        auto prop = backwardPropagators[dist - 1];
        auto sum = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          sum += to.vec[a] * gsl_matrix_get(prop, a, b);
        }
        from.vec[b] = sum;
      }
    }

  }
  
  const(PSMCmodel) getPSMC() const {
    return psmc;
  }
  
  @property size_t forwardStateSize() const {
    return psmc.nrStates;
  }

  @property size_t backwardStateSize() const {
    return psmc.nrStates;
  }
  
  State_t newForwardState() const {
    return new State_t(psmc.nrStates, 0, 0);
  }

  State_t newBackwardState() const {
    return new State_t(psmc.nrStates, 0, 0);
  }

  State_t newForwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, 0, 0, stateAllocator);
  }

  State_t newBackwardState(StateVecAllocator stateAllocator) const {
    return new State_t(psmc.nrStates, 0, 0, stateAllocator);
  }
  
  void initialState(State_t s) const {
    foreach(a; 0 .. psmc.nrStates) {
      auto val = psmc.equilibriumProb(a);
      s.vec[a] = val;
    }
  }
  
  void setState(State_t s, double x, in SegSite_t segsite) const {
    foreach(a; 0 .. psmc.nrStates)
      s.vec[a] = x;
  }
  
  void getTransitionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][] eMat) const
  {
    foreach(a; 0 .. psmc.nrStates)
      foreach(b; 0 .. psmc.nrStates)
        eMat[a][b] = f.vec[b] * gsl_matrix_get(transitionMatrix, a, b) * b.vec[a] * fullE(to_segsite, a);
  }

  void getEmissionExpectation(State_t f, State_t b, in SegSite_t to_segsite, double[][] eMat) const
  {
    foreach(i; 0 .. 2) {
      if(to_segsite.obs.canFind(i + 1)) {
        auto norm = 0.0;
        foreach(a; 0 .. psmc.nrStates) {
          auto n = f.vec[a] * b.vec[a] / cast(double)(to_segsite.obs.length);
          eMat[i][a] = n;
          norm += n;
        }
        foreach(a; 0 .. psmc.nrStates)
          eMat[i][a] /= norm;
      }
      else
        eMat[i][] = 0.0;
    }
  }
  
  @property size_t maxDistance() const {
    return cast(size_t)forwardPropagators.length;
  }
  
}

