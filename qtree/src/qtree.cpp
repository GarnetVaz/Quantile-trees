/*******************************************************************
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Copyright (C) <2014>  <Garnet J. Vaz>
*******************************************************************/

/*
  R based quantile tree code.
*/
#include "qtree.h"

#include<Rmath.h>
#include<R.h>

#include<algorithm>
#include<queue>
#include<vector>
#include<functional>
#include<cmath>
#include<iostream>

typedef priority_queue< double, vector<double> > minHeap;
typedef priority_queue< double, vector<double>, std::greater<double> > maxHeap;

void getLeftQad(const vector<double>& y, vector<double>& qad,
		double tau, int ylen, double& quant);
void getRightQad(const vector<double>& y, vector<double>& qad,
		 double tau, int ylen);

void getQad(const vector<double>& x,
	    const vector<double>& y,
	    vector<double>& qad,
	    double tau,
	    int minCut,
	    int ylen,
	    double &cut,
	    double &minQad,
	    double &quant) {
  // Compute left and right qad's
  getLeftQad(y, qad, tau, ylen, quant);
  getRightQad(y, qad, tau, ylen);
  double min = qad[minCut];
  unsigned int minInd = minCut;
  for(int i=minCut+1; i <= (ylen-minCut); ++i) {
    if((qad[i] < min) && (x[i-1] < x[i])) {
      min = qad[i];
      minInd = i;
    }
  }
  cut = 0.5*(x[minInd-1] + x[minInd]);
  minQad = min;
}

void getLeftQad(const vector<double>& ys,
		vector<double>& qad,
		double tau,
		int ylen,
		double& quant) {
  minHeap low;
  maxHeap high;
  low.push(ys[0]);
  int nl, nr, nn;		// nn = total points, nl = points in left
  int nlold, nrold;
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;
  qad[0] = qad[1] = 0.0;
  double qp, qr;
  qr = 0.0;
  qp = ys[0];
  double k, i, j, l;
  int test;			// Used to test a shift from one heap to another.

  // qd is ahead of y by 1 value.
  for(int ii = 1; ii < ylen; ++ii) {
    nlold = nl;
    nrold = nr;
    k = i = j = 0.;
    if(ys[ii] <= low.top()) {
      low.push(ys[ii]);
      k = 1.;
      nl++;
    } else {
      high.push(ys[ii]);
      nr++;
    }
    ++nn;

    // Balance the heaps to maintain the quantile
    test = ceil((nn-1)*tau);
    if(nl > test) {
      // Push from low to high.
      high.push(low.top());
      l = low.top();
      low.pop();
      i = 1.0;
      j = 1.0;
      nl--;
      nr++;
    } else if(nl < test) {
      // Push from high to low.
      low.push(high.top());
      l = high.top();
      high.pop();
      i = 1.0;
      j = -1.0;
      nr--;
      nl++;
    }

    // Compute the new median using linear interpolation
    qr = low.top() + (high.top()-low.top())*(tau*(nn -1.) -(test-1.));

    // Get new QAD from old one.
    qad[ii+1] = qad[ii] + (qp - qr) * ((tau-1.)*nlold + tau*nrold) \
      + i*j*(l-qr) + k*(tau-1.)*(ys[ii]-qr) + (1.-k)*tau*(ys[ii]-qr);

    // Do updates
    qp = qr;
  }
  quant = qr;			// Value of quantile with all points.
}

void getRightQad(const vector<double>& y,
		 vector<double>& qad,
		 double tau, int ylen) {
  minHeap low;
  maxHeap high;
  double k, i, j, l;
  int test;			// Used to test a shift from one heap to another.
  double qadp, qadr, qp, qr;
  int nl, nr, nn, nlold, nrold;
  low.push(y[ylen-1]);
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;
  qadp = qadr = 0.0;
  qp = y[ylen-1];
  for(int ii=ylen-2; ii >= 0; --ii) {
    nlold = nl;
    nrold = nr;
    qadp = qadr;
    k = i = j = 0.;
    if(y[ii] <= low.top()) {
      low.push(y[ii]);
      k = 1.;
      nl++;
    } else {
      high.push(y[ii]);
      nr++;
    }
    ++nn;

    test = ceil((nn-1)*tau);
    // Balance heaps
    if(nl > test) {
      // Push from low to high.
      high.push(low.top());
      l = low.top();
      low.pop();
      i = 1.0;
      j = 1.0;
      nl--;
      nr++;
    } else if(nl < test) {
      // Push from high to low.
      low.push(high.top());
      l = high.top();
      high.pop();
      i = 1.0;
      j = -1.0;
      nr--;
      nl++;
    }

    // Compute the new median using linear interpolation
    qr = low.top() + (high.top()-low.top())*(tau*(nn -1.) -(test-1.));

    // Get new QAD from old one.
    qadr = qadp + (qp - qr) * ((tau-1.)*nlold + tau*nrold) \
      + i*j*(l-qr) + k*(tau-1.)*(y[ii]-qr) + (1.-k)*tau*(y[ii]-qr);
    qad[ii] += qadr;

    // Do updates
    qp = qr;
  }
}

void splitNode(vector< unsigned int>& indices,
	       NumericVector yvals,
	       NumericMatrix Xmat,
	       vector< unsigned int>& left,
	       vector< unsigned int>& right,
	       double sroot,
	       double mindev,
	       double minCut,
	       double &val,
	       double &quantile,
	       double &sold,
	       bool &empty,
	       unsigned int &indx,
	       double tau)
{
  unsigned int nNode = indices.size();
  unsigned int nPredictors = Xmat.ncol();
  double stemp = sroot;
  unsigned int ui, uj;
  double cut, bestCut, minQad, quant, bestQuant;
  unsigned int nLeft = 0;
  bestCut = cut = bestQuant = 0.0;

  vector<double> x(nNode), xCopy(nNode), y(nNode), ySort(nNode);
  vector<unsigned int> cutLeft, cutRight;
  vector<unsigned int> leftTemp, rightTemp;
  vector<double> qd(nNode+1);
  vector<int> index(nNode);

  // Copy yvals at this node over to vector y.
  for(ui=0; ui<nNode; ++ui) y[ui] = yvals(indices[ui]);

  // Start loop over every column
  for (ui=0; ui<nPredictors; ++ui)  {
    // Initialize values for i'th column
    for(uj=0; uj<nNode; ++uj) {
      x[uj] = xCopy[uj] = Xmat(indices[uj],ui);
      index[uj] = (int) uj;
      qd[uj] = 0.0;
    }
    rsort_with_index(&x[0],&index[0],nNode);
    // R_qsort_I(&x[0], &index[0], 1, nNode);
    bool diff = false;		// used to test if node is pure
    for(unsigned int um=0; um<nNode; ++um) {
      ySort[um] = y[index[um]];
      x[um] = xCopy[index[um]];
    }
    for(unsigned int um=1; um<nNode; ++um) {
      if(y[um-1] != y[um]) {
	diff = true;
	break;
      }
    }

    cut = minQad = 0.0;
    quant = x[0];
    // If node is pure
    if(!diff) {
      quant = y[0];
      qd[0] = 0.0;
      break;
    } else {
      getQad(x, ySort, qd, tau, minCut, (int) nNode, cut, minQad, quant);
    }

    if(minQad < stemp) {
      stemp = minQad;
      indx = ui;
      unsigned int jLeft, jRight;
      bestCut = cut;
      // bestQuant = quant;
      jLeft = jRight = 0;
      cutLeft.clear();
      cutRight.clear();
      for(unsigned int ii=0; ii < nNode; ++ii) {
	if(xCopy[ii] <= bestCut)  {
	  cutLeft.push_back(ii);
	  jLeft++;
	} else {
	  cutRight.push_back(ii);
	}
      }
      nLeft = jLeft;
    }
  }
  quantile = quant;
  sold = qd[0];

  // qd[0] represents sold here.
  if ((sold-stemp) > (mindev*sroot))
  {
    left.resize(nLeft);
    right.resize(nNode-nLeft);
    for(ui=0; ui<nLeft; ++ui) left[ui] = indices[cutLeft[ui]];
    for(ui=0; ui<nNode-nLeft; ++ui) right[ui] = indices[cutRight[ui]];
    // leftTemp.resize(nLeft);
    // rightTemp.resize(nNode-nLeft);
    // for(ui=0; ui<nLeft; ++ui) leftTemp[ui] = indices[cutLeft[ui]];
    // for(ui=0; ui<nNode-nLeft; ++ui) rightTemp[ui] = indices[cutRight[ui]];
    val = bestCut;
    empty = false;
  }
  // if(!empty) {
  //   left = leftTemp;
  //   right = rightTemp;
  // }
}

void getQuantileAndQAD(const NumericVector& ys, double& quant, double& qad, const double tau) {
  int size = ys.size();
  if(size == 1) {
    qad = 0.0;
    quant = ys(0);
    return;
  }
  vector<double> y(size);
  for(int ui=0; ui < size; ++ui) y[ui] = ys[ui];
  int low = ceil((size-1)*tau);
  std::partial_sort(y.begin(),y.begin()+low+1,y.end());
  quant = y[low-1] + (y[low]-y[low-1])*(tau*(size -1.) -(low-1.));
  qad = 0.0;
  for(NumericVector::iterator it=ys.begin(); it!= ys.end(); ++it) {
    if(*it < quant) qad += (tau-1.)*(*it-quant);
    else qad += tau*(*it-quant);
  }
}

// [[Rcpp::export]]
List qtreeCPP(NumericMatrix pred,
	      NumericVector resp,
	      double minDev,
	      int minCut,
	      int minSize,
	      double tau)
{
  unsigned int nSamples = pred.nrow();
  unsigned int nPredictors = pred.ncol();
  vector<double> yhat(nSamples);
  vector<unsigned int> indices;
  double sroot, sold;
  vector<string> xlevels, var;
  vector<double> val;
  vector<double> dev;
  vector<double>  yval;
  vector< vector<unsigned int> > activelist;
  vector< vector<unsigned int> > leaflist;
  vector<int> nodeID, nodeIDList;
  vector<int> n, valguide;
  string s1;
  unsigned int ui;
  double quant;

  getQuantileAndQAD(resp, quant, sroot, tau);
  yhat.assign(yhat.size(), quant);
  // Need to move this stuff to R.
  xlevels.push_back("<leaf>");
  for (ui=0; ui<nPredictors; ui++)
  {
    s1 = "X";
    stringstream tempstream;
    tempstream << (ui+1);
    s1 += tempstream.str();
    xlevels.push_back(s1);
  }

  indices.resize(nSamples);
  for(ui=0; ui<nSamples; ++ui) indices[ui] = ui;
  activelist.push_back(indices);
  nodeIDList.push_back(1);

  while (activelist.size() > 0)
  {
    // indices of 1st node in active list
    indices = activelist[0];
    unsigned int nNode = indices.size();

    // check if we do not need to split current partition at all
    if (nNode < (unsigned int) minSize)
    {
      NumericVector yvec(nNode);
      for(ui=0; ui<nNode; ++ui) yvec(ui) = resp(indices[ui]);
      getQuantileAndQAD(yvec, quant, sold, tau);
      yval.push_back(quant);
      nodeID.push_back(nodeIDList[0]);
      n.push_back(nNode);
      dev.push_back(sold);
      for(ui=0; ui<nNode; ++ui) yhat[ui] = quant;
      leaflist.push_back(indices);

      // delete this node from the list of active nodes
      activelist.erase(activelist.begin());
      nodeIDList.erase(nodeIDList.begin());
      var.push_back(xlevels[0]);
      val.push_back(0);
      valguide.push_back(0); // 0 is for empty
    }
    else
    {
      double value, quantile, sold;
      bool empty = true;
      unsigned int indx = 0;
      vector<unsigned int> left, right;
      value = quantile = sold = 0.0;
      splitNode(indices, resp, pred, left, right, sroot, minDev, minCut, value, quantile, sold,	empty, indx, tau);
      for(ui=0; ui<nNode; ++ui) yhat[ui] = quantile;
      dev.push_back(sold);
      if (!empty)
      {
	unsigned int nodeNum = nodeIDList[0];
        activelist.erase(activelist.begin());
	nodeID.push_back(nodeNum);
        nodeIDList.erase(nodeIDList.begin());
        activelist.insert(activelist.begin(), 1, right);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum+1);
        activelist.insert(activelist.begin(), 1, left);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum);
        var.push_back(xlevels[indx+1]);
        val.push_back(value);
        valguide.push_back(1);    // 1 is for full
	n.push_back(nNode);
	yval.push_back(quantile);
      }
      else
      {
	yval.push_back(quantile);
	nodeID.push_back(nodeIDList[0]);
	n.push_back(nNode);
	leaflist.push_back(indices);
	activelist.erase(activelist.begin());
	nodeIDList.erase(nodeIDList.begin());
	var.push_back(xlevels[0]);
	val.push_back(0);
	valguide.push_back(0); // 0 is for empty
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("val") = val,
    Rcpp::Named("valguide") = valguide,
    Rcpp::Named("var") = var,
    Rcpp::Named("xlevels") = xlevels,
    Rcpp::Named("n") = n,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("yval") = yval,
    Rcpp::Named("leaflist") = leaflist,
    Rcpp::Named("nodeID") = nodeID
  );

}
