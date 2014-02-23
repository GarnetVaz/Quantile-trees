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

void getLeftQad(double *y, double *qad,	double tau, int ylen, double& quant);
void getRightQad(double *y, double *qad, double tau, int ylen);

void getQad(double *x,
	    double *y,
	    double *qd,
	    double tau,
	    int minCut,
	    int ylen,
	    double &cut,
	    double &minQad,
	    double &quant,
	    int &nleft) {
  // Compute left and right qad's
  getLeftQad(y, qd, tau, ylen, quant);
  getRightQad(y, qd, tau, ylen);
  double min = *(qd+minCut);
  unsigned int minInd = minCut;
  for(int i=minCut+1; i <= (ylen-minCut); ++i) {
    if((qd[i] < min) && (x[i-1] < x[i])) {
      min = qd[i];
      minInd = i;
    }
  }
  cut = 0.5*(x[minInd-1] + x[minInd]);
  minQad = min;
  nleft = minInd;
}

void getLeftQad(double *y,
		double *qad,
		double tau,
		int ylen,
		double& quant) {
  minHeap low;
  maxHeap high;
  double *ys, *qd;
  ys = y;
  qd = qad;
  double first = ys[0];
  low.push(first);
  double qp, qr; 	// Quantile based value here. Name might be misleading.
  qp = first;
  qr = 0.0;
  qd[0] = qd[1] = 0.0;
  int nl, nr, nn;		// nn = total points, nl = points in left
  int nlold, nrold;
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;

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
    qr = low.top() + (high.top()-low.top())*(tau*((double)nn -1.) -((double)test-1.));

    // Get new QAD from old one.
    qd[ii+1] = qd[ii] + (qp - qr) * ((tau-1.)*nlold + tau*nrold) \
      + i*j*(l-qr) + k*(tau-1.)*(ys[ii]-qr) + (1.-k)*tau*(ys[ii]-qr);

    // Do updates
    qp = qr;
  }
  quant = qr;			// Value of quantile with all points.

}

void getRightQad(double *y,
		 double *qad,
		 double tau,
		 int ylen) {
  minHeap low;
  maxHeap high;
  double *ys, *qd, *ypt, *qpt;
  double k, i, j, l, last;
  int test;			// Used to test a shift from one heap to another.
  double qadp, qadr, qp, qr;
  int nl, nr, nn, nlold, nrold;
  ys = y;
  qd = qad;
  ypt = ys+ylen-1;
  qpt = qd+ylen+1;
  last = ypt[0];
  qp = last;
  qr = 0.0;
  low.push(last);
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;
  qadp = qadr = 0.0;
  qpt -= 3;
  ypt--;
  for(int ii=ylen-1; ii > 0; --ii, --ypt, --qpt) {
    nlold = nl;
    nrold = nr;
    qadp = qadr;
    k = i = j = 0.;
    if(*ypt <= low.top()) {
      low.push(*ypt);
      k = 1.;
      nl++;
    } else {
      high.push(*ypt);
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
    qr = low.top() + (high.top()-low.top())*(tau*((double)nn -1.) -((double)test-1.));

    // Get new QAD from old one.
    qadr = qadp + (qp - qr) * ((tau-1.)*nlold + tau*nrold) \
      + i*j*(l-qr) + k*(tau-1.)*(*ypt-qr) + (1.-k)*tau*(*ypt-qr);
    *qpt += qadr;

    // Do updates
    qp = qr;
  }
}

nodeStruct splitNode(vector< unsigned int>& indices,
		     NumericVector& yvals,
		     NumericMatrix& Xmat,
		     double sroot,
		     double mindev,
		     double minCut,
		     double tau)
{
  unsigned int nNode = indices.size();
  unsigned int nPredictors = Xmat.ncol();
  nodeStruct output;
  output.empty = true;
  double stemp = sroot;
  unsigned int indx, ui, uj;
  indx = 0;
  // int i, j;
  double cut, minQad, quant;
  vector<unsigned int> cutLeft, cutRight;
  int nLeft = 0;

  // If there are only two points we might need this.
  // if(y.n_elem == 2) {
  //   output.li = Xmat.at(0,0) < Xmat.at(0,1) ? 0 : 1;
  //   output.ri = Xmat.at(0,0) < Xmat.at(0,1) ? 1 : 0;
  //   output.i = 0;
  //   output.val = 0.5*(Xmat.at(0,0) + Xmat(0,1));
  //   output.empty = false;
  //   output.quantile = y(0) + (y(1) - y(0))*(tau);
  //   output.sold = 0.0;
  //   // cout << "used weird case" << endl;
  //   return output;
  // }

  double *x = new double[nNode];
  double *xCopy = new double[nNode]; // Need to find cuts.
  double *y = new double[nNode];
  double *ySort = new double[nNode];
  double *qd = new double[nNode+1];
  int *index = new int[nNode];


  // Copy y values over to array.
  for(ui=0; ui<nNode; ++ui) y[ui] = yvals(indices[ui]);

  // Start loop over every column
  for (ui=0; ui<nPredictors; ++ui)  {
    // Initialize values for i'th column
    for(uj=0; uj<nNode; ++uj) {
      x[uj] = Xmat(uj,ui);
      xCopy[uj] = Xmat(uj,ui);
      index[uj] = uj;
      qd[uj] = 0.0;
    }
    cut = quant = minQad = 0.0;
    nLeft = 0;
    R_qsort_I(x, index, 0, nNode);
    for(unsigned int um=0; um<nNode; ++um) ySort[um] = y[index[um]];

    getQad(x, ySort, qd, tau, minCut, nNode, cut, minQad, quant, nLeft);

    if(minQad < stemp) {
      stemp = minQad;
      indx = ui;
      cutLeft.resize(nLeft);
      cutRight.resize(nNode-nLeft);
      unsigned int jLeft, jRight;
      jLeft = jRight = 0;
      for(unsigned int ii=0; ii < nNode; ++ii) {
	if(xCopy[ii] <= cut) cutLeft[jLeft++] = ii;
	else cutRight[jRight++] = ii;
      }
    }
  }

  // qd[0] represents sold here.
  if ((qd[0]-stemp) > (mindev*sroot))
  {
    output.li.resize(nLeft);
    for(ui=0; ui<nLeft; ++ui) output.li[ui] = indices[cutLeft[ui]];
    output.ri.resize(nNode-nLeft);
    for(ui=0; ui<nNode-nLeft-1; ++ui) output.ri[ui] = indices[cutRight[ui]];
    output.i = indx;
    output.val = cut;
    output.empty = false;
    output.quantile = quant;
    output.sold = qd[0];
  } else {
    output.quantile = quant;
    output.sold = qd[0];
  }
  delete [] x;
  delete [] xCopy;
  delete [] y;
  delete [] ySort;
  delete [] qd;
  delete [] index;
  return output;
}

void getQuantileAndQAD(const NumericVector& ys, double& quant, double& qad, const double tau) {
  int size = ys.size();
  if(size == 1) {
    qad = 0.0;
    quant = ys(0);
    return;
  }
  int low = ceil((size-1)*tau);
  std::partial_sort(ys.begin(),ys.begin()+low+1,ys.end());
  quant = ys(low-1) + (ys(low)-ys(low-1))*(tau*(size -1.) -(low-1.));
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
  vector< vector<int> > leaflist;
  vector<int> nodeID, nodeIDList;
  vector<int> n, valguide;
  string s1;
  unsigned int ui;
  int i;
  double quant;
  nodeStruct splitOut;

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
    if (nNode <= (unsigned int) minSize)
    {
      NumericVector yvec(nNode);
      for(ui=0; ui<nNode; ++ui) yvec(ui) = resp(indices[ui]);
      getQuantileAndQAD(yvec, quant, sold, tau);
      yval.push_back(quant);
      nodeID.push_back(nodeIDList[0]);
      n.push_back(nNode);
      dev.push_back(sold);
      for(ui=0; ui<nNode; ++ui) yhat[ui] = quant;
      leaflist.push_back(vector<int>(indices.begin(),indices.end()));

      // delete this node from the list of active nodes
      activelist.erase(activelist.begin());
      nodeIDList.erase(nodeIDList.begin());
      var.push_back(xlevels[0]);
      val.push_back(0);
      valguide.push_back(0); // 0 is for empty
    }
    else
    {
      splitOut = splitNode(indices, resp, pred, sroot, minDev, minCut, tau);
      for(ui=0; ui<nNode; ++ui) yhat[ui] = splitOut.quantile;
      dev.push_back(splitOut.sold);
      if (! splitOut.empty)
      {
        activelist.erase(activelist.begin());
	nodeID.push_back(nodeIDList[0]);
	unsigned int nodeNum = nodeIDList[0];
        nodeIDList.erase(nodeIDList.begin());
        activelist.insert(activelist.begin(), 1, splitOut.ri);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum+1);
        activelist.insert(activelist.begin(), 1, splitOut.li);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum);
        var.push_back(xlevels[(splitOut.i)+1]);
        val.push_back(splitOut.val);
        valguide.push_back(1);    // 1 is for full
	n.push_back(nNode);
	yval.push_back(splitOut.quantile);
      }
      else
      {
	yval.push_back(splitOut.quantile);
	nodeID.push_back(nodeIDList[0]);
	n.push_back(nNode);
	leaflist.push_back(vector<int>(indices.begin(),indices.end()));
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
