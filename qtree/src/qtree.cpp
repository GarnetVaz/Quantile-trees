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

using namespace std;

typedef priority_queue< double, vector<double> > minHeap;
typedef priority_queue< double, vector<double>, std::greater<double> > maxHeap;

void getLeftQad(double *y, double *qad,	double tau, int ylen, double& quant);
void getRightQad(double *y, double *qad, double tau, int ylen);

void getQad(double *x,
	    double *y,
	    double *qd,
	    double tau,
	    int minSize,
	    int ylen,
	    double cut,
	    double minQad,
	    double quant,
	    int nleft) {
  // Compute left and right qad's
  getLeftQad(y, qd, tau, ylen, quant);
  getRightQad(y, qd, tau, ylen);
  double min = *(qd+minSize);
  unsigned int minInd = minSize;
  for(int i=minSize+1; i <= (ylen-minSize); ++i) {
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

ourVector myfun(arma::uvec& indices,
			 arma::vec& yvals,
			 arma::mat& Xmat,
			 double sroot,
			 double mindev,
			 double minSize,
			 double tau)
{
  ourVector output;
  output.empty = true;
  double stemp = sroot;
  unsigned int indx = 0;
  uint i, j;
  double cut, minQad, quant;
  arma::uvec cutLeft, cutRight;
  uint nleft = 0;
  double v;

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
  uint ylen = yvals.n_elem;
  double *x = new double[ylen];
  double *xCopy = new double[ylen]; // Need to find cuts.
  double *y = new double[ylen];
  double *ySort = new double[ylen];
  double *qd = new double[ylen+1];
  int *index = new int[ylen];


  // Copy y values over to array.
  for(i=0; i<ylen; ++i) y[i] = yvals(i);

  // Start loop over every column
  for (i=0; i<Xmat.n_cols; i++)  {
    // Initialize values for i'th column
    for(j=0; j<ylen; ++j) {
      x[j] = Xmat(j,i);
      xCopy[j] = Xmat(j,i);
      index[j] = j;
      qd[j] = 0.0;
    }
    cut = quant = minQad = 0.0;
    nleft = 0;
    R_qsort_I(x, index, 0, ylen);
    for(int m=0; m<ylen; ++m) ySort[m] = y[index[m]];

    getQad(x, y, qd, tau, minSize, (int) ylen, cut, minQad, quant, nleft);
    if((minQad < stemp) && (nleft > minSize)) {
      stemp = minQad;
      v = cut;
      indx = i;
      cutLeft.resize(nleft);
      cutRight.resize(ylen-nleft);
      uint jLeft, jRight;
      jLeft = jRight = 0;
      for(uint ii=0; ii < nleft; ++ii) {
	if(*(xCopy+ii) <= v) cutLeft(jLeft++) = ii;
	else cutRight(jRight++) = ii;
      }
    }
  }
  // *qd represents sold here.
  if (((*(qd)-stemp) > (mindev*sroot)) &&
      (cutLeft.n_elem >= minSize) && // Gauranteed! change this!
      (cutRight.n_elem >= minSize))
  {
    output.li = indices.elem(cutLeft);
    output.ri = indices.elem(cutRight);
    output.i = indx;
    output.val = v;
    output.empty = false;
    output.quantile = quant;
    output.sold = *qd;
  } else {
    output.quantile = quant;
    output.sold = *qd;
  }
  delete [] x;
  delete [] xCopy;
  delete [] y;
  delete [] ySort;
  delete [] qd;
  delete [] index;
  return output;
}

void getQuantileAndQAD(const arma::vec& ys, double& quant, double& qad, const double tau) {
  if(ys.n_elem == 1) {
    qad = 0.0;
    quant = ys(0);
    return;
  }
  arma::vec sorty = arma::sort(ys);
  unsigned int size = sorty.n_elem;
  qad = 0.0;
  int low = ceil((size-1)*tau);
  quant = sorty(low-1) + (sorty(low)-sorty(low-1))*(tau*((double)size -1.) -(low-1.));
  qad = 0.0;

  for(arma::vec::const_iterator it=ys.begin(); it!= ys.end(); ++it) {
    if(*it < quant) qad += (tau-1.)*(*it-quant);
    else qad += tau*(*it-quant);
  }

}

SEXP qtreeCPP(SEXP s_mypred,
                  SEXP s_myresp,
                  SEXP s_mindev,
                  SEXP s_mincut,
                  SEXP s_minsize,
                  SEXP s_mytau)
{
  Rcpp::NumericMatrix rs_mypred(s_mypred);
  arma::mat mypred(rs_mypred.begin(), rs_mypred.nrow(), rs_mypred.ncol(), false);

  Rcpp::NumericVector rs_myresp(s_myresp);
  arma::vec myresp(rs_myresp.begin(), rs_myresp.size(), false);

  double mindev = Rcpp::as<double>(s_mindev);
  double mincut = Rcpp::as<double>(s_mincut);
  double minsize = Rcpp::as<double>(s_minsize);
  double mytau = Rcpp::as<double>(s_mytau);
  arma::vec yhat(myresp.n_elem); // Predicted value for each training point.
  arma::uvec indices;
  double sroot, sold;
  vector<string> xlevels, var;
  vector<double> val;
  vector<double> dev;
  vector<double>  yval;
  vector< arma::uvec > activelist;
  vector< vector<int> > leaflist;
  vector<int> nodeID, nodeIDList;
  vector<int> n, valguide;
  string s1;
  unsigned int i;
  ourVector splitout;

  double quant;
  getQuantileAndQAD(myresp, quant, sroot, mytau);
  yhat.fill(quant);
  xlevels.push_back("<leaf>");
  for (i=0; i<mypred.n_cols; i++)
  {
    s1 = "X";
    stringstream tempstream;
    tempstream << (i+1);
    s1 += tempstream.str();
    xlevels.push_back(s1);
  }

  // arma::uvec allcols = arma::linspace<arma::uvec>(0,mypred.n_cols);
  indices = arma::linspace<arma::uvec>(0, (mypred.n_rows-1), mypred.n_rows);
  activelist.push_back(indices);
  nodeIDList.push_back(1);

  while (activelist.size() > 0)
  {
    // indices of 1st node in active list
    indices = activelist[0];

    arma::vec yvec = myresp.elem(indices);

    // check if we do not need to split current partition at all
    if (indices.n_elem <= mincut)
    {
      getQuantileAndQAD(yvec, quant, sold, mytau);
      yval.push_back(quant);
      nodeID.push_back(nodeIDList[0]);
      n.push_back(indices.size());
      dev.push_back(sold);
      (yhat.elem(indices)).fill(quant);
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
      arma::mat Xmat = mypred.rows(indices);
      splitout = myfun(indices, yvec, Xmat, sroot, mindev, minsize, mytau);
      (yhat.elem(indices)).fill(splitout.quantile);
      dev.push_back(splitout.sold);
      if (! splitout.empty)
      {
        activelist.erase(activelist.begin());
	nodeID.push_back(nodeIDList[0]);
	unsigned int nodeNum = nodeIDList[0];
        nodeIDList.erase(nodeIDList.begin());
        activelist.insert(activelist.begin(), 1, splitout.ri);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum+1);
        activelist.insert(activelist.begin(), 1, splitout.li);
        nodeIDList.insert(nodeIDList.begin(), 1, 2*nodeNum);
        var.push_back(xlevels[(splitout.i)+1]);
        val.push_back(splitout.val);
        valguide.push_back(1);    // 1 is for full
	n.push_back(indices.size());
	yval.push_back(splitout.quantile);
      }
      else
      {
	yval.push_back(splitout.quantile);
	nodeID.push_back(nodeIDList[0]);
	n.push_back(indices.size());
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
