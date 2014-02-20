/*    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    Copyright (C) <2014>  <Garnet J. Vaz>
*/

#include "qtree.h"

#include<algorithm>
#include<queue>
#include<vector>
#include<functional>
#include<cmath>

using namespace std;

typedef priority_queue< double, vector<double> > minHeap;
typedef priority_queue< double, vector<double>, std::greater<double> > maxHeap;

void getLeftQad(const double *ys, double *qd,
		const double tau, uint ylen,
		double& quant);
void getRightQad(const double *ys, double *qd,
		 const double tau, uint ylen);

void getQad(const arma::vec& xs, const arma::vec& yvals, arma::vec& qad, const double tau, uint minSize, double& cut, double& minQad, double& quant) {
  qad.zeros();
  arma::uvec xsortix = arma::sort_index(xs);
  unsigned int minInd = minSize;
  unsigned int n = qad.n_elem - minSize;
  arma::vec xis = xs.elem(xsortix);
  arma::vec ys = yvals.elem(xsortix);

  double *y = new double[ys.n_elem];
  double *ypt = y;
  for(unsigned int i=0; i<ys.n_elem; ++i) y[i] = ys.at(i);
  double *qd = new double[qad.n_elem];
  double *qpt = qd;
  getLeftQad(ypt, qpt, tau, ys.n_elem, quant);
  ypt = y;
  qpt = qd;
  getRightQad(ypt, qpt, tau, ys.n_elem);

  double min = *(qd+minSize);
  for(unsigned int i=minSize+1;  i<n; ++i) {
    if((*(qd+i) < min) && (xis(i-1) < xis(i))) {
      min = *(qd+i);
      minInd = i;
    }
  }

  cut = 0.5*(xis(minInd-1) + xis(minInd));
  minQad = min;

  for(unsigned int i = 0; i<qad.n_elem; ++i) qad(i) = *(qd+i);
  delete[] y;
  delete[] qd;

}

void getLeftQad(const double *ys, double *qd,
		const double tau, uint ylen,
		double& quant) {
  minHeap low;
  maxHeap high;
  double first = *ys;
  low.push(first);
  double qp, qr; 	// Quantile based value here. Name might be misleading.
  qp = first;
  qr = 0.0;
  *qd++ = 0.0;
  *qd++ = 0.0;
  int nl, nr, nn;		// nn = total points, nl = points in left
  int nlold, nrold;
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;

  double k, i, j, l;
  int test;			// Used to test a shift from one heap to another.

  ys++;
  for(unsigned int ii = 1; ii < ylen; ++ii, ys++, qd++) {
    nlold = nl;
    nrold = nr;
    k = i = j = 0.;
    if(*ys <= low.top()) {
      low.push(*ys);
      k = 1.;
      nl++;
    } else {
      high.push(*ys);
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
    *qd = *(qd-1) + (qp - qr) * ((tau-1.)*nlold + tau*nrold) \
      + i*j*(l-qr) + k*(tau-1.)*(*ys-qr) + (1.-k)*tau*(*ys-qr);

    // Do updates
    qp = qr;
  }
  quant = qr;			// Value of quantile with all points.

}

void getRightQad(const double *ys, double *qd,
		 const double tau, const uint ylen) {
  minHeap low;
  maxHeap high;
  const double *ypt = ys+ylen-1;
  double *qpt = qd+ylen+1;
  double last = *ypt;
  low.push(last);
  double qp, qr;
  qp = last;
  qr = 0.0;
  int nlold, nrold;
  int nl, nr, nn;
  nn = 1;
  nlold = nl = 1;
  nrold = nr = 0;
  double k, i, j, l;
  int test;			// Used to test a shift from one heap to another.
  double qadp, qadr;
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
			 arma::vec& y,
			 arma::mat& Xmat,
			 double sroot,
			 double mindev,
			 double minsize,
			 double tau)
{
  ourVector output;
  output.empty = true;
  double stemp = sroot;
  unsigned int indx = 0;
  unsigned int i;
  arma::vec Xmatcoli;
  arma::vec qad(y.n_elem+1);
  double cut, minQad, quant;
  arma::uvec cutLeft, cutRight;
  double v;

  if(y.n_elem == 2) {
    output.li = Xmat.at(0,0) < Xmat.at(0,1) ? 0 : 1;
    output.ri = Xmat.at(0,0) < Xmat.at(0,1) ? 1 : 0;
    output.i = 0;
    output.val = 0.5*(Xmat.at(0,0) + Xmat(0,1));
    output.empty = false;
    output.quantile = y(0) + (y(1) - y(0))*(tau);
    output.sold = 0.0;
    // cout << "used weird case" << endl;
    return output;
  }

  for (i=0; i<Xmat.n_cols; i++)
  {

    Xmatcoli = Xmat.col(i);
    cut = quant = minQad = 0.0;

    getQad(Xmatcoli, y, qad, tau, minsize, cut, minQad, quant);

    if(minQad < stemp) {
      stemp = minQad;
      v = cut;
      indx = i;
      cutRight = arma::find(Xmatcoli > cut);
      cutLeft = arma::find(Xmatcoli <= cut);
    }

  }

  if (((qad(0)-stemp) > (mindev*sroot)) &&
      (cutLeft.n_elem >= minsize) &&
      (cutRight.n_elem >= minsize))
  {
    output.li = indices.elem(cutLeft);
    output.ri = indices.elem(cutRight);
    output.i = indx;
    output.val = v;
    output.empty = false;
    output.quantile = quant;
    output.sold = qad[0];
    return output;
  }
  else {
    output.quantile = quant;
    output.sold = qad(0);
    return output;
  }
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
