/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <algorithm>
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"

#define OUTPI

// change to hadron EoS (e.g. Laine) to calculate v,T,mu at the surface
#define SWAP_EOS

using namespace std;

namespace output{  // a namespace containing all the output streams
  ofstream fkw, fkw_dim, fxvisc, fyvisc, fdiagvisc, fx,
     fy, fdiag, fz, faniz, f2d, ffreeze, maniz;
}

// returns the velocities in cartesian coordinates, fireball rest frame.
// Y=longitudinal rapidity of fluid
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb,
                            double &nq, double &ns, double &vx, double &vy,
                            double &Y) {
 double p, vz;
 c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
 double eta = getZ(c->getZ());
 //	Y = eta + TMath::ATanH(vz) ;
 Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
 vx = vx * cosh(Y - eta) / cosh(Y);
 vy = vy * cosh(Y - eta) / cosh(Y);
}

Fluid::Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny,
             int _nz, double _minx, double _maxx, double _miny, double _maxy,
             double _minz, double _maxz, double _dt, double eCrit) {
 eos = _eos;
 eosH = _eosH;
 trcoeff = _trcoeff;
 nx = _nx;
 ny = _ny;
 nz = _nz;
 minx = _minx;
 maxx = _maxx;
 miny = _miny;
 maxy = _maxy;
 minz = _minz;
 maxz = _maxz;
 dx = (maxx - minx) / (nx - 1);
 dy = (maxy - miny) / (ny - 1);
 dz = (maxz - minz) / (nz - 1);
 dt = _dt;

 cell = new Cell[nx * ny * nz];

 cell0 = new Cell;
 cell0->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0.,
                   0.);  // tau is not important here, since *0
 cell0->setAllM(0.);

 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    getCell(ix, iy, iz)->setPrev(X_, getCell(ix - 1, iy, iz));
    getCell(ix, iy, iz)->setNext(X_, getCell(ix + 1, iy, iz));
    getCell(ix, iy, iz)->setPrev(Y_, getCell(ix, iy - 1, iz));
    getCell(ix, iy, iz)->setNext(Y_, getCell(ix, iy + 1, iz));
    getCell(ix, iy, iz)->setPrev(Z_, getCell(ix, iy, iz - 1));
    getCell(ix, iy, iz)->setNext(Z_, getCell(ix, iy, iz + 1));
    getCell(ix, iy, iz)->setPos(ix, iy, iz);
   }

 output_nt = 0;
 output_nx = 0;
 output_ny = 0;

 //---- Cornelius init
 double arrayDx[4] = {dt, dx, dy, dz};
 cornelius = new Cornelius;
 cornelius->init(4, eCrit, arrayDx);
 ecrit = eCrit;
 vEff = 0.;
 EtotSurf = 0.0;
}

Fluid::~Fluid() {
 delete[] cell;
 delete cell0;
}

void Fluid::initOutput(const char *dir, double tau0) {
 char command[255];
 sprintf(command, "mkdir -p %s", dir);
 int return_mkdir = system(command);
 cout << "mkdir returns: " << return_mkdir << endl;
 string outx = dir;
 outx.append("/outx.dat");
 string outxvisc = dir;
 outxvisc.append("/outx.visc.dat");
 string outyvisc = dir;
 outyvisc.append("/outy.visc.dat");
 string outdiagvisc = dir;
 outdiagvisc.append("/diag.visc.dat");
 string outy = dir;
 outy.append("/outy.dat");
 string outdiag = dir;
 outdiag.append("/outdiag.dat");
 string outz = dir;
 outz.append("/outz.dat");
 string outaniz = dir;
 outaniz.append("/out.aniz.dat");
 string outmaniz = dir;
 outmaniz.append("/out.maniz.dat");
 string out2d = dir;
 out2d.append("/out2D.dat");
 string outfreeze = dir;
 outfreeze.append("/freezeout.dat");
 output::fx.open(outx.c_str());
 output::fy.open(outy.c_str());
 output::fz.open(outz.c_str());
 output::fdiag.open(outdiag.c_str());
 output::f2d.open(out2d.c_str());
 output::fxvisc.open(outxvisc.c_str());
 output::fyvisc.open(outyvisc.c_str());
 output::fdiagvisc.open(outdiagvisc.c_str());
 output::faniz.open(outaniz.c_str());
 output::maniz.open(outmaniz.c_str());
 output::ffreeze.open(outfreeze.c_str());
 //################################################################
 // important remark. for correct diagonal output, nx=ny must hold.
 //################################################################
 outputGnuplot(tau0);
 output::faniz << "#  tau  <<v_T>>  e_p  e'_p  (to compare with SongHeinz)\n";
// output::maniz <<  setw(5) << "tau" << setw(18) << "eps1" << setw(18) << "eps2" << setw(18) << "eps8" << setw(18) << "eps4" << setw(18) << "eps11" <<  setw(18) << "eps12" <<  setw(18) << "eps3" << setw(18) << "eps5" << setw(18) << "eps6" << setw(18) << "eps7" << setw(18)
//  << "eps9" << setw(18) << "eps10"  << setw(18) << "eps13"  << setw(18) << "eps14"  << setw(18) << "eps15"  << setw(18) << "eps16"  << setw(18) << "eps17"  << setw(18) << "eps18" << setw(18) << "eps19" << setw(18) << "eps20" << endl;
}

void Fluid::correctImagCells(void) {
 double Q[7];
 // Z
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   // left boundary
   getCell(ix, iy, 2)->getQ(Q);
   getCell(ix, iy, 1)->setQ(Q);
   getCell(ix, iy, 0)->setQ(Q);
   // right boundary
   getCell(ix, iy, nz - 3)->getQ(Q);
   getCell(ix, iy, nz - 2)->setQ(Q);
   getCell(ix, iy, nz - 1)->setQ(Q);
  }
 // Y
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   getCell(ix, 2, iz)->getQ(Q);
   getCell(ix, 1, iz)->setQ(Q);
   getCell(ix, 0, iz)->setQ(Q);
   // right boundary
   getCell(ix, ny - 3, iz)->getQ(Q);
   getCell(ix, ny - 2, iz)->setQ(Q);
   getCell(ix, ny - 1, iz)->setQ(Q);
  }
 // X
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   getCell(2, iy, iz)->getQ(Q);
   getCell(1, iy, iz)->setQ(Q);
   getCell(0, iy, iz)->setQ(Q);
   // right boundary
   getCell(nx - 3, iy, iz)->getQ(Q);
   getCell(nx - 2, iy, iz)->setQ(Q);
   getCell(nx - 1, iy, iz)->setQ(Q);
  }
}

void Fluid::correctImagCellsFull(void) {
 double Q[7], _pi[4][4], _Pi;
 // Z
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   // left boundary
   getCell(ix, iy, 2)->getQ(Q);
   getCell(ix, iy, 1)->setQ(Q);
   getCell(ix, iy, 0)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, iy, 2)->getpi(i, j);
   _Pi = getCell(ix, iy, 2)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, iy, 0)->setpi(i, j, _pi[i][j]);
     getCell(ix, iy, 1)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, iy, 0)->setPi(_Pi);
   getCell(ix, iy, 1)->setPi(_Pi);
   // right boundary
   getCell(ix, iy, nz - 3)->getQ(Q);
   getCell(ix, iy, nz - 2)->setQ(Q);
   getCell(ix, iy, nz - 1)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(ix, iy, nz - 3)->getpi(i, j);
   _Pi = getCell(ix, iy, nz - 3)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, iy, nz - 2)->setpi(i, j, _pi[i][j]);
     getCell(ix, iy, nz - 1)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, iy, nz - 2)->setPi(_Pi);
   getCell(ix, iy, nz - 1)->setPi(_Pi);
  }
 // Y
 for (int ix = 0; ix < nx; ix++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   getCell(ix, 2, iz)->getQ(Q);
   getCell(ix, 1, iz)->setQ(Q);
   getCell(ix, 0, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, 2, iz)->getpi(i, j);
   _Pi = getCell(ix, 2, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, 0, iz)->setpi(i, j, _pi[i][j]);
     getCell(ix, 1, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, 0, iz)->setPi(_Pi);
   getCell(ix, 1, iz)->setPi(_Pi);
   // right boundary
   getCell(ix, ny - 3, iz)->getQ(Q);
   getCell(ix, ny - 2, iz)->setQ(Q);
   getCell(ix, ny - 1, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(ix, ny - 3, iz)->getpi(i, j);
   _Pi = getCell(ix, ny - 3, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(ix, ny - 2, iz)->setpi(i, j, _pi[i][j]);
     getCell(ix, ny - 1, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(ix, ny - 2, iz)->setPi(_Pi);
   getCell(ix, ny - 1, iz)->setPi(_Pi);
  }
 // X
 for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
   // left boundary
   getCell(2, iy, iz)->getQ(Q);
   getCell(1, iy, iz)->setQ(Q);
   getCell(0, iy, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) _pi[i][j] = getCell(2, iy, iz)->getpi(i, j);
   _Pi = getCell(2, iy, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(0, iy, iz)->setpi(i, j, _pi[i][j]);
     getCell(1, iy, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(0, iy, iz)->setPi(_Pi);
   getCell(1, iy, iz)->setPi(_Pi);
   // right boundary
   getCell(nx - 3, iy, iz)->getQ(Q);
   getCell(nx - 2, iy, iz)->setQ(Q);
   getCell(nx - 1, iy, iz)->setQ(Q);
   for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
     _pi[i][j] = getCell(nx - 3, iy, iz)->getpi(i, j);
   _Pi = getCell(nx - 3, iy, iz)->getPi();

   for (int i = 0; i < 4; i++)
    for (int j = 0; j <= i; j++) {
     getCell(nx - 2, iy, iz)->setpi(i, j, _pi[i][j]);
     getCell(nx - 1, iy, iz)->setpi(i, j, _pi[i][j]);
    }
   getCell(nx - 2, iy, iz)->setPi(_Pi);
   getCell(nx - 1, iy, iz)->setPi(_Pi);
  }
}

void Fluid::updateM(double tau, double dt) {
 for (int ix = 0; ix < getNX(); ix++)
  for (int iy = 0; iy < getNY(); iy++)
   for (int iz = 0; iz < getNZ(); iz++) {
    Cell *c = getCell(ix, iy, iz);
    c->setDM(X_, 0.);
    c->setDM(Y_, 0.);
    c->setDM(Z_, 0.);
    if (getCell(ix, iy, iz)->getMaxM() < 1.) {
     if (getCell(ix + 1, iy, iz)->getM(X_) >= 1. ||
         getCell(ix - 1, iy, iz)->getM(X_) >= 1.)
      c->setDM(X_, dt / dx);
     if (getCell(ix, iy + 1, iz)->getM(Y_) >= 1. ||
         getCell(ix, iy - 1, iz)->getM(Y_) >= 1.)
      c->setDM(Y_, dt / dy);
     if (getCell(ix, iy, iz + 1)->getM(Z_) >= 1. ||
         getCell(ix, iy, iz - 1)->getM(Z_) >= 1.)
      c->setDM(Z_, dt / dz / tau);

     if (c->getDM(X_) == 0. && c->getDM(Y_) == 0.) {
      if (getCell(ix + 1, iy + 1, iz)->getMaxM() >= 1. ||
          getCell(ix + 1, iy - 1, iz)->getMaxM() >= 1. ||
          getCell(ix - 1, iy + 1, iz)->getMaxM() >= 1. ||
          getCell(ix - 1, iy - 1, iz)->getMaxM() >= 1.) {
       c->setDM(X_, 0.707 * dt / dx);
       c->setDM(Y_, 0.707 * dt / dy);
      }
     }
    }  // if
   }

 for (int ix = 0; ix < getNX(); ix++)
  for (int iy = 0; iy < getNY(); iy++)
   for (int iz = 0; iz < getNZ(); iz++) {
    Cell *c = getCell(ix, iy, iz);
    c->addM(X_, c->getDM(X_));
    c->addM(Y_, c->getDM(Y_));
    c->addM(Z_, c->getDM(Z_));
   }
}

void Fluid::outputGnuplot(double tau) {
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;

 // X direction
 for (int ix = 0; ix < nx; ix++) {
  double x = getX(ix);
  Cell *c = getCell(ix, ny / 2, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  output::fx << setw(14) << tau << setw(14) << x << setw(14) << vx << setw(14) << vy
        << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14) << mub;
  output::fx << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fx << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fx << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fx << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fx << endl;

 // Y direction
 for (int iy = 0; iy < ny; iy++) {
  double y = getY(iy);
  Cell *c = getCell(nx / 2, iy, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  output::fy << setw(14) << tau << setw(14) << y << setw(14) << vy << setw(14) << vx
        << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14) << mub;
  output::fy << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fy << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fy << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fy << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fy << endl;

 // diagonal
 for (int ix = 0; ix < nx; ix++) {
  double x = getY(ix);
  Cell *c = getCell(ix, ix, nz / 2);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  output::fdiag << setw(14) << tau << setw(14) << sqrt(2.) * x << setw(14) << vx
           << setw(14) << vy << setw(14) << e << setw(14) << nb << setw(14) << t
           << setw(14) << mub << endl;
  output::fdiag << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
           << setw(14) << c->getpi(0, 2);
  output::fdiag << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
           << setw(14) << c->getpi(1, 2);
  output::fdiag << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
           << setw(14) << c->getpi(2, 3);
  output::fdiag << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
           << c->getViscCorrCutFlag() << endl;
 }
 output::fdiag << endl;

 // Z direction
 for (int iz = 0; iz < nz; iz++) {
  double z = getZ(iz);
  Cell *c = getCell(nx / 2, ny / 2, iz);
  getCMFvariables(getCell(nx / 2, ny / 2, iz), tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  output::fz << setw(14) << tau << setw(14) << z << setw(14) << vz << setw(14) << vx
        << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14) << mub;
  output::fz << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1) << setw(14)
        << c->getpi(0, 2);
  output::fz << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1) << setw(14)
        << c->getpi(1, 2);
  output::fz << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2) << setw(14)
        << c->getpi(2, 3);
  output::fz << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
        << c->getViscCorrCutFlag() << endl;
 }
 output::fz << endl;
}

// unput: geom. rapidity + velocities in Bjorken frame, --> output: velocities
// in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz) {
 const double Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
 vx = vx * cosh(Y - eta) / cosh(Y);
 vy = vy * cosh(Y - eta) / cosh(Y);
 vz = tanh(Y);
}

void Fluid::outputSurface(double tau) {
 static double nbSurf = 0.0;
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7];
 double E = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
        vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
        txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0.,
        eps_p = 0. ;
 double eta = 0;
 int nelements = 0, nsusp = 0;  // all surface emenents and suspicious ones
 int nCoreCells = 0,
     nCoreCutCells = 0;  // cells with e>eCrit and cells with cut visc.corr.
                         //-- Cornelius: allocating memory for corner points
 double ****ccube = new double ***[2];
 for (int i1 = 0; i1 < 2; i1++) {
  ccube[i1] = new double **[2];
  for (int i2 = 0; i2 < 2; i2++) {
   ccube[i1][i2] = new double *[2];
   for (int i3 = 0; i3 < 2; i3++) {
    ccube[i1][i2][i3] = new double[2];
   }
  }
 }
//----end Cornelius
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 output::f2d << endl;
 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    Cell *c = getCell(ix, iy, iz);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    c->getQ(Q);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    double s = eos->s(e, nb, nq, ns);
    eta = getZ(iz);
    const double cosh_int = (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
    const double sinh_int = (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
    E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
             (cosh_int - tanh(vz) * sinh_int) -
         tau * p * cosh_int;
    Nb1 += Q[NB_];
    Nb2 += tau * nb * (cosh_int - tanh(vz) * sinh_int) /
           sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    //---- inf check
    if (std::isinf(E)) {
     cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
          << setw(14) << vy << setw(14) << vz << endl;
     exit(1);
    }
    //--------------
    Efull += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh(eta) - tanh(vz) * sinh(eta)) -
             tau * p * cosh(eta);
    if (trcoeff->isViscous())
     Efull +=
         tau * c->getpi(0, 0) * cosh(eta) + tau * c->getpi(0, 3) * sinh(eta);
    if (e > ecrit) {
     nCoreCells++;
     if (c->getViscCorrCutFlag() < 0.9) nCoreCutCells++;
    }
    // -- noneq. corrections to entropy flux
    const double gmumu[4] = {1., -1., -1., -1.};
    double deltas = 0.;
    if (trcoeff->isViscous())
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
       deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
    if (t > 0.02) {
     s += 1.5 * deltas / ((e + p) * t);
     S += tau * s * (cosh_int - tanh(vz) * sinh_int) /
          sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    }
    Px += tau * (e + p) * vx / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    if (iz > nz/2-3 and iz < nz/2+3) {
      vxvy_num += e * (fabs(vx) - fabs(vy));
      vxvy_den += e;
      vt_den += e / sqrt(1. - vx * vx - vy * vy);
      vt_num += e / sqrt(1. - vx * vx - vy * vy) * sqrt(vx * vx + vy * vy);
      txxyy_num += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (vx * vx - vy * vy);
      txxyy_den += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (vx * vx + vy * vy) +
                 2. * p;
    }
    pi0x_num += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                fabs(c->getpi(0, 1));
    pi0x_den += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    //----- Cornelius stuff
    double QCube[2][2][2][2][7];
    double piSquare[2][2][2][10], PiSquare[2][2][2];
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       double _p, _nb, _nq, _ns, _vx, _vy, _vz;
       Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
       cc->getPrimVar(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
       cc->getQ(QCube[1][jx][jy][jz]);
       ccube[1][jx][jy][jz] = e;
       cc->getPrimVarPrev(eos, tau - dt, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
       cc->getQprev(QCube[0][jx][jy][jz]);
       ccube[0][jx][jy][jz] = e;
       // ---- get viscous tensor
       for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj <= ii; jj++)
         piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
       PiSquare[jx][jy][jz] = cc->getPi();
      }
    cornelius->find_surface_4d(ccube);
    const int Nsegm = cornelius->get_Nelements();
    for (int isegm = 0; isegm < Nsegm; isegm++) {
     nelements++;
     output::ffreeze.precision(15);
     output::ffreeze << setw(24) << tau + cornelius->get_centroid_elem(isegm, 0)
             << setw(24) << getX(ix) + cornelius->get_centroid_elem(isegm, 1)
             << setw(24) << getY(iy) + cornelius->get_centroid_elem(isegm, 2)
             << setw(24) << getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     // ---- interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., piC[10], PiC = 0., nbC = 0.,
            nqC = 0.;  // values at the centre, to be interpolated
     double QC[7] = {0., 0., 0., 0., 0., 0., 0.};
     double eC = 0., pC = 0.;
     for (int ii = 0; ii < 10; ii++) piC[ii] = 0.0;
     double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / dt,
                        cornelius->get_centroid_elem(isegm, 0) / dt};
     double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm, 1) / dx,
                        cornelius->get_centroid_elem(isegm, 1) / dx};
     double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
                        cornelius->get_centroid_elem(isegm, 2) / dy};
     double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / dz,
                        cornelius->get_centroid_elem(isegm, 3) / dz};
     for (int jt = 0; jt < 2; jt++)
      for (int jx = 0; jx < 2; jx++)
       for (int jy = 0; jy < 2; jy++)
        for (int jz = 0; jz < 2; jz++)
         for (int i = 0; i < 7; i++) {
          QC[i] += QCube[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
     for (int i = 0; i < 7; i++)
      QC[i] = QC[i] / (tau + cornelius->get_centroid_elem(isegm, 0));
     double _ns = 0.0;
     transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     if (TC > 0.4 || fabs(mubC) > 0.85) {
      cout << "#### Error (surface): high T/mu_b ####\n";
     }
     if (eC > ecrit * 2.0 || eC < ecrit * 0.5) nsusp++;
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        for (int ii = 0; ii < 10; ii++)
         piC[ii] +=
             piSquare[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        PiC += PiSquare[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
       }
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);

     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     const double tauC = tau + cornelius->get_centroid_elem(isegm, 0);
     double dsigma[4];
     // ---- transform dsigma to lab.frame :
     const double ch = cosh(etaC);
     const double sh = sinh(etaC);
     dsigma[0] = tauC * (ch * cornelius->get_normal_elem(0, 0) -
                         sh / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[3] = tauC * (-sh * cornelius->get_normal_elem(0, 0) +
                         ch / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[1] = tauC * cornelius->get_normal_elem(0, 1);
     dsigma[2] = tauC * cornelius->get_normal_elem(0, 2);
     double dVEff = 0.0;
     for (int ii = 0; ii < 4; ii++)
      dVEff += dsigma[ii] * uC[ii];  // normalize for Delta eta=1
     vEff += dVEff;
     for (int ii = 0; ii < 4; ii++) output::ffreeze << setw(24) << dsigma[ii];
     for (int ii = 0; ii < 4; ii++) output::ffreeze << setw(24) << uC[ii];
     output::ffreeze << setw(24) << TC << setw(24) << mubC << setw(24) << muqC
             << setw(24) << musC;
#ifdef OUTPI
     double picart[10];
     /*pi00*/ picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                      2. * ch * sh * piC[index44(0, 3)] +
                                      sh * sh * piC[index44(3, 3)];
     /*pi01*/ picart[index44(0, 1)] =
         ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
     /*pi02*/ picart[index44(0, 2)] =
         ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
     /*pi03*/ picart[index44(0, 3)] =
         ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC[index44(0, 3)];
     /*pi11*/ picart[index44(1, 1)] = piC[index44(1, 1)];
     /*pi12*/ picart[index44(1, 2)] = piC[index44(1, 2)];
     /*pi13*/ picart[index44(1, 3)] =
         sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
     /*pi22*/ picart[index44(2, 2)] = piC[index44(2, 2)];
     /*pi23*/ picart[index44(2, 3)] =
         sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
     /*pi33*/ picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                      ch * ch * piC[index44(3, 3)] +
                                      2. * sh * ch * piC[index44(0, 3)];
     for (int ii = 0; ii < 10; ii++) output::ffreeze << setw(24) << picart[ii];
     output::ffreeze << setw(24) << PiC << endl;
#else
     output::ffreeze << setw(24) << dVEff << endl;
#endif
     double dEsurfVisc = 0.;
     for (int i = 0; i < 4; i++)
      dEsurfVisc += picart[index44(0, i)] * dsigma[i];
     EtotSurf += (eC + pC) * uC[0] * dVEff - pC * dsigma[0] + dEsurfVisc;
     nbSurf += nbC * dVEff;
    }
    // if(cornelius->get_Nelements()>1) cout << "oops, Nelements>1\n" ;
    //----- end Cornelius
   }
 E = E * dx * dy * dz;
 Efull = Efull * dx * dy * dz;
 S = S * dx * dy * dz;
 Nb1 *= dx * dy * dz;
 Nb2 *= dx * dy * dz;
 eps_p = txxyy_num / txxyy_den ;
 output::faniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
           << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << endl;
 cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
      << nbSurf << setw(13) << S << setw(10) << nelements << setw(10) << nsusp
      << setw(13) << (float)(nCoreCutCells) / (float)(nCoreCells) << setw(13)
      << eps_p << setw(13) << (double)vt_num/vt_den << setw(13) << (double)vxvy_num/vxvy_den << endl;
 //-- Cornelius: all done, let's free memory
 for (int i1 = 0; i1 < 2; i1++) {
  for (int i2 = 0; i2 < 2; i2++) {
   for (int i3 = 0; i3 < 2; i3++) {
    delete[] ccube[i1][i2][i3];
   }
   delete[] ccube[i1][i2];
  }
  delete[] ccube[i1];
 }
 delete[] ccube;
//----end Cornelius
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 if (nelements == 0) exit(0);
}

void Fluid::outputManiz(double tau) {
  double eps_p_num1 = 0., eps_p_den1 = 0., eps_p_num2 = 0., eps_p_den2 = 0., eps_p_num3 = 0., eps_p_den3 = 0., eps_p_num4 = 0., eps_p_den4 = 0., eps_p_num5 = 0., eps_p_den5 = 0., eps_p_num6 = 0., eps_p_den6 = 0., eps_p_num7 = 0., eps_p_den7 = 0., eps_p_num8 = 0., eps_p_den8 = 0., eps_p_num9 = 0., eps_p_den9 = 0.,
   eps_p_num10 = 0., eps_p_den10 = 0.,  eps_p_num11 = 0., eps_p_den11 = 0.,  eps_p_num12 = 0., eps_p_den12 = 0.,  eps_p_num13 = 0., eps_p_den13 = 0.,  eps_p_num14 = 0., eps_p_den14 = 0.,  eps_p_num15 = 0., eps_p_den15 = 0.,  eps_p_num16 = 0., eps_p_den16 = 0.,   eps_p_num17 = 0., eps_p_den17 = 0., eps_p_num18 = 0.,
   eps_p_den18 = 0.,  eps_p_num19 = 0., eps_p_den19 = 0.,  eps_p_num20 = 0., eps_p_den20 = 0., eps_p_num21 = 0., eps_p_den21 = 0., eps_p_num22 = 0., eps_p_den22 = 0., eps_p_num23 = 0., eps_p_den23 = 0., eps_p_num24 = 0., eps_p_den24 = 0.,  eps_p_num25 = 0., eps_p_den25 = 0.,  eps_p_num26 = 0., eps_p_den26 = 0.,
     eps_p_num27 = 0., eps_p_den27 = 0.,  eps_p_num28 = 0., eps_p_den28 = 0.,   eps_p_num29 = 0., eps_p_den29 = 0.,  eps_p_num30 = 0., eps_p_den30 = 0.,   eps_p_num31 = 0., eps_p_den31 = 0., eps_p_num32 = 0., eps_p_den32 = 0.,  eps_p_num33 = 0., eps_p_den33 = 0.,  eps_p_num34 = 0., eps_p_den34 = 0.,  eps_p_num35 = 0., eps_p_den35 = 0.,
         phi = 0., order1 = 2., order2 = 3., order3 = 4.,
        q_1 = 0., q_2 = 0., q_3 = 0., q_4 = 0., q_5 = 0., q_6 = 0., q_7 = 0., q_8 = 0., q_9 = 0., q_10 = 0.,  q_11 = 0., q_12 = 0., q_13 = 0.,  q_14 = 0., q_15 = 0., q_16 = 0., q_17 = 0. , q_18 = 0.,  q_19 = 0. , q_20 = 0., x = 0., y = 0., z = 0.; //Tomas variables
  double e, nb, nq, ns, vx, vy, vz, t, mub, muq, mus, p;
  double psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8, psi9, psi10;

  cout << "initiated Momentum Anizotropy computation as Maniz routine" << endl;


  //Space averaging of Q's
  for (int ix = 2; ix < nx - 2; ix++)
   for (int iy = 2; iy < ny - 2; iy++)
    for (int iz = 2; iz < nz - 2; iz++) {
  Cell *c = getCell(ix, iy, iz);
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
// index T^{i1} i=1
  x = getX(ix) ;
  y = getY(iy) ;
  z = getZ(iz) ;
  phi = atan2( y , x );
  //cout << phi << setw(10) << "phi" << q_1 << setw(10) << "q_1" << endl;




  q_1 += ( vx * vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) + p ) * cos( order1 * phi);
  q_2 += ( vx * vy * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * sin( order1 * phi );

  if( tau < 0.8 ){

  if( abs(z) < 0.5 )
  {
    q_3 += ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * cos( order2 * phi);
    q_4 += ( vy * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * sin( order2 * phi );
    q_5 += ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * cos( order1 * phi);
    q_6 += ( vy * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * sin( order1 * phi );
    q_7 += ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) )  ) * cos( order3 * phi);
    q_8 += ( vy * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * sin( order3 * phi );

    // Eccentricity definition by Niemi, m = 2, n = 1
    q_19 +=  ( x * x + y * y  ) * sin( 1 * phi) * e ;
    q_20 += ( x * x + y * y  ) * cos( 1 * phi) * e;

    // Eccentricity definition by Niemi, m = 2, n = 2
    q_9 +=  ( x * x + y * y  ) * sin( 2 * phi) * e ;
    q_10 += ( x * x + y * y  ) * cos( 2 * phi) * e;

    // Eccentricity definition by Niemi, m = 3, n = 2
    q_11 +=  ( x * x + y * y  ) * sqrt( x * x + y * y  )  * sin( 2 * phi) * e ;
    q_12 += ( x * x + y * y  ) * sqrt( x * x + y * y  )  * cos( 2 * phi) * e;

    // Eccentricity definition by Niemi, m = 2, n = 3
    q_13 +=  ( x * x + y * y  ) * sin( 3 * phi) * e ;
    q_14 += ( x * x + y * y  ) * cos( 3 * phi) * e;

    // Eccentricity definition by Niemi, m = 3, n = 3
    q_15 +=  ( x * x + y * y  ) * sqrt( x * x + y * y  )  * sin( 3 * phi) * e ;
    q_16 += ( x * x + y * y  ) * sqrt( x * x + y * y  )  * cos( 3 * phi) * e;

    // Eccentricity definition by Niemi, m = 2, n = 4
    q_17 +=  ( x * x + y * y  )  * sin( 4 * phi) * e ;
    q_18 += ( x * x + y * y  ) * cos( 4 * phi) * e;
  }

//cout << "just q_1" << setw(10) << q_1 << endl;
  //cout << (1. - vx * vx - vy * vy - vz*vz) << setw(10) << "gamma" << setw(10) << e << setw(10) << "energie" << setw(10) << (vx*vy*(e+p))  << setw(10) << "v_z" << setw(10) << q_1 << setw(10) << "those are Q's"  <<  setw(10) << q_2 << setw(10) << cos(order*phi) <<  setw(10) << sin(order*phi) << endl;
  }
}
  psi1 = atan2( q_2 , q_1 );
  //cout << psi << setw(10) << "<- psi" << setw(10) << q_1 << setw(10) << "<- q_1" << setw(10) << q_2 << setw(10) << "<- q_2" << endl;
  //Using phasefactor psi in space averaging of anizotropies esp_p_num, resp. esp_p_den. MidRap restriction
  psi2 = atan2( q_4 , q_3 ); //Order 2 = 3
  psi3 = atan2( q_6 , q_5 ); //Order 1 = 2
  psi4 = atan2( q_8 , q_7 ); // Order 3 = 4

  //Eccentricity psis for Neimi's definition
  psi10 = atan2( q_20 , q_19 )/2 + 3.1415/2;   //(m,n)=(2,1)
  psi5 = atan2( q_10 , q_9 )/2 + 3.1415/2;   //(m,n)=(2,2)
  psi6 = atan2( q_12 , q_11 )/2 + 3.1415/2;   //(m,n)=(3,2)
  psi7 = atan2( q_14 , q_13 )/3 + 3.1415/3;   //(m,n)=(2,3)
  psi8 = atan2( q_16 , q_15 )/3 + 3.1415/3;   //(m,n)=(3,3)
  psi9 = atan2( q_18 , q_17 )/4 + 3.1415/4;   //(m,n)=(2,4)

   std::cout << setw(15) << "Eccentricity psi's" << setw(15) << tau << setw(15) << psi5 << setw(15) << psi6 << setw(15) << psi7 << setw(15) << psi8 << setw(15) << "Eccentricity psi's" +  '\n';


    for (int ix = 2; ix < nx - 2; ix++)
     for (int iy = 2; iy < ny - 2; iy++)
      for (int iz = 2; iz < nz - 2; iz++) {
        Cell *c = getCell(ix, iy, iz);
        getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
        eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  // index T^{i1} i=1 , [vz or as upwards tanh(vz)?]
  x = getX(ix) ;
  y = getY(iy) ;
  z = getZ(iz) ;
  phi = atan2( y , x );
  //cout << "this is phi"  <<  setw(10) << phi << "this is psi"  <<  setw(10) << psi << endl;

  //Midrap restrikce
  if( abs(z) < 0.5 ) {

    //Maniz(1)[tj. T00 normalizace] m=2, n=1
    eps_p_num31 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
    ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
    (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  ) *
    cos( 1 * ( phi - psi10 ) ) ;
    eps_p_den31 += ( x * x + y * y  ) * (( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p);

  //Maniz(1)[tj. T00 normalizace] m=2, n=2
  eps_p_num1 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  ) *
  cos( order1 * ( phi - psi5 ) );
  eps_p_den1 += ( x * x + y * y  ) * (( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p);

  //Maniz(1)[tj. T00 normalizace] m=2, n=3

  eps_p_num2 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  ) *
  cos( order2 * ( phi - psi7 ) );
  eps_p_den2 += (( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p) * ( x * x + y * y  );

  //Maniz(1)[tj. T00 normalizace] m=2, n=4
  eps_p_num8  += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  ) *
  cos( order3 * ( phi - psi9 ) );
  eps_p_den8 += (( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p) * ( x * x + y * y  );

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=2, n=2

  eps_p_num3 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  ) *
  cos( order1 * ( phi - psi5 ) );
  eps_p_den3 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * ( x * x + y * y  );


  //Maniz(1)*r^2 [tj. T00*r^2 normalizace] m=1, n=2
  eps_p_num4 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  ) *
  cos( order1 * ( phi - psi5 ) );
  eps_p_den4 += ( ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p ) * sqrt( x * x + y * y  );

  //Maniz(1)*r^2 [tj. T00*r^2 normalizace] m=1, n=3
  eps_p_num11 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  ) *
  cos( order2 * ( phi - psi7 ) );
  eps_p_den11 += ( ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p ) * sqrt( x * x + y * y  );

  //Maniz(1)*r^2 [tj. T00*r^2 normalizace] m=1, n=4
  eps_p_num12 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  ) *
  cos( order3 * ( phi - psi9 ) );
  eps_p_den12 += ( ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) - p ) * sqrt( x * x + y * y  ) ;

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=2, n=3

  eps_p_num5 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  ( x * x + y * y  ) * cos( order2 * ( phi - psi7 ) );
  eps_p_den5 += ( x * x + y * y  ) * sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) );

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=2, n=4
  eps_p_num6 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  ( x * x + y * y  ) * cos( order3 * ( phi - psi9 ) );
  eps_p_den6 += ( x * x + y * y  ) * sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) );

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=2, n=1
  eps_p_num32 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  ( x * x + y * y  ) * cos( 1 * ( phi - psi10 ) );
  eps_p_den32 += ( x * x + y * y  ) * sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) );


  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=1, n=2
  eps_p_num7 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  sqrt( x * x + y * y  ) * cos( order1 * ( phi - psi5 ) );
  eps_p_den7 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  );

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=1, n=3
  eps_p_num9 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  sqrt( x * x + y * y  ) * cos( order2 * ( phi - psi7 ) );
  eps_p_den9 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  );

  //Maniz(2)[tj. T10^2+T02^2 normalizace] m=1, n=4
  eps_p_num10 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) *
  sqrt( x * x + y * y  ) * cos( order3 * ( phi - psi9 ) );
  eps_p_den10 += sqrt( ( vx * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) *
  ( vx *( e + p )/(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz))) + ( vy * ( e + p )/
  (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) * ( vy * ( e + p ) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) ) ) * sqrt( x * x + y * y  );

  //Maniz(3 - Iurii) m=2, n=1
  eps_p_num13 +=  (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  ) * cos( 1 * (  phi - psi10 ));
  eps_p_den13 += (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  );

  //Maniz(3 - Iurii) m=2, n=2
  eps_p_num13 +=  (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  ) * cos( order1 * (  phi - psi5 ));
  eps_p_den13 += (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  );

  //Maniz(3 - Iurii) m=2, n=3
  eps_p_num14 +=  (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  ) * cos( order2 *  (  phi - psi7 ) );
  eps_p_den14 += (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  );

  //Maniz(3 - Iurii) m=2, n=4
  eps_p_num15 +=  (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  ) * cos( order3 *  (  phi - psi9 ) );
  eps_p_den15 += (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * ( x * x + y * y  );

  //Maniz(3 - Iurii) m=1, n=2
  eps_p_num16 +=  ( sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * cos( order1 *  (  phi - psi5 ) );
  eps_p_den16 += (  sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );

  //Maniz(3 - Iurii) m=1, n=3
  eps_p_num17 +=  ( sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * cos( order2 *  (  phi - psi7 ) );
  eps_p_den17 += ( sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );

  //Maniz(3 - Iurii) m=1, n=4
  eps_p_num18 +=  ( sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) ) * cos( order3 *  (  phi - psi9 ) );
  eps_p_den18 += ( sqrt( x * x + y * y  ) * ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );

  //Maniz(Txx - Tyy) should be same as Maniz(3 - Iurii) m=0, n=2 [probe]
  eps_p_num19 +=  (  ( vx * vx - vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );
  eps_p_den19 += (  ( vx * vx + vy * vy ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );

  //Maniz(T0x - T0y) should be same as Maniz(3 - Iurii) m=0, n=2 [probe]
  eps_p_num20 +=  (  ( abs(vx)  - abs(vy) ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );
  eps_p_den20 += (  ( abs(vx)  + abs(vy)  ) * ( e + p ) / ( 1. - vx * vx - vy * vy - tanh(vz) * tanh(vz) ) );

  //Eccentricity definition Neimi (m,n)=(2,1)
  eps_p_num34 += ( ( x * x + y * y  ) * cos( 1 * (phi - psi10 ) ) * e );
  eps_p_den34 += ( ( x * x + y * y  ) * e);

  //Eccentricity definition Neimi (m,n)=(2,2) 26 ve vypisu
  eps_p_num21 += ( ( x * x + y * y  ) * cos( 2 * (phi - psi5 ) ) * e );
  eps_p_den21 += ( ( x * x + y * y  ) * e);

  //Eccentricity definition Neimi (m,n)=(3,2)
  eps_p_num22 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * cos( 2 * ( phi - psi6 ) ) * e );
  eps_p_den22 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e);

  //Eccentricity definition Neimi (m,n)=(2,3) 28 ve vypisu
  eps_p_num23 += ( ( x * x + y * y  ) * cos( 3 * (phi - psi7 ) ) * e );
  eps_p_den23 += ( ( x * x + y * y  ) * e);

  //Eccentricity definition Neimi (m,n)=(3,3)
  eps_p_num24 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * cos( 3 * (phi - psi8 ) ) * e );
  eps_p_den24 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e);

  //Eccentricity definition Neimi (m,n)=(2,4)
  eps_p_num29 += ( ( x * x + y * y  ) * cos( 4 * (phi - psi9 ) ) * e );
  eps_p_den29 += ( ( x * x + y * y  ) * e);

  //Eccentricity simple definition as in Huichao (m,n)=(2,2)
  eps_p_num25 += -( ( x * x + y * y  ) * e * cos( 2 * (  phi - psi5 ) ) );
  eps_p_den25 += ( ( x * x + y * y  ) * e );

  //Eccentricity simple definition as in Huichao (m,n)=(2,2)
  eps_p_num35 += -( ( x * x + y * y  ) * e * cos( 1 * (  phi - psi10 ) ) );
  eps_p_den35 += ( ( x * x + y * y  ) * e );

  //Eccentricity simple definition as in Huichao (m,n)=(3,2)
  eps_p_num26 += -( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e * cos( 2 * (  phi - psi5 ) ) );
  eps_p_den26 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e );

  //Eccentricity simple definition as in Huichao (m,n)=(2,3)
  eps_p_num27 += -( ( x * x + y * y  ) * e * cos( 3 * (  phi - psi7 ) ) );
  eps_p_den27 += ( ( x * x + y * y  ) * e );

  //Eccentricity simple definition as in Huichao (m,n)=(3,3)
  eps_p_num28 += -( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e * cos( 3 * (  phi - psi7 ) ) );
  eps_p_den28 += ( ( x * x + y * y  ) * sqrt( x * x + y * y  ) * e );

  //Eccentricity simple definition as in Huichao (m,n)=(2,4)
  eps_p_num30 += -( ( x * x + y * y  )  * e * cos( 4 * (  phi - psi9 ) ) );
  eps_p_den30 += ( ( x * x + y * y  ) * e );

}
  }


  output::maniz <<  setw(5) << tau << setw(18) << eps_p_num31/eps_p_den31 << setw(18) << eps_p_num1/eps_p_den1 << setw(18) << eps_p_num2/eps_p_den2 << setw(18) << eps_p_num8/eps_p_den8 << setw(18) << eps_p_num4/eps_p_den4 << setw(18) <<
   eps_p_num11/eps_p_den11 <<  setw(18) << eps_p_num12/eps_p_den12 << setw(18) << eps_p_num32/eps_p_den32 <<  setw(18) << eps_p_num3/eps_p_den3 << setw(18) << eps_p_num5/eps_p_den5 << setw(18) << eps_p_num6/eps_p_den6 << setw(18) <<
   eps_p_num7/eps_p_den7 << setw(18) << eps_p_num9/eps_p_den9 << setw(18) << eps_p_num10/eps_p_den10 << setw(18) << eps_p_num33/eps_p_den33  << setw(18) << eps_p_num13/eps_p_den13 << setw(18) << eps_p_num14/eps_p_den14 << setw(18) << eps_p_num15/eps_p_den15 << setw(18) <<
   eps_p_num16/eps_p_den16 << setw(18) << eps_p_num17/eps_p_den17 << setw(18) << eps_p_num18/eps_p_den18 << setw(18) << eps_p_num19/eps_p_den19 << setw(18) << eps_p_num20/eps_p_den20 << setw(18) << eps_p_num34/eps_p_den34  << setw(18) <<
   eps_p_num21/eps_p_den21 << setw(18) << eps_p_num22/eps_p_den22 << setw(18) << eps_p_num23/eps_p_den23 << setw(18) << eps_p_num24/eps_p_den24 <<  setw(18) << eps_p_num29/eps_p_den29 << setw(18) << eps_p_num35/eps_p_den35  << setw(18) << eps_p_num25/eps_p_den25 << setw(18) << eps_p_num26/eps_p_den26
   << setw(18) << eps_p_num27/eps_p_den27 << setw(18) << eps_p_num28/eps_p_den28  << setw(18) << eps_p_num30/eps_p_den30 << setw(18) << psi10 << setw(18) << psi5 << setw(18) << psi7 << setw(18) << psi9 << endl;
}


void Fluid::outputCorona(double tau) {
 static double nbSurf = 0.0;
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7];
 double E = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
        vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
        txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0.;
 double eta = 0;
 int nelements = 0;

#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 output::f2d << endl;
 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    Cell *c = getCell(ix, iy, iz);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    c->getQ(Q);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    double s = eos->s(e, nb, nq, ns);
    eta = getZ(iz);
    const double cosh_int = (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
    const double sinh_int = (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
    E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
             (cosh_int - tanh(vz) * sinh_int) -
         tau * p * cosh_int;
    Nb1 += Q[NB_];
    Nb2 += tau * nb * (cosh_int - tanh(vz) * sinh_int) /
           sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    //---- inf check
    if (std::isinf(E)) {
     cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
          << setw(14) << vy << setw(14) << vz << endl;
     exit(1);
    }
    //--------------
    Efull += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh(eta) - tanh(vz) * sinh(eta)) -
             tau * p * cosh(eta);
    if (trcoeff->isViscous())
     Efull +=
         tau * c->getpi(0, 0) * cosh(eta) + tau * c->getpi(0, 3) * sinh(eta);
    // -- noneq. corrections to entropy flux
    const double gmumu[4] = {1., -1., -1., -1.};
    double deltas = 0.;
    if (trcoeff->isViscous())
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
       deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
    if (t > 0.02) {
     s += 1.5 * deltas / ((e + p) * t);
     S += tau * s * (cosh_int - tanh(vz) * sinh_int) /
          sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    }
    Px += tau * (e + p) * vx / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
    vt_num += e / sqrt(1. - vx * vx - vy * vy) * sqrt(vx * vx + vy * vy);
    vt_den += e / sqrt(1. - vx * vx - vy * vy);
    vxvy_num += e * (fabs(vx) - fabs(vy));
    vxvy_den += e;
    txxyy_num += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (vx * vx - vy * vy);
    txxyy_den += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (vx * vx + vy * vy) +
                 2. * p;
    pi0x_num += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                fabs(c->getpi(0, 1));
    pi0x_den += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));

    //----- Cornelius stuff
    bool isCorona = true, isTail = true;
    double QCube[2][2][2][7];
    double piSquare[2][2][2][10], PiSquare[2][2][2];
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       double _p, _nb, _nq, _ns, _vx, _vy, _vz;
       Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
       cc->getPrimVar(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
       cc->getQ(QCube[jx][jy][jz]);
       if (e > ecrit) isCorona = false;
       if (e > 0.01) isTail = false;
       // ---- get viscous tensor
       for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj <= ii; jj++)
         piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
       PiSquare[jx][jy][jz] = cc->getPi();
      }
    if (isCorona && !isTail) {
     nelements++;
     output::ffreeze.precision(15);
     output::ffreeze << setw(24) << tau << setw(24) << getX(ix) + 0.5 * dx << setw(24)
             << getY(iy) + 0.5 * dy << setw(24) << getZ(iz) + 0.5 * dz;
     // ---- interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., piC[10], PiC = 0., nbC = 0.,
            nqC = 0.;  // values at the centre, to be interpolated
     double QC[7] = {0., 0., 0., 0., 0., 0., 0.};
     double eC = 0., pC = 0.;
     for (int ii = 0; ii < 10; ii++) piC[ii] = 0.0;
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++)
        for (int i = 0; i < 7; i++) {
         QC[i] += QCube[jx][jy][jz][i] * 0.125;
        }
     for (int i = 0; i < 7; i++) QC[i] = QC[i] / tau;
     double _ns = 0.0;
     transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     if (TC > 0.4 || fabs(mubC) > 0.85) {
      cout << "#### Error (surface): high T/mu_b ####\n";
     }
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        for (int ii = 0; ii < 10; ii++)
         piC[ii] += piSquare[jx][jy][jz][ii] * 0.125;
        PiC += PiSquare[jx][jy][jz] * 0.125;
       }
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = getZ(iz) + 0.5 * dz;
     transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);

     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     const double tauC = tau;
     double dsigma[4];
     // ---- transform dsigma to lab.frame :
     const double ch = cosh(etaC);
     const double sh = sinh(etaC);
     dsigma[0] = tauC * (ch * dx * dy * dz);
     dsigma[3] = tauC * (-sh * dx * dy * dz);
     dsigma[1] = 0.0;
     dsigma[2] = 0.0;
     double dVEff = 0.0;
     for (int ii = 0; ii < 4; ii++)
      dVEff += dsigma[ii] * uC[ii];  // normalize for Delta eta=1
     vEff += dVEff;
     for (int ii = 0; ii < 4; ii++) output::ffreeze << setw(24) << dsigma[ii];
     for (int ii = 0; ii < 4; ii++) output::ffreeze << setw(24) << uC[ii];
     output::ffreeze << setw(24) << TC << setw(24) << mubC << setw(24) << muqC
             << setw(24) << musC;
#ifdef OUTPI
     double picart[10];
     /*pi00*/ picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                      2. * ch * sh * piC[index44(0, 3)] +
                                      sh * sh * piC[index44(3, 3)];
     /*pi01*/ picart[index44(0, 1)] =
         ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
     /*pi02*/ picart[index44(0, 2)] =
         ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
     /*pi03*/ picart[index44(0, 3)] =
         ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC[index44(0, 3)];
     /*pi11*/ picart[index44(1, 1)] = piC[index44(1, 1)];
     /*pi12*/ picart[index44(1, 2)] = piC[index44(1, 2)];
     /*pi13*/ picart[index44(1, 3)] =
         sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
     /*pi22*/ picart[index44(2, 2)] = piC[index44(2, 2)];
     /*pi23*/ picart[index44(2, 3)] =
         sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
     /*pi33*/ picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                      ch * ch * piC[index44(3, 3)] +
                                      2. * sh * ch * piC[index44(0, 3)];
     for (int ii = 0; ii < 10; ii++) output::ffreeze << setw(24) << picart[ii];
     output::ffreeze << setw(24) << PiC << endl;
#else
     output::ffreeze << setw(24) << dVEff << endl;
#endif
     double dEsurfVisc = 0.;
     for (int i = 0; i < 4; i++)
      dEsurfVisc += picart[index44(0, i)] * dsigma[i];
     EtotSurf += (eC + pC) * uC[0] * dVEff - pC * dsigma[0] + dEsurfVisc;
     nbSurf += nbC * dVEff;
    }
    //----- end Cornelius
   }
 E = E * dx * dy * dz;
 Efull = Efull * dx * dy * dz;
 S = S * dx * dy * dz;
 Nb1 *= dx * dy * dz;
 Nb2 *= dx * dy * dz;
 output::faniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
           << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << endl;
 cout << setw(10) << "tau" << setw(13) << "E" << setw(13) << "Efull" << setw(13)
      << "Nb" << setw(13) << "Sfull" << setw(10) << "elements" << setw(10)
      << "susp." << setw(13) << "\%cut" << endl;
 cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
      << nbSurf << setw(13) << S << endl;
#ifdef SWAP_EOS
 swap(eos, eosH);
#endif
 cout << "corona elements : " << nelements << endl;
}


void Fluid::InitialAnisotropies(double tau0) {
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;

 double xcm = 0.0, xcm_nom = 0.0, xcm_denom = 0.0 ;
 double ycm = 0.0, ycm_nom = 0.0, ycm_denom = 0.0 ;

 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   for (int iz = 0; iz < nz; iz++) {
    if(fabs(getZ(iz)) < 0.5) {
     Cell* c = getCell(ix, iy, iz);
     double x = getX(ix) ;
     double y = getY(iy) ;
     getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
     xcm_nom += x * e ;
     xcm_denom += e ;
     ycm_nom += y * e ;
     ycm_denom += e ;
    }
   }
  }
 }

 xcm = xcm_nom / xcm_denom ;
 ycm = ycm_nom / ycm_denom ;

 double eps2 = 0.0, eps2_nom_real = 0.0, eps2_nom_imag = 0.0, eps2_denom = 0.0 ;
 double eps3 = 0.0, eps3_nom_real = 0.0, eps3_nom_imag = 0.0, eps3_denom = 0.0 ;

 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   for (int iz = 0; iz < nz; iz++) {
    if(fabs(getZ(iz)) < 0.5) {
     Cell* c = getCell(ix, iy, iz);
     double x = getX(ix) ;
     double y = getY(iy) ;
     x = x - xcm ;
     y = y - ycm ;
     double r = sqrt(x*x + y*y) ;
     double phi = atan2(y, x) ;
     getCMFvariables(c, tau0, e, nb, nq, ns, vx, vy, vz);
     eps2_denom += pow(r, 2) * e ;
     eps2_nom_real += pow(r, 2) * cos(2 * phi) * e ;
     eps2_nom_imag += pow(r, 2) * sin(2 * phi) * e ;
     eps3_denom += pow(r, 3) * e ;
     eps3_nom_real += pow(r, 3) * cos(3 * phi) * e ;
     eps3_nom_imag += pow(r, 3) * sin(3 * phi) * e ;
    }
   }
  }
 }

 eps2 = sqrt(pow(eps2_nom_real, 2) + pow(eps2_nom_imag, 2)) / eps2_denom ;
 eps3 = sqrt(pow(eps3_nom_real, 2) + pow(eps3_nom_imag, 2)) / eps3_denom ;

 cout << "epsilon2 = " << eps2 << endl;
 cout << "epsilon3 = " << eps3 << endl;

 exit(1) ;
}
