/*
 *
 * This file accompanied with the header file specfun.h is a partial
 * C translation of the Fortran code by Zhang and Jin following
 * original description:
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *       COMPUTATION OF SPECIAL FUNCTIONS
 *
 *          Shanjie Zhang and Jianming Jin
 *
 *       Copyrighted but permission granted to use code in programs.
 *       Buy their book:
 *
 *          Shanjie Zhang, Jianming Jin,
 *          Computation of Special Functions,
 *          Wiley, 1996,
 *          ISBN: 0-471-11963-6,
 *          LC: QA351.C45.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *       Scipy changes:
 *       - Compiled into a single source file and changed REAL To DBLE throughout.
 *       - Changed according to ERRATA.
 *       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
 *       - Made functions return sf_error codes in ISFER variables instead
 *         of printing warnings. The codes are
 *         - SF_ERROR_OK        = 0: no error
 *         - SF_ERROR_SINGULAR  = 1: singularity encountered
 *         - SF_ERROR_UNDERFLOW = 2: floating point underflow
 *         - SF_ERROR_OVERFLOW  = 3: floating point overflow
 *         - SF_ERROR_SLOW      = 4: too many iterations required
 *         - SF_ERROR_LOSS      = 5: loss of precision
 *         - SF_ERROR_NO_RESULT = 6: no result obtained
 *         - SF_ERROR_DOMAIN    = 7: out of domain
 *         - SF_ERROR_ARG       = 8: invalid input parameter
 *         - SF_ERROR_OTHER     = 9: unclassified error
 *       - Improved initial guesses for roots in JYZO.
 *
 *
 */

/*
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "../config.h"
#include <memory>

namespace xsf {
namespace specfun {

    // The Status enum is the return type of a few private, low-level functions
    // defined here.  Currently the only use is by functions that allocate
    // memory internally.  If the allocation fails, the function returns
    // Status::NoMemory.

    enum class Status { OK = 0, NoMemory, Other };

    void airyb(double, double *, double *, double *, double *);
    void bjndd(double, int, double *, double *, double *);

    void cerzo(int, std::complex<double> *);

    void cyzo(int, int, int, std::complex<double> *, std::complex<double> *);

    void cerf(std::complex<double>, std::complex<double> *, std::complex<double> *);
    std::complex<double> cgama(std::complex<double>, int);
    double chgubi(double, double, double, int *);
    double chguit(double, double, double, int *);
    double chgul(double, double, double, int *);
    double chgus(double, double, double, int *);
    void cpbdn(int, std::complex<double>, std::complex<double> *, std::complex<double> *);
    std::complex<double> cpdla(int, std::complex<double>);
    std::complex<double> cpdsa(int, std::complex<double>);
    double cv0(double, double, double);
    double cvf(int, int, double, double, int);
    double cvql(int, int, double);
    double cvqm(int, double);
    double gaih(double);
    double gam0(double);
    double gamma2(double);

    template <typename T>
    void jynbh(int, int, T, int *, T *, T *);

    void jyndd(int, double, double *, double *, double *, double *, double *, double *);

    double lpmv0(double, int, double);
    int msta1(double, int);
    int msta2(double, int, int);
    double psi_spec(double);
    double refine(int, int, double, double);

    template <typename T>
    void sckb(int, int, T, T *, T *);

    template <typename T>
    Status sdmn(int, int, T, T, int, T *);

    template <typename T>
    void sphj(T, int, int *, T *, T *);

    template <typename T>
    void sphy(T, int, int *, T *, T *);

    template <typename T>
    Status aswfa(T x, int m, int n, T c, int kd, T cv, T *s1f, T *s1d) {

        // ===========================================================
        // Purpose: Compute the prolate and oblate spheroidal angular
        //          functions of the first kind and their derivatives
        // Input :  m  --- Mode parameter,  m = 0,1,2,...
        //          n  --- Mode parameter,  n = m,m+1,...
        //          c  --- Spheroidal parameter
        //          x  --- Argument of angular function, |x| < 1.0
        //          KD --- Function code
        //                 KD=1 for prolate;  KD=-1 for oblate
        //          cv --- Characteristic value
        // Output:  S1F --- Angular function of the first kind
        //          S1D --- Derivative of the angular function of
        //                  the first kind
        // Return value:
        //          Status::OK
        //              Normal return.
        //          Status::NoMemory
        //              An internal memory allocation failed.
        //
        // Routine called:
        //          SDMN for computing expansion coefficients df
        //          SCKB for computing expansion coefficients ck
        // ===========================================================

        int ip, k, nm, nm2;
        T a0, d0, d1, r, su1, su2, x0, x1;
        auto ck = std::unique_ptr<T[]>{new (std::nothrow) T[200]()};
        auto df = std::unique_ptr<T[]>{new (std::nothrow) T[200]()};
        if (ck == nullptr || df == nullptr) {
            return Status::NoMemory;
        }
        const T eps = 1e-14;
        x0 = x;
        x = fabs(x);
        ip = ((n - m) % 2 == 0 ? 0 : 1);
        nm = 40 + (int)((n - m) / 2 + c);
        nm2 = nm / 2 - 2;
        if (sdmn(m, n, c, cv, kd, df.get()) == Status::NoMemory) {
            return Status::NoMemory;
        }
        sckb(m, n, c, df.get(), ck.get());
        x1 = 1.0 - x * x;
        if ((m == 0) && (x1 == 0.0)) {
            a0 = 1.0;
        } else {
            a0 = pow(x1, 0.5 * m);
        }
        su1 = ck[0];
        for (k = 1; k <= nm2; k++) {
            r = ck[k] * pow(x1, k);
            su1 += r;
            if ((k >= 10) && (fabs(r / su1) < eps)) {
                break;
            }
        }
        *s1f = a0 * pow(x, ip) * su1;
        if (x == 1.0) {
            if (m == 0) {
                *s1d = ip * ck[0] - 2.0 * ck[1];
            } else if (m == 1) {
                *s1d = -1e100;
            } else if (m == 2) {
                *s1d = -2.0 * ck[0];
            } else if (m >= 3) {
                *s1d = 0.0;
            }
        } else {
            d0 = ip - m / x1 * pow(x, ip + 1.0);
            d1 = -2.0 * a0 * pow(x, ip + 1.0);
            su2 = ck[1];
            for (k = 2; k <= nm2; k++) {
                r = static_cast<T>(k) * ck[k] * pow(x1, (k - 1.0));
                su2 += r;
                if ((k >= 10) && (fabs(r / su2) < eps)) {
                    break;
                }
            }
            *s1d = d0 * a0 * su1 + d1 * su2;
        }
        if ((x0 < 0.0) && (ip == 0)) {
            *s1d = -*s1d;
        }
        if ((x0 < 0.0) && (ip == 1)) {
            *s1f = -*s1f;
        }
        x = x0;
        return Status::OK;
    }
