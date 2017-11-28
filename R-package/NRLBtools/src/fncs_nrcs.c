/*<<<
 * ------------------------------------------------------------------------
 * Copyright (C) 2006-2017 The Trustees of Columbia University in the City
 * of New York. All rights reserved.
 *
 * The REDUCE Suite of software programs has been developed and actively
 * maintained by Dr. Xiang-Jun Lu at the Bussemaker Laboratory at the
 * Department of Biological Sciences, Columbia University. The suite has
 * its origin from the linear regression-based REDUCE algorithm of
 * Dr. Harmen Bussemaker to discover cis-regulatory elements. Dr. Barret
 * Foat extended REDUCE by adding an optimization procedure to get the
 * so-called Position Specific Affinity Matrix, and implemented his
 * algorithms in the MatrixREDUCE pipeline using Perl and GNU Scientific
 * Library. Dr. Lu has completely rewritten and significantly enhanced the
 * REDUCE/MatrixREDUCE codebase using ANSI C to make the programs efficient
 * and self-contained. The REDUCE Suite (now at v2.2) consists of a dozen
 * standalone, yet interconnected programs.
 *
 * This suite of software programs is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Click
 * Software License file (at the $REDUCE_SUITE/license directory, in either
 * Word, text, or PDF format) for more details.
 *
 * Please check http://reducesuite.bussemakerlab.org/ for more information.
 * ------------------------------------------------------------------------
>>>*/

#include "REDUCE_Suite.h"

void NRC_avevar(double data[], long n, double *ave, double *var)
{
    long j;
    double s, ep;

    for (*ave = 0.0, j = 1; j <= n; j++)
        *ave += data[j];
    *ave /= n;

    *var = ep = 0.0;
    for (j = 1; j <= n; j++) {
        s = data[j] - (*ave);
        ep += s;
        *var += s * s;
    }

    *var = (*var - ep * ep / n) / (n - 1);
}

double NRC_betai(double a, double b, double x)
/* both-sides, corresponding to R: (1 - pt(fabs(tval), df)) * 2 */
{
    double bt;

    if (x < 0.0 || x > 1.0)
        nrerror("bad x in routine betai");

    if (x == 0.0 || x == 1.0)
        bt = 0.0;
    else
        bt = exp(NRC_gammln(a + b) - NRC_gammln(a) - NRC_gammln(b) +
                 a * log(x) + b * log(1.0 - x));

    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * NRC_betacf(a, b, x) / a;
    else
        return 1.0 - bt * NRC_betacf(b, a, 1.0 - x) / b;
}

double NRC_gammln(double xx)
{
    long j;
    double x, y, tmp, ser;
    static double cof[6] = { 76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5
    };

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;

    for (j = 0; j <= 5; j++)
        ser += cof[j] / ++y;

    return -tmp + log(2.5066282746310005 * ser / x);
}

double NRC_betacf(double a, double b, double x)
{
    long m, m2;
    double aa, c, d, del, h, qab, qam, qap;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;

    c = 1.0;
    d = 1.0 - qab * x / qap;

    if (fabs(d) < NR_FPMIN)
        d = NR_FPMIN;
    d = 1.0 / d;
    h = d;

    for (m = 1; m <= ITMAX; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));

        d = 1.0 + aa * d;
        if (fabs(d) < NR_FPMIN)
            d = NR_FPMIN;

        c = 1.0 + aa / c;
        if (fabs(c) < NR_FPMIN)
            c = NR_FPMIN;

        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

        d = 1.0 + aa * d;
        if (fabs(d) < NR_FPMIN)
            d = NR_FPMIN;

        c = 1.0 + aa / c;
        if (fabs(c) < NR_FPMIN)
            c = NR_FPMIN;

        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < NR_EPS)
            break;
    }

    if (m > ITMAX)
        nrerror("a or b too big, or ITMAX too small in NRC_betacf()");

    return h;
}

double NRC_erff(double x)
{
    return (x < 0.0) ? -NRC_gammp(0.5, x * x) : NRC_gammp(0.5, x * x);
}

double NRC_gammp(double a, double x)
{
    double gamser = 0, gammcf = 0, gln = 0;

    if (x < 0.0 || a <= 0.0)
        nrerror("invalid arguments in routine NRC_gammp()");

    if (x < (a + 1.0)) {
        NRC_gser(&gamser, a, x, &gln);
        return gamser;

    } else {
        NRC_gcf(&gammcf, a, x, &gln);
        return 1.0 - gammcf;
    }
}

double NRC_gammq(double a, double x)
{
    double gamser = 0, gammcf = 0, gln = 0;

    if (x < 0.0 || a <= 0.0)
        nrerror("invalid arguments in routine NRC_gammq()");

    if (x < (a + 1.0)) {
        NRC_gser(&gamser, a, x, &gln);
        return 1.0 - gamser;

    } else {
        NRC_gcf(&gammcf, a, x, &gln);
        return gammcf;
    }
}

void NRC_gcf(double *gammcf, double a, double x, double *gln)
{
    long i;
    double an, b, c, d, del, h;

    *gln = NRC_gammln(a);

    b = x + 1.0 - a;
    c = 1.0 / NR_FPMIN;
    d = 1.0 / b;
    h = d;

    for (i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;

        d = an * d + b;
        if (fabs(d) < NR_FPMIN)
            d = NR_FPMIN;

        c = b + an / c;
        if (fabs(c) < NR_FPMIN)
            c = NR_FPMIN;

        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < NR_EPS)
            break;
    }

    if (i > ITMAX)
        nrerror("a too large, ITMAX too small in NRC_gcf()");

    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}

void NRC_gser(double *gamser, double a, double x, double *gln)
{
    long n;
    double sum, del, ap;

    *gln = NRC_gammln(a);
    if (x <= 0.0) {
        if (x < 0.0)
            nrerror("x less than 0 in routine gser");
        *gamser = 0.0;
        return;

    } else {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1; n <= ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * NR_EPS) {
                *gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine NRC_gser()");
        return;
    }
}

void NRC_fit(double x[], double y[], long ndata, double *a, double *b, double *siga,
             double *sigb, double *chisq)
{
    long i;
    double t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat, tmp;

    *b = 0.0;
    for (i = 1; i <= ndata; i++) {
        sx += x[i];
        sy += y[i];
    }
    ss = ndata;

    sxoss = sx / ss;
    for (i = 1; i <= ndata; i++) {
        t = x[i] - sxoss;
        st2 += t * t;
        *b += t * y[i];
    }

    *b /= st2;
    *a = (sy - sx * (*b)) / ss;
    *siga = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    *sigb = sqrt(1.0 / st2);

    *chisq = 0.0;
    for (i = 1; i <= ndata; i++) {
        tmp = y[i] - (*a) - (*b) * x[i];
        *chisq += tmp * tmp;
    }

    sigdat = sqrt((*chisq) / (ndata - 2));
    *siga *= sigdat;
    *sigb *= sigdat;
}

void NRC_lfit(double **X, double *y, long ndat, long ma, double *a, double **covar,
              double *chisq)
{
    long i, *ia, j, k, l, m, mfit = 0;
    double ym, wt, sum, tmp, **beta;

    ia = lvector(1, ma);
    init_lvector(ia, 1, ma, 1);  /* fit all parameters */

    mfit = ma;
    if (mfit == 0)
        nrerror("NRC_lfit(): no parameters to fit");

    beta = dmatrix(1, ma, 1, 1);
    for (j = 1; j <= mfit; j++) {
        for (k = 1; k <= mfit; k++)
            covar[j][k] = 0.0;
        beta[j][1] = 0.0;
    }

    for (i = 1; i <= ndat; i++) {
        ym = y[i];
        if (mfit < ma) {
            for (j = 1; j <= ma; j++)
                if (!ia[j])
                    ym -= a[j] * X[i][j];
        }

        for (j = 0, l = 1; l <= ma; l++) {
            if (ia[l]) {
                wt = X[i][l];
                for (j++, k = 0, m = 1; m <= l; m++)
                    if (ia[m])
                        covar[j][++k] += wt * X[i][m];
                beta[j][1] += ym * wt;
            }
        }
    }

    for (j = 2; j <= mfit; j++)
        for (k = 1; k < j; k++)
            covar[k][j] = covar[j][k];

    NRC_gaussj(covar, mfit, beta, 1);
    for (j = 0, l = 1; l <= ma; l++)
        if (ia[l])
            a[l] = beta[++j][1];

    *chisq = 0.0;
    for (i = 1; i <= ndat; i++) {
        for (sum = 0.0, j = 1; j <= ma; j++)
            sum += a[j] * X[i][j];
        tmp = y[i] - sum;
        *chisq += tmp * tmp;
    }

    NRC_covsrt(covar, ma, ia, mfit);

    free_dmatrix(beta, 1, DUMMY, 1, DUMMY);
    free_lvector(ia, 1, DUMMY);
}

void NRC_my_svdfit(double **X, double y[], long ndat, long ma, double a[],
                   double **covar, double *chisq)
/* customized NRC_svdfit(), called in the same style as NRC_lfit() above */
{
    double **u, **v, *w;

    u = dmatrix(1, ndat, 1, ma);
    v = dmatrix(1, ma, 1, ma);
    w = dvector(1, ma);

    NRC_svdfit(X, y, ndat, ma, a, u, v, w, chisq);
    NRC_svdvar(v, ma, w, covar);

    free_dmatrix(u, 1, DUMMY, 1, DUMMY);
    free_dmatrix(v, 1, DUMMY, 1, DUMMY);
    free_dvector(w, 1, DUMMY);
}

void NRC_svdfit(double **X, double y[], long ndat, long ma, double a[],
                double **u, double **v, double w[], double *chisq)
{
    long j, i;
    double tmp, dsum;

    copy_dmatrix(X, 1, ndat, 1, ma, u);

    NRC_svdcmp(u, ndat, ma, w, v);
    NRC_svdzwi(w, ma);
    NRC_svbksb(u, w, v, ndat, ma, y, a);

    *chisq = 0.0;
    for (i = 1; i <= ndat; i++) {
        dsum = 0.0;
        for (j = 1; j <= ma; j++)
            dsum += a[j] * X[i][j];
        tmp = y[i] - dsum;
        *chisq += tmp * tmp;
    }
}

void NRC_svdvar(double **v, long ma, double w[], double **covar)
/* called after NRC_svdfit(), with w[] zeroed out ... */
{
    long i, j, k;
    double sum, *wti;

    wti = dvector(1, ma);

    for (i = 1; i <= ma; i++) {
        wti[i] = 0.0;
        if (w[i])
            wti[i] = 1.0 / (w[i] * w[i]);
    }

    for (i = 1; i <= ma; i++) {
        for (j = 1; j <= i; j++) {
            for (sum = 0.0, k = 1; k <= ma; k++)
                sum += v[i][k] * v[j][k] * wti[k];
            covar[j][i] = covar[i][j] = sum;
        }
    }

    free_dvector(wti, 1, DUMMY);
}

void NRC_gaussj(double **a, long n, double **b, long m)
{
    long i, icol = 0, irow = 0, j, k, l, ll;
    long *indxc, *indxr, *ipiv;
    double big, dum, pivinv;

    indxc = lvector(1, n);
    indxr = lvector(1, n);
    ipiv = lvector(1, n);

    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if (ipiv[j] != 1)
                for (k = 1; k <= n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }

        ++(ipiv[icol]);
        if (irow != icol) {
            for (l = 1; l <= n; l++)
                dval_swap(&a[irow][l], &a[icol][l]);
            for (l = 1; l <= m; l++)
                dval_swap(&b[irow][l], &b[icol][l]);
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] <= TOLX)
            nrerror("NRC_gaussj(): singular matrix");

        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++)
            a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++)
            b[icol][l] *= pivinv;

        for (ll = 1; ll <= n; ll++)
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++)
                    a[ll][l] -= a[icol][l] * dum;
                for (l = 1; l <= m; l++)
                    b[ll][l] -= b[icol][l] * dum;
            }
    }

    for (l = n; l >= 1; l--) {
        if (indxr[l] != indxc[l])
            for (k = 1; k <= n; k++)
                dval_swap(&a[k][indxr[l]], &a[k][indxc[l]]);
    }

    free_lvector(ipiv, 1, DUMMY);
    free_lvector(indxr, 1, DUMMY);
    free_lvector(indxc, 1, DUMMY);
}

void NRC_covsrt(double **covar, long ma, long ia[], long mfit)
{
    long i, j, k;

    for (i = mfit + 1; i <= ma; i++)
        for (j = 1; j <= i; j++)
            covar[i][j] = covar[j][i] = 0.0;

    k = mfit;
    for (j = ma; j >= 1; j--) {
        if (ia[j]) {
            for (i = 1; i <= ma; i++)
                dval_swap(&covar[i][k], &covar[i][j]);
            for (i = 1; i <= ma; i++)
                dval_swap(&covar[k][i], &covar[j][i]);
            k--;
        }
    }
}

void NRC_pearsn(double x[], double y[], long n, double *r, double *prob, double *z,
                double *tval)
{
    long j;
    double yt, xt, t, df;
    double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;

    for (j = 1; j <= n; j++) {
        ax += x[j];
        ay += y[j];
    }

    ax /= n;
    ay /= n;

    for (j = 1; j <= n; j++) {
        xt = x[j] - ax;
        yt = y[j] - ay;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }

    *r = sxy / (sqrt(sxx * syy) + NR_EPS);
    *z = 0.5 * log((1.0 + (*r) + NR_EPS) / (1.0 - (*r) + NR_EPS));

    df = n - 2;
    t = (*r) * sqrt(df / ((1.0 - (*r) + NR_EPS) * (1.0 + (*r) + NR_EPS)));

    *prob = NRC_betai(0.5 * df, 0.5, df / (df + t * t));
    *tval = t;
}

void NRC_ludcmp(double **a, long n, long *indx, double *d)
{
    long i, imax = 0, j, k;
    double big, dum, sum, temp;
    double *vv;

    vv = dvector(1, n);

    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        if (big <= TOLX)
            nrerror("singular matrix in routine NRC_ludcmp()");
        vv[i] = 1.0 / big;
    }

    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }

        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = NR_EPS;
        if (j != n) {
            dum = 1.0 / (a[j][j]);
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }

    free_dvector(vv, 1, DUMMY);
}

void NRC_lubksb(double **a, long n, long *indx, double b[])
{
    long i, ii = 0, ip, j;
    double sum;

    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];

        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }

    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void NRC_lu_inv(double **a, long n, double **y)
{
    long i, j, *indx;
    double d, *col;

    col = dvector(1, n);
    indx = lvector(1, n);

    NRC_ludcmp(a, n, indx, &d);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            col[i] = 0.0;
        col[j] = 1.0;

        NRC_lubksb(a, n, indx, col);

        for (i = 1; i <= n; i++)
            y[i][j] = col[i];
    }

    free_lvector(indx, 1, DUMMY);
    free_dvector(col, 1, DUMMY);
}

void NRC_mprove(double **a, double **alud, long n, long indx[], double b[], double x[])
{
    long i, j;
    double sdp, *r;

    r = dvector(1, n);

    for (i = 1; i <= n; i++) {
        sdp = -b[i];
        for (j = 1; j <= n; j++)
            sdp += a[i][j] * x[j];
        r[i] = sdp;
    }

    NRC_lubksb(alud, n, indx, r);

    for (i = 1; i <= n; i++)
        x[i] -= r[i];

    free_dvector(r, 1, DUMMY);
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void NRC_svdcmp(double **a, long m, long n, double w[], double **v)
{
    long flag, i, its, j, jj, k, l = -1, nm = -1;
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

    rv1 = dvector(1, n);
    g = scale = anorm = 0.0;

    for (i = 1; i <= n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;

        if (i <= m) {
            for (k = i; k <= m; k++)
                scale += fabs(a[k][i]);
            if (scale) {
                for (k = i; k <= m; k++) {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = i; k <= m; k++)
                        s += a[k][i] * a[k][j];
                    f = s / h;
                    for (k = i; k <= m; k++)
                        a[k][j] += f * a[k][i];
                }
                for (k = i; k <= m; k++)
                    a[k][i] *= scale;
            }
        }

        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n) {
            for (k = l; k <= n; k++)
                scale += fabs(a[i][k]);
            if (scale) {
                for (k = l; k <= n; k++) {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; k++)
                    rv1[k] = a[i][k] / h;
                for (j = l; j <= m; j++) {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[j][k] * a[i][k];
                    for (k = l; k <= n; k++)
                        a[j][k] += s * rv1[k];
                }
                for (k = l; k <= n; k++)
                    a[i][k] *= scale;
            }
        }
        anorm = dval_max(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    for (i = n; i >= 1; i--) {
        if (i < n) {
            if (g) {
                for (j = l; j <= n; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[i][k] * v[k][j];
                    for (k = l; k <= n; k++)
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j <= n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    for (i = lval_min(m, n); i >= 1; i--) {
        l = i + 1;
        g = w[i];
        for (j = l; j <= n; j++)
            a[i][j] = 0.0;
        if (g) {
            g = 1.0 / g;
            for (j = l; j <= n; j++) {
                for (s = 0.0, k = l; k <= m; k++)
                    s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for (k = i; k <= m; k++)
                    a[k][j] += f * a[k][i];
            }
            for (j = i; j <= m; j++)
                a[j][i] *= g;
        } else
            for (j = i; j <= m; j++)
                a[j][i] = 0.0;
        ++a[i][i];
    }

    for (k = n; k >= 1; k--) {
        for (its = 1; its <= ITMAX; its++) {
            flag = 1;
            for (l = k; l >= 1; l--) {
                nm = l - 1;
                if ((double) (fabs(rv1[l]) + anorm) == anorm) {
                    flag = 0;
                    break;
                }
                if ((double) (fabs(w[nm]) + anorm) == anorm)
                    break;
            }

            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((double) (fabs(f) + anorm) == anorm)
                        break;
                    g = w[i];
                    h = NRC_pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++) {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }

            z = w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j = 1; j <= n; j++)
                        v[j][k] = -v[j][k];
                }
                break;
            }

            if (its == ITMAX)
                fatal("no convergence in ITMAX (%ld) NRC_svdcmp() iterations\n", ITMAX);

            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = NRC_pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;

            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = NRC_pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++) {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = NRC_pythag(f, h);
                w[j] = z;
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 1; jj <= m; jj++) {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }

            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    free_dvector(rv1, 1, DUMMY);
}

void NRC_svbksb(double **u, double w[], double **v, long m, long n, double b[],
                double x[])
/* assuming NRC_svdzwi() has been called to zero-out small singular values */
{
    long jj, j, i;
    double s, *tmp;

    tmp = dvector(1, n);

    for (j = 1; j <= n; j++) {
        s = 0.0;
        if (w[j]) {
            for (i = 1; i <= m; i++)
                s += u[i][j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }

    for (j = 1; j <= n; j++) {
        s = 0.0;
        for (jj = 1; jj <= n; jj++)
            s += v[j][jj] * tmp[jj];
        x[j] = s;
    }

    free_dvector(tmp, 1, DUMMY);
}

void NRC_svdzwi(double *w, long n)
/* zero out singular values */
{
    long i;
    double dval = -1.0;

    for (i = 1; i <= n; i++)
        if (w[i] > dval)
            dval = w[i];

    dval *= TOLX;  /* relative to the largest */

    for (i = 1; i <= n; i++)
        if (w[i] < dval)
            w[i] = 0.0;
}

void NRC_svd_inv(double **u, double w[], double **v, long m, long n, double **ginv)
/* assuming NRC_svdzwi() has been called to zero-out small singular values */
{
    long i, j, k;

    for (i = 1; i <= n; i++)
        for (j = 1; j <= m; j++) {
            ginv[i][j] = 0.0;
            for (k = 1; k <= n; k++)
                if (w[k])
                    ginv[i][j] += v[i][k] * (u[j][k] / w[k]);
        }
}

double NRC_pythag(double a, double b)
{
    double absa, absb;

    absa = fabs(a);
    absb = fabs(b);

    if (absa > absb)
        return absa * sqrt(1.0 + dval_sqr(absb / absa));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + dval_sqr(absa / absb)));
}

long NRC_chol_posdef(double **a, long n, double p[])
/* this routine also has the functionality of NRC_choldc() */
{
    long i, j, k;
    double sum;

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            for (sum = a[i][j], k = i - 1; k >= 1; k--)
                sum -= a[i][k] * a[j][k];
            if (i == j) {
                if (sum <= TOLX)  /* 0.0 causing trouble! */
                    return FALSE;
                p[i] = sqrt(sum);
            } else
                a[j][i] = sum / p[i];
        }
    }

    return TRUE;
}

void NRC_choldc(double **a, long n, double p[])
{
    long i, j, k;
    double sum;

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            for (sum = a[i][j], k = i - 1; k >= 1; k--)
                sum -= a[i][k] * a[j][k];
            if (i == j) {
                if (sum <= 0.0)
                    nrerror("NRC_choldc() failed");
                p[i] = sqrt(sum);
            } else
                a[j][i] = sum / p[i];
        }
    }
}

void NRC_cholsl(double **a, long n, double p[], double b[], double x[])
{
    long i, k;
    double sum;

    for (i = 1; i <= n; i++) {
        for (sum = b[i], k = i - 1; k >= 1; k--)
            sum -= a[i][k] * x[k];
        x[i] = sum / p[i];
    }

    for (i = n; i >= 1; i--) {
        for (sum = x[i], k = i + 1; k <= n; k++)
            sum -= a[k][i] * x[k];
        x[i] = sum / p[i];
    }
}

void NRC_chol_inv(double **a, long n, double p[], double **Ainv)
{
    long i, j;
    double *b, *x;

    b = dvector(1, n);
    x = dvector(1, n);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            b[i] = 0.0;
        b[j] = 1.0;

        NRC_cholsl(a, n, p, b, x);

        for (i = 1; i <= n; i++)
            Ainv[i][j] = x[i];
    }

    free_dvector(b, 1, DUMMY);
    free_dvector(x, 1, DUMMY);
}

void NRC_chol_L(double **a, long n, double p[], double **L)
{
    long i, j;

    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            L[i][j] = ((i > j) ? a[i][j] : (i == j ? p[i] : 0.0));
}

void NRC_chol_Linv(double **a, long n, double p[], double **Linv)
{
    long i, j, k;
    double sum;

    NRC_chol_L(a, n, p, Linv);

    for (i = 1; i <= n; i++) {
        Linv[i][i] = 1.0 / p[i];
        for (j = i + 1; j <= n; j++) {
            sum = 0.0;
            for (k = i; k < j; k++)
                sum -= Linv[j][k] * Linv[k][i];
            Linv[j][i] = sum / p[j];
        }
    }
}

void NRC_chol_Ainv(double **a, long n, double p[], double **Ainv)
/* A = L * L'; A-inv = L'-inv * L-inv; L'-inv = (L-inv)'
 * the covariance matrix in ls-fitting can be solved this way */
{
    long i, j, k;
    double sum, **Linv;

    Linv = dmatrix(1, n, 1, n);

    NRC_chol_Linv(a, n, p, Linv);

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= i; j++) {
            sum = 0.0;
            for (k = i; k <= n; k++)
                sum += Linv[k][i] * Linv[k][j];
            Ainv[i][j] = sum;
        }
    }

    for (i = 1; i <= n; i++)
        for (j = i + 1; j <= n; j++)
            Ainv[i][j] = Ainv[j][i];

    free_dmatrix(Linv, 1, DUMMY, 1, DUMMY);
}

void NRC_qrdcmp(double **a, long n, double *c, double *d, long *sing)
{
    long i, j, k;
    double scale, sigma, sum, tau;

    *sing = 0;
    for (k = 1; k < n; k++) {
        scale = 0.0;
        for (i = k; i <= n; i++)
            scale = dval_max(scale, fabs(a[i][k]));

        if (scale == 0.0) {
            *sing = 1;
            c[k] = d[k] = 0.0;

        } else {
            for (i = k; i <= n; i++)
                a[i][k] /= scale;
            for (sum = 0.0, i = k; i <= n; i++)
                sum += a[i][k] * a[i][k];

            sigma = SIGN(sqrt(sum), a[k][k]);
            a[k][k] += sigma;
            c[k] = sigma * a[k][k];
            d[k] = -scale * sigma;

            for (j = k + 1; j <= n; j++) {
                for (sum = 0.0, i = k; i <= n; i++)
                    sum += a[i][k] * a[i][j];
                tau = sum / c[k];
                for (i = k; i <= n; i++)
                    a[i][j] -= tau * a[i][k];
            }
        }
    }

    d[n] = a[n][n];
    if (d[n] == 0.0)
        *sing = 1;
}

#undef SIGN

static void NRC_rsolv(double **a, long n, double d[], double b[])
{
    long i, j;
    double sum;

    b[n] /= d[n];

    for (i = n - 1; i >= 1; i--) {
        for (sum = 0.0, j = i + 1; j <= n; j++)
            sum += a[i][j] * b[j];
        b[i] = (b[i] - sum) / d[i];
    }
}

void NRC_qrsolv(double **a, long n, double c[], double d[], double b[])
{
    long i, j;
    double sum, tau;

    for (j = 1; j < n; j++) {
        for (sum = 0.0, i = j; i <= n; i++)
            sum += a[i][j] * b[i];
        tau = sum / c[j];

        for (i = j; i <= n; i++)
            b[i] -= tau * a[i][j];
    }

    NRC_rsolv(a, n, d, b);
}

void NRC_qr_inv(double **a, long n, double c[], double d[], double **ainv)
{
    long i, j;
    double *b;

    b = dvector(1, n);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            b[i] = 0.0;
        b[j] = 1.0;

        NRC_qrsolv(a, n, c, d, b);

        for (i = 1; i <= n; i++)
            ainv[i][j] = b[i];
    }

    free_dvector(b, 1, DUMMY);
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double NRC_ran1(long *idum)
{
    long j, k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0)
                *idum += IM;
            if (j < NTAB)
                iv[j] = *idum;
        }

        iy = iv[0];
    }

    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;

    if (*idum < 0)
        *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;

    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double NRC_ran2(long *idum)
{
    long j, k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);

        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;

    if (*idum < 0)
        *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;

    if (idum2 < 0)
        idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;

    if (iy < 1)
        iy += IMM1;
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double NRC_ran3(long *idum)
{
    static long inext, inextp;
    static long ma[56];
    static long iff = 0;
    long mj, mk;
    long i, ii, k;

    if (*idum < 0 || iff == 0) {
        iff = 1;
        mj = labs(MSEED - labs(*idum));
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;

        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }

        for (k = 1; k <= 4; k++)
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }

        inext = 0;
        inextp = 31;
        *idum = 1;
    }

    if (++inext == 56)
        inext = 1;
    if (++inextp == 56)
        inextp = 1;
    mj = ma[inext] - ma[inextp];

    if (mj < MZ)
        mj += MBIG;
    ma[inext] = mj;

    return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double NRC_gasdev(long *idum)
{
    static long iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (*idum < 0)
        iset = 0;

    if (iset == 0) {
        do {
            v1 = 2.0 * NRC_ran1(idum) - 1.0;
            v2 = 2.0 * NRC_ran1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;

        return v2 * fac;

    } else {
        iset = 0;
        return gset;
    }
}

void lux_Axsolution(long m, double **A0, double *g, double mu, double *h)
/* solve (A + mu * I)h = g, using SVD if necessary */
{
    long i, j;
    double **A, *p;  /* for Cholesky decomposition */

    A = dmatrix(1, m, 1, m);
    p = dvector(1, m);

    for (i = 1; i <= m; i++) {
        for (j = 1; j <= m; j++)
            A[i][j] = A0[i][j];
        A[i][i] += mu;
    }

    if (NRC_chol_posdef(A, m, p))
        NRC_cholsl(A, m, p, g, h);

    else {
        double **u, **v, *w;

        log_msg("\t\t*not* positive definite -- solved with svd");

        /* restore A to its original form */
        for (i = 2; i <= m; i++)
            for (j = 1; j < i; j++)
                A[i][j] = A[j][i];

        u = A;  /* for consistence */
        v = dmatrix(1, m, 1, m);
        w = dvector(1, m);

        NRC_svdcmp(u, m, m, w, v);
        NRC_svdzwi(w, m);
        NRC_svbksb(u, w, v, m, m, g, h);

        free_dmatrix(v, 1, DUMMY, 1, DUMMY);
        free_dvector(w, 1, DUMMY);
    }

    free_dmatrix(A, 1, DUMMY, 1, DUMMY);
    free_dvector(p, 1, DUMMY);
}

void get_covar(long m, double **A0, double **covar)
{
    long i, j;
    double **A, *p;  /* for Cholesky decomposition */

    A = dmatrix(1, m, 1, m);
    p = dvector(1, m);

    copy_dmatrix(A0, 1, m, 1, m, A);  /* to keep original values in A0 untouched */

    if (NRC_chol_posdef(A, m, p))
        NRC_chol_inv(A, m, p, covar);
    else {
        double **u, **v, *w;

        log_msg("\tget_covar(): not positive definite -- solved with svd");

        /* restoring A to its original form */
        for (i = 2; i <= m; i++)
            for (j = 1; j < i; j++)
                A[i][j] = A[j][i];

        u = A;  /* for consistence */
        v = dmatrix(1, m, 1, m);
        w = dvector(1, m);

        NRC_svdcmp(u, m, m, w, v);
        NRC_svdzwi(w, m);
        NRC_svd_inv(u, w, v, m, m, covar);

        free_dmatrix(v, 1, DUMMY, 1, DUMMY);
        free_dvector(w, 1, DUMMY);
    }

    free_dmatrix(A, 1, DUMMY, 1, DUMMY);
    free_dvector(p, 1, DUMMY);
}

static void NRC_rotate(double **a, long i, long j, long k, long l, double *g,
                       double *h, double s, double tau)
{
    *g = a[i][j];
    *h = a[k][l];
    a[i][j] = *g - s * (*h + *g * tau);
    a[k][l] = *h + s * (*g - *h * tau);
}

void NRC_eigsrt(double *d, double **v, long n)
/* sort eigenvalues into DESCENDING order, and rearrange eigenvectors */
{
    long i, j, k;
    double p;

    for (i = 1; i < n; i++) {
        p = d[k = i];
        for (j = i + 1; j <= n; j++)
            if (d[j] >= p)
                p = d[k = j];

        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; j++)
                dval_swap(&v[j][i], &v[j][k]);
        }
    }
}

void NRC_jacobi(double **a, long n, double *d, double **v)
{
    long i, j, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

    b = dvector(1, n);
    z = dvector(1, n);
    identity_dmatrix(v, n);
    for (ip = 1; ip <= n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }

    for (i = 1; i <= ITMAX; i++) {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++)
                sm += fabs(a[ip][iq]);
        }

        if (sm < TOLX) {
            free_dvector(z, 1, DUMMY);
            free_dvector(b, 1, DUMMY);
            NRC_eigsrt(d, v, n);
            return;
        }

        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;

        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip])
                    && (fabs(d[iq]) + g) == fabs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if ((fabs(h) + g) == fabs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++)
                        NRC_rotate(a, j, ip, j, iq, &g, &h, s, tau);
                    for (j = ip + 1; j <= iq - 1; j++)
                        NRC_rotate(a, ip, j, j, iq, &g, &h, s, tau);
                    for (j = iq + 1; j <= n; j++)
                        NRC_rotate(a, ip, j, iq, j, &g, &h, s, tau);
                    for (j = 1; j <= n; j++)
                        NRC_rotate(v, j, ip, j, iq, &g, &h, s, tau);
                }
            }
        }

        for (ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    fatal("too many iterations than ITMAX (%ld)\n", ITMAX);
}
