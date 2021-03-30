#include <bits/stdc++.h>
#include <math.h>

using namespace std;

const double cf = 2;
const double x0 = 0;
const double v0 = 700;
/// a = 10e-5 b = 1 c = 10e5
double U_x(double x){
    return 2.39520 * x * cosh(x*x);
}

void Euler(int N, double h){
    vector <double> t(N + 1);
    vector <double> v(N + 1);
    vector <double> x(N + 1);
    t[0] = 0; x[0] = x0, v[0] = v0;
    for(int i = 1; i <= N; i++){
        t[i] = i * h;
        v[i] = v[i - 1] + h * U_x(t[i]);
        x[i] = x[i - 1] + h * v[i - 1];
    }

    ofstream fout("euler_d");
    for(int i = 0; i <= N; i++){
        fout << x[i] << " " << t[i] << "\n";
    }
    fout.close();
}

void Adams(int N, double h){
    vector <double> t(N + 1);
    vector <double> v(N + 1);
    vector <double> x(N + 1);
    t[0] = 0;
    v[0] = v0;
    x[0] = x0;
    t[1] = h;
    v[1] = v[0] + h * U_x(t[0]);
    x[1] = x[0] + h * v[t[0]];
    for(int i = 2; i <= N; i++){
        t[i] = i * h;
        v[i] = v[i - 1] + (3.0 * U_x(t[i]) -  U_x(t[i - 1])) * h / 2.;
        x[i] = x[i - 1] + (3.0 * v[i] - v[i - 1]) * h / 2.;
    }
    ofstream fout("adams_d");
    for(int i = 0; i <= N; i++){
        fout << x[i] << " " << t[i] << "\n";
    }
    fout.close();

}
void RungeKutte(int N, double h){
    vector <double> t(N + 1);
    vector <double> v(N + 1);
    vector <double> x(N + 1);
    t[0] = 0;
    v[0] = v0;
    x[0] = x0;
    vector<double> tv(4), tx(4);
    for(int i = 1; i <= N; i++){
        t[i] = i * h;

        tv[0] = U_x(t[i]);
        tv[1] = U_x(t[i]) + tv[0] * h / 2.;
        tv[2] = U_x(t[i]) + tv[1] * h / 2.;
        tv[3] = U_x(t[i]) + tv[2] * h;
        v[i] = v[i - 1] + h / 6. * (tv[0] + 2 * tv[1] + 2 * tv[2] + tv[3]);

        tx[0] = v[i];
        tx[1] = v[i] + tx[0] * h / 2.;
        tx[2] = v[i] + tx[2] * h / 2.;
        tx[3] = v[i] + tx[3] * h;
        x[i] = x[i - 1] + h / 6. *(tx[0] + 2 * tx[1] + 2 * tx[2] + tx[3]);
    }
}
/// u'' + p(x)*u' + q(x)*y = f(x), u = u(x)
/// a_0 * u(a) + a_1 * u'(a) = A, b_0*u(b) + b_1 * u'(b) = B
/// a_0 = 0, a_1 = 1, b_0 = 1, b_1 = 0 p(x) = -7, q(x) = 12, f(x) = 5, A = 1, B = 2;
/// u'' - 7 * u' + 12 * y = 5, u = u(x)
/// u'(a) = 1, u(b) = 2


void triDiag(vector <double>& a, vector <double>& b, vector <double>& c, vector <double>& d, vector <double>& y){
    int n = c.size() - 1;
    vector <double> eps(n + 2);
    vector <double> nu(n + 2);
    eps[0] = -c[0] / b[0];
    nu[0] = d[0] / b[0];
    for(int i = 1; i <= n + 1; i++){
        eps[i] = - c[i]/(b[i] + a[i] * eps[i - 1]);
        nu[i] = (d[i] - a[i] * nu[i - 1])/(b[i] + a[i] * eps[i - 1]);
    }
    y[n] = nu[n + 1];
    for(int i = n  ; i >= 0; i--){
        y[i] = eps[i] * y[i + 1] + nu[i];
    }

}
void mkch(int N, double h){
    vector <double> a(N + 1), b(N + 1), c(N + 1), d(N + 1), u(n + 1);

    /*..*/
    a[0] = 0;
    b[0] = -1;
    c[0] = 1;
    d[0] = h;
    a[N] = 0;
    b[N] = 1;
    c[N] = 0;
    d[N] = 2;

    for(int i = 1 ; i < N; i++){
        a[i] = 1 - (-7) * h / 2.
        b[i] = -2 + h*h * 12;
        c[i] = 1 + 7 * h / 2;
        d[i] = h*h * 5;
    }

}



int main()
{

    return 0;
}
