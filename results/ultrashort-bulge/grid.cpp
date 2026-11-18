// Profile an independent source flux and baseline per stream for each template.
// Inputs are time-sorted. Outside |(t-t0)/tE|=32 the excess magnification is set to zero.
#include <algorithm>
#include <cmath>
#include <vector>

extern "C" void evaluate(int n, const double *time, const double *flux, const double *weight,
                         const int *group, int ng, const double *source_cap,
                         int nt, const double *centers, int nd, const double *scales,
                         int nm, int half, const double *bank, double *out) {
    const double factor = half / std::log1p(32.);
    const int width = 2*half + 1;
    std::vector<double> sum_w(ng), mean(ng), x(ng), xx(ng), xy(ng);
    for (int i=0;i<n;++i) {
        sum_w[group[i]] += weight[i];
        mean[group[i]] += weight[i]*flux[i];
    }
    for (int g=0;g<ng;++g) mean[g] /= sum_w[g];
    std::vector<int> index(n);
    std::vector<double> fraction(n), centered(n);
    for (int i=0;i<n;++i) centered[i] = weight[i]*(flux[i]-mean[group[i]]);
    for (int a=0;a<nt;++a) for (int b=0;b<nd;++b) {
        double t0=centers[a], te=scales[b];
        int left=std::lower_bound(time,time+n,t0-32*te)-time;
        int right=std::upper_bound(time,time+n,t0+32*te)-time;
        for (int i=left;i<right;++i) {
            double tau=(time[i]-t0)/te;
            double k=half+std::copysign(factor*std::log1p(std::abs(tau)),tau);
            index[i]=std::clamp(int(k),0,width-2);
            fraction[i]=std::clamp(k-index[i],0.,1.);
        }
        for (int m=0;m<nm;++m) {
            std::fill(x.begin(),x.end(),0.);
            std::fill(xx.begin(),xx.end(),0.);
            std::fill(xy.begin(),xy.end(),0.);
            const double *model=bank+m*width;
            for (int i=left;i<right;++i) {
                const int j=index[i], g=group[i];
                double amp=model[j]+fraction[i]*(model[j+1]-model[j]);
                x[g]+=weight[i]*amp;
                xx[g]+=weight[i]*amp*amp;
                xy[g]+=centered[i]*amp;
            }
            double gain=0.;
            for (int g=0;g<ng;++g) {
                double var=xx[g]-x[g]*x[g]/sum_w[g];
                if (var<=0.) continue;
                double fs=std::clamp(xy[g]/var,0.,source_cap[g]);
                gain+=fs*(2*xy[g]-fs*var);
            }
            out[(a*nd+b)*nm+m]=gain;
        }
    }
}
