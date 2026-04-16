#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <TF1.h>
#include <TGraph.h>
#include "WF.h"

WF::WF() 
    : t(), V(), ch(0), timestamp("") {
}
WF::WF(const std::vector<double> &t, const std::vector<double> &V, int ch, const std::string &timestamp)
    : t(t) , V(V), ch(ch), timestamp(timestamp){}

double WF::Baseline(int N) const{

    double BL = 0;

    for(int i = 0; i < N; i++){ BL += V[i]; }

        BL /= N;

        return BL;

}
double WF::Amp(int N) const{

    double amp = 0;

    for(size_t i = 1; i < V.size(); i++){
        double new_amp = abs(V[i] - Baseline(N));
        if (new_amp > amp) {amp = new_amp;}

    }

        return amp;

}
double WF::DropTime(const double &w, int N) const{

    double thr = (Baseline(N) - Amp(N) * w);

    for (size_t i = 1; i < V.size(); i++) {

        if (V[i-1] > thr && V[i] <= thr) {

            double t1 = t[i-1];
            double t2 = t[i];

            double v1 = V[i-1];
            double v2 = V[i];

            double frac = (thr - v1) / (v2 - v1);

            return t1 + frac * (t2 - t1);
        }
    }

    return -1; 
}

double WF::RiseTime(int N, int M) const{

    size_t index = 0;

    double amp = 0.;

    for(size_t i = 1; i < V.size(); i++){

        double new_amp = abs(V[i] - Baseline(N));

        if (new_amp > amp) {

            amp = new_amp;

            index = i;

        }

    }

    double tau_sum = 0.0;

    for (int i = 1; i < M; i++){tau_sum += (t[index] - t[index - i])/(log(V[index]/V[index - i]));}

    tau_sum /= M;

    return tau_sum;
}

void WF::SetTimeStamp(const std::string &date){
    timestamp = date;
}

void WF::SetChannel(int n){
    ch = n;
}

void WF::PushBack(double t_val, double V_val) {

    t.push_back(t_val);

    V.push_back(V_val);

}

WF Event::GetChannel(int i) const {

    return waveforms[i];

}
