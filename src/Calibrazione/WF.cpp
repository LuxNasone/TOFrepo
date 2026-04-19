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
double WF::CFT(const double &w, int N) const{

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
double WF::ARC(const double &k, const double &dt) const {

    double time_step = t[1] - t[0];

    int step = static_cast<int>(std::round(dt / time_step));

    if (step <= 0 || step >= (int)V.size()){return -1;}

    std::vector<double> arc(V.size(), 0.0);

    for (size_t i = step; i < V.size(); i++){ arc[i] = V[i - step] - k * V[i]; }

    auto max_it = std::min_element(V.begin(), V.end());
    size_t min_idx = std::distance(V.begin(), max_it);


    for (size_t i = step; i < min_idx; ++i) {
        if (std::abs(arc[i]) < 1e-12){return t[i];}

        if (arc[i] * arc[i + 1] < 0) {

            double t0 = t[i];
            double t1 = t[i + 1];
            double y0 = arc[i];
            double y1 = arc[i + 1];

            return t0 - y0 * (t1 - t0) / (y1 - y0);

        }
    }

    return -1;

}
double WF::maxDer() const{
    
    double max_der = 0;

    size_t index = 0;

    double dt = t[1] - t[0];

    for(size_t i = 1; i < V.size(); i++){
        double new_der = std::abs((V[i + 1] - V[i - 1])/dt);
        if (new_der > max_der) {
            max_der = new_der;
            index = i;
        }
    }

    double drop_time = t[index];

    return drop_time;
    
}
bool WF::IsOsc(int N) const{

    double BL = Baseline(N);

    double A = Amp(N);

    int count = 0;
    double threshold = 0.25 * A;

    for (size_t i = 0; i < V.size(); i++) {
        if (std::abs(V[i] - BL) > threshold) {
            count++;
        }
    }

    return count > 0.1 * V.size();
}
bool WF::IsClipped(double Vmax, double Vmin) const {

    const double eps = 1e-3 * (Vmax - Vmin); 

    int countHigh = 0;
    int countLow = 0;

    int maxRun = 0;
    int currentRun = 0;

    for (size_t i = 0; i < V.size(); i++) {

        if (V[i] >= Vmax - eps) {
            countHigh++;
            currentRun++;
        }
        else if (V[i] <= Vmin + eps) {
            countLow++;
            currentRun++;
        }
        else {
            currentRun = 0;
        }

        if (currentRun > maxRun)
            maxRun = currentRun;
    }

    bool enoughSaturation = (countHigh + countLow) > 0.05 * V.size();
    bool longPlateau = maxRun > 5; 

    return enoughSaturation || longPlateau;
}
bool WF::IsBadWF(int N, double Vmax, double Vmin) const{
    return IsOsc(N) || IsClipped(Vmax,Vmin);
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
const WF& Event::GetChannel(int i) const{

    return waveforms[i];

}
bool Event::IsBad(int N, double Vmax, double Vmin) const{

    bool bad_1 = this->GetChannel(0).IsBadWF(N, Vmax, Vmin);

    bool bad_2 = this->GetChannel(1).IsBadWF(N, Vmax, Vmin);

    bool bad_3 = this->GetChannel(2).IsBadWF(N, Vmax, Vmin);

    return bad_1 || bad_2 || bad_3;

}
