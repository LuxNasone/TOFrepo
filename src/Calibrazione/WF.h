#ifndef WF_H
#define WF_H

class WF{

    private:

        std::vector<double> t;
        std::vector<double> V;
        int ch;
        std::string timestamp;

    public:

        WF();

        WF(const std::vector<double>& t, const std::vector<double>& V, int ch, const std::string& timestamp);

        double Baseline(int N) const;

        double Amp(int N) const;

        double DropTime(const double &w, int N) const;

        double RiseTime(int N, int M) const;

        void SetTimeStamp(const std::string &date);

        void SetChannel(int n);

        void PushBack(double t_val, double V_val);

};

class Event {

    public:

        std::vector<WF> waveforms; 

        Event() : waveforms(3) {} 

        WF GetChannel(int i) const;
};

#endif