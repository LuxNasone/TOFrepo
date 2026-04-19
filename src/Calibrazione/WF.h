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

        double CFT(const double &w, int N) const;

        double ARC(const double &k, const double &dt) const;

        double maxDer() const;

        bool IsOsc(int N) const;

        bool IsClipped(double Vmax, double Vmin) const;

        bool IsBadWF(int N, double Vmax, double Vmin) const;

        void SetTimeStamp(const std::string &date);

        void SetChannel(int n);

        void PushBack(double t_val, double V_val);

};

class Event {

    public:

        std::vector<WF> waveforms; 

        Event() : waveforms(3) {} 

        const WF& GetChannel(int i) const;

        bool IsBad(int N, double Vmax, double Vmin) const;
};

#endif