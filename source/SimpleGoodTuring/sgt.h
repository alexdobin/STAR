#ifndef SGT_H
#define SGT_H
// Simple Good-Turing estimation
//
// Copyright (c) David Elworthy 2004.

// A class for implementing simple Good-Turing re-estimation, as described by
// Geoff Sampson in the book Empirical Linguistics (2001), and on the web at
// http://www.grsampson.net/RGoodTur.html. The code here is a C++ adaptation
// of the published code by Sampson and Gale, with the bug fix issued in
// 2000. It is encapsulated as a class to allow it to be incorporated into
// other programs. An additional coding change is that the data can be
// presented in any order, whereas the original code required the data to be
// in ascending order.
//
// Copyright (c) David Elworthy 2004.
// All rights reserved.
//
// Redistribution and use in source and binary forms for any purpose, with or
// without modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions, and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions, and the disclaimer that follows 
//    these conditions in the documentation and/or other materials 
//    provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// You may contact me at david@friendlymoose.com.

#include <map>
#include <vector>
#include <cmath>
using namespace std;

// Simple Good-Turing class.
// To use the class, create an SGT object and data to it by calling add() with
// each data point. A data point consists of the observed value and the
// frequency of the the observation (what Sampson and Gale refer to as the
// frequency, and the frequency of the the frequency). When you have added all
// the data points call analyse(). You can then call estimate() with an
// observed value as argument to get the estimated frequency for that value,
// or call iterate() to iterate over the data points. There is one special
// case to estimate(). If called with an argument of zero, it delivers the
// estimated frequency for unseen events. This is not delivered from pair().
// To get back from the estimate value to the smoothed value of the
// observation, multiply by total();
//
// In the original Sampson and Gale version, the observation was an integer.
// For this version, we make the code be a template over the observation type.
// However, it must always be some suitable numeric type, such as int or double.
//
// The code is implemented using the Standard Template Library (STL).

template <class ObsType> class SGT
{
    private:
    // Data block, holding the frequency and estimate. The estimate is set up
    // by analyse().
    struct Data
    {
        Data(unsigned int f) : freq(f), estimate(0) {}
        
        unsigned int freq;
        double estimate;
    };
    
    // Internal representation, as a map from observations to frequencies.
    // After calling analyse(), it provides the estimates as well.
    typedef map<ObsType, Data, less<ObsType> > DataMap;

    // Minimum number of data points for a valid analysis
#ifdef _WIN32
#define MinInput (5)
#else
    static const unsigned int MinInput = 5;
#endif
    
    template <class T> double sq(T d) { return ((double) d)*d; }

    double smoothed(ObsType i, double intercept, double slope)
    { return (exp(intercept + slope * log((double) i))); }
    
    public:
    // Iterator type for iterate();
    typedef typename DataMap::const_iterator iterator;
    
    // Construct a SGT object.
    SGT() : totalObs(0) {}

    // Destroy SGT object.
    ~SGT() {}

    // Add a data point.
    // If an observation with the same value has already been supplied, this adds
    // to its frequency.
    void add(ObsType observation, unsigned int frequency)
    {
        typename DataMap::iterator i = data.find(observation);
        if (i == data.end())
            data.insert(make_pair(observation, Data(frequency)));
        else
            (*i).second.freq += frequency;

        totalObs += observation * frequency;
    }

    // Get total number of observations (= sum of obs*freq)
    ObsType total() const { return totalObs; }

    // Analyse the data.
    // Returns false if there is not enough data for a valid analysis.
    // In this case, the estimate is set to the original value.
    bool analyse()
    {
        if (data.size() < MinInput)
            return false;

        // The code which follows is based on S and G's analyseInput()
        ObsType bigN = 0;
        unsigned int rows = data.size();

        // j could be declared in each for statement, but has to be here for
        // Visual C++, which disobeys the ANSI standard on variable scope.
        typename DataMap::iterator j;
        for (j = data.begin(); j != data.end(); ++j)
            bigN += (*j).first * (*j).second.freq;
    
        // Find the frequency for observation of value 1, if any
        iterator row1 = row(1, data.begin());
        PZero = (row1 == data.end()) ? 0 : (*row1).second.freq / (double) bigN;

        // Set up internal arrays
        vector<double> log_obs(rows);
        vector<double> log_Z(rows);
        vector<double> rStar(rows);

        double XYs = 0, Xsquares = 0, meanX = 0, meanY = 0;
        ObsType prevObs = 0;
        unsigned int r = 0;
    
        for (j = data.begin(); j != data.end(); ++r)
        {
            ObsType obs = (*j).first;
            Data &d = (*j).second;
        
            double k = (++j == data.end())
                ? (double) (2 * obs - prevObs) : (double) (*j).first;

            double Z   = 2 * d.freq / (k - prevObs);
            log_obs[r] = log((double) obs);
            log_Z[r]   = log(Z);

            meanX += log_obs[r];
            meanY += log_Z[r];

            prevObs = obs;
        }

        // Find the line with the best fit.
        meanX /= rows;
        meanY /= rows;

        for (r = 0; r < rows; ++r)
        {
            XYs += (log_obs[r] - meanX) * (log_Z[r] - meanY);
            Xsquares += sq(log_obs[r] - meanX);
        }
        double slope = XYs / Xsquares;
        double intercept = meanY - slope * meanX;
    
        // Now construct the estimates smoothing using the fitted line.
        bool indiffValsSeen = false;
        
        for (j = data.begin(), r = 0; j != data.end(); ++j, ++r)
        {
            ObsType obs = (*j).first;
            Data &d = (*j).second;

            ObsType obs1 = obs + 1;
            double y = obs1 * smoothed(obs1, intercept, slope)
                / smoothed(obs, intercept, slope);

            iterator nextRow = row(obs1, j);
            if (nextRow == data.end())
            {
                indiffValsSeen = true;
            }
            else if (!indiffValsSeen)
            {
                unsigned int next_n = (*nextRow).second.freq;
                unsigned int freq   = d.freq;

                double x = obs1 * next_n / (double) freq;
                if (fabs(x - y) <= 1.96 * sqrt(sq(obs1) * next_n
                        / (sq(freq)) * (1 + next_n / (double) freq)))
                {
                    indiffValsSeen = true;
                }
                else
                {
                    rStar[r] = x;
                }
            }
            
            if (indiffValsSeen)
            {
                rStar[r] = y;
            }
        }
    
        double bigNprime = 0.0;
        for (j = data.begin(), r = 0; j != data.end(); ++j, ++r)
            bigNprime += (*j).second.freq * rStar[r];

        for (j = data.begin(), r = 0; j != data.end(); ++j, ++r)
            (*j).second.estimate = (1 - PZero) * rStar[r] / bigNprime;
        
        return true;
    }

    // Analyze the data.
    // This just calls analyse(), and is included as a concession to speakers
    // of debased dialects of English.
    void analyze() { analyse(); }

    // Get the estimate for an observation.
    // If there was no such observation, return false.
    // Otherwise return true and yield the estimate.
    bool estimate(ObsType observation, double &estimate) const
    {
        if (observation == 0)
        {
            estimate = PZero;
            return true;
        }
        
        iterator rownum = row(observation, data.begin());
        if (rownum == data.end())
        {
            return false;
        }

        estimate = (*rownum).second.estimate;
        return true;
    }
    
    // Get start and end iterators over the data map.
    // You do not derefence these iterators directly, but instead used the
    // access functions, obs, freq and estimate.
    pair<iterator, iterator> iterate() const
    { return make_pair(data.begin(), data.end()); }
    
    // Get the observation from an iterator.
    ObsType obs(iterator i) const { return (*i).first; }

    // Get the frequency from an iterator (as supplied by add).
    unsigned int freq(iterator i) const { return (*i).second.freq; }

    // Get the estimated relative frequency from an iterator.
    double estimate(iterator i) const { return (*i).second.estimate; }

    private:
    // The data points
    DataMap data;

    // Zero estimate (only valid after a call to analyse()).
    double PZero;

    // Total number of observations
    ObsType totalObs;

    // Find the last row of the data which has a value equals to obs.
    // If there is no such value, return data.end().
    // start is a hint about where to start searching.
    iterator row(ObsType obs, iterator start) const
    {
        iterator j = start;
        
        while (j != data.end() && (*j).first < obs)
            ++j;

        return ((j != data.end() && (*j).first == obs) ? j : data.end());
    }
};

#endif //SGT_H
