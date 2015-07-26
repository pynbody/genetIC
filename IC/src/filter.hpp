//
// Created by Andrew Pontzen on 26/07/15.
//

#ifndef IC_FILTER_HPP
#define IC_FILTER_HPP

#include <math.h>
#include <stdexcept>


template<typename T>
class Filter {
public:
    virtual ~Filter() {}
    
    virtual T operator()(T x) const {
        return 1.0;
    }
};

template<typename T>
void tabulateFilter(Filter<T> *pF) {
    std::cerr << "tabulateFilter:" << std::endl;
    for(T k=0.1;k<3.0;k+=0.1) {
        std::cerr << " " << k << " " << (*pF)(k) << std::endl;
    }
}

template<typename T>
class LowPassFermiFilter : public Filter<T> {
private:
    T kcut;
    T temperature;

public:
    LowPassFermiFilter(T kcut_in) : kcut(kcut_in) {
        // The following arbitrary choice for the temperature (i.e.
        // sharpness of the cut-off) has been made by numerical experiments
        // trading off the various pros and cons of a sharp cut discussed
        // in the notes.
        temperature = kcut/10;
    };

    T operator()(T k) const override {
        return 1./(1.+exp((k-kcut)/temperature));
    }
};

template<typename T>
class CovarianceFilterAdaptor : public Filter<T> {
private:
    Filter<T> *pUnderlying;
public:
    CovarianceFilterAdaptor(Filter<T> *pOriginal) : pUnderlying(pOriginal) { };
    T operator()(T k) const override {
        T denFilter = (*pUnderlying)(k);
        return denFilter*denFilter;
    }
};

template<typename T>
class ComplementaryFilterAdaptor : public Filter<T> {
private:
    Filter<T> *pUnderlying;
public:
    ComplementaryFilterAdaptor(Filter<T> *pOriginal) : pUnderlying(pOriginal) { };
    T operator()(T k) const override {
        return 1.-(*pUnderlying)(k);
    }
};



template<typename T>
class FilterFamily {
    Filter<T> filter;
public:
    virtual ~FilterFamily() { };

    virtual Filter<T> & getFilterForDensityOnLevel(int level) {
        if(level==0)
            return filter;
        throw std::runtime_error("FilterFamily cannot produce filter for requested level");
    }
    virtual Filter<T> & getFilterForCovarianceOnLevel(int level) {
        if(level==0)
            return filter;
        throw std::runtime_error("FilterFamily cannot produce filter for requested level");
    }
    virtual Filter<T> & getFilterForCovarianceResidualOnLevel(int level) {
        throw std::runtime_error("FilterFamily cannot produce filter for requested level");
    }
};

template<typename T>
class TwoLevelFilterFamily : public FilterFamily<T> {
private:
    Filter<T> * pDenFilterLow;
    Filter<T> * pDenFilterHigh;
    Filter<T> * pCovFilterLow;
    Filter<T> * pCovFilterHighPartOfLowRes;
    Filter<T> * pCovFilterHigh;
public:

    TwoLevelFilterFamily(T k_cut) {
        pDenFilterLow = new LowPassFermiFilter<T>(k_cut);
        pDenFilterHigh = new ComplementaryFilterAdaptor<T>(pDenFilterLow);
        pCovFilterLow = new CovarianceFilterAdaptor<T>(pDenFilterLow);
        pCovFilterHigh = new ComplementaryFilterAdaptor<T>(pCovFilterLow);
        pCovFilterHighPartOfLowRes = new CovarianceFilterAdaptor<T>(pDenFilterHigh);
    }

    ~TwoLevelFilterFamily() {
        delete pDenFilterLow;
        delete pDenFilterHigh;
        delete pCovFilterLow;
        delete pCovFilterHighPartOfLowRes;
        delete pCovFilterHigh;
    }

    Filter<T> & getFilterForDensityOnLevel(int level) override {
        switch(level) {
            case 0:
                return *pDenFilterLow;
            case 1:
                return *pDenFilterHigh;
            default:
                throw std::runtime_error("TwoLevelFilterFamily cannot produce filter for requested level");
        }
    }

    Filter<T> & getFilterForCovarianceOnLevel(int level) override {
        switch(level) {
            case 0:
                return *pCovFilterLow;
            case 1:
                return *pCovFilterHigh;
            default:
                throw std::runtime_error("TwoLevelFilterFamily cannot produce filter for requested level");
        }
    }

    Filter<T> & getFilterForCovarianceResidualOnLevel(int level) override {
        switch(level) {
            case 0:
                return *pCovFilterHighPartOfLowRes;
            default:
                throw std::runtime_error("TwoLevelFilterFamily cannot produce filter for requested level");
        }
    }
};


#endif //IC_FILTER_HPP
