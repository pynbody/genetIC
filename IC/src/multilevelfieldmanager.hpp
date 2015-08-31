#ifndef _MULTILEVELFIELDMANAGER_HPP
#define _MULTILEVELFIELDMANAGER_HPP

#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <cassert>
#include <vector>

//#include <Eigen/Dense>


template<typename T>
class MultiLevelFieldManager {
private:
    std::vector<std::shared_ptr<Grid<T>>> pGrid;
    std::vector<std::vector<T>> C0s;
    std::vector<T> weights;
    std::shared_ptr<FilterFamily<T>> pFilters;
    std::vector<complex<T>> pField_k_0_high; // high-k modes for level 0

protected:
    std::vector<size_t> Ns;
    std::vector<size_t> cumu_Ns;
    size_t Ntot;
    size_t nLevels;
    bool enforceExactPowerSpectrum;

    void mapIdToLevelId(size_t i, size_t &level, size_t &level_id) {
        level =0;
        while(level < nLevels && i>=Ns[level]) {
            i-=Ns[level];
            level++;
        }
        if(i>=Ns[level])
            throw std::runtime_error("ID out of range when mapping into underlying fields in mapIdToLevelId");
        level_id = i;
    }

    MultiLevelFieldManager(size_t N) {
        nLevels =1;
        Ns.push_back(N);
        Ntot = N;
    }

    void setupFilters() {
        if(nLevels ==2) {
            const T FRACTIONAL_K_SPLIT = 0.3;
            pFilters = std::make_shared<TwoLevelFilterFamily<T>>(
                    ((T) pGrid[0]->size) * FRACTIONAL_K_SPLIT * 2. * M_PI / pGrid[0]->boxsize);
        } else if(nLevels ==1) {
            pFilters = std::make_shared<FilterFamily<T>>();
        } else {
            throw(std::runtime_error("Don't know how to make appropriate filters for multi-level field"));
        }
    }

public:
    MultiLevelFieldManager() {
        nLevels =0;
        Ntot=0;
        enforceExactPowerSpectrum = false;
    }

    virtual ~MultiLevelFieldManager() { }

    void addLevel(std::vector<T> C0, std::shared_ptr<Grid<T>> pG, size_t N) {
        C0s.push_back(std::move(C0));
        if(pGrid.size()==0) {
            weights.push_back(1.0);
        } else {
            weights.push_back(pow(pG->dx/pGrid[0]->dx,3.0));
        }
        pGrid.push_back(pG);
        Ns.push_back(N);
        cumu_Ns.push_back(Ntot);
        Ntot+=N;
        nLevels +=1;
        setupFilters();
    }

    void setExactPowerSpectrumEnforcement(bool value) {
        enforceExactPowerSpectrum = true;
    }

    size_t getNumCells() const {
        return Ntot;
    }

    size_t getNumLevels() const {
        return nLevels;
    }

    Grid<T> &getGridForLevel(size_t level) {
        return *pGrid[level];
    }

    std::vector<T> &getCovariance(size_t level) {
        return C0s[level];
    }


    virtual std::complex<T> cov(const std::vector<std::complex<T>> &vec, size_t i) {
        // returns Sum_j C[i,j] vec[j]
        // Here, though, C0 is diagonal
        size_t comp;
        size_t j;
        mapIdToLevelId(i, comp, j);
        return vec[i]*C0s[comp][j]*weights[comp];
    }

    virtual std::complex<T> cov(const std::vector<std::complex<T>> &vec, size_t comp, size_t j) {
        // returns Sum_i C[j,i] vec[i]
        // Here, though, C0 is diagonal
        size_t overall_j = j+cumu_Ns[comp];
        return C0s[comp][j]*vec[overall_j]*weights[comp];
    }

    virtual void splitLevel0() {
        pField_k_0_high = pGrid[0]->getFieldFourier();
        applySpectrumToField(pField_k_0_high,
                             C0s[0],
                             pFilters->getFilterForCovarianceResidualOnLevel(0),
                             *pGrid[0]);

    }

    virtual void recombineLevel0() {

        auto & pField_k = pGrid[0]->getFieldFourier();
        size_t size=this->pGrid[0]->size3;

        for(size_t i=0; i<size; i++)
            pField_k[i]+=pField_k_0_high[i];

        pField_k_0_high.clear();

    }

    vector<complex<T>> createEmptyFieldForLevel(size_t level) const {
        vector<complex<T>> ar(pGrid[level]->size3);
        return ar;
    }

    vector<complex<T>> combineDataOnLevelsIntoOneVector(const vector<vector<complex<T>>> &dataOnLevel) const {

        assert(dataOnLevel.size()==pGrid.size());

        for(size_t i=0; i<pGrid.size(); ++i) {
            assert(dataOnLevel[i].size()==pGrid[i]->size3);
        }

        vector<complex<T>> returnValue(getNumCells());

        size_t return_i=0;

        for(size_t level=0; level<pGrid.size(); ++level) {
            const size_t levelSize = dataOnLevel[level].size();
            const auto & dataOnThisLevel = dataOnLevel[level];
            for(size_t i=0; i<levelSize; ++i) {
                returnValue[return_i]=dataOnThisLevel[i];
                return_i++;
            }
        }

        return returnValue;
    }

    vector<vector<complex<T>>> generateMultilevelFromHighResField(vector<complex<T>> &&data) {
        assert(data.size()==pGrid.back()->size3);
        vector<vector<complex<T>>> dataOnLevels;
        for(size_t level=0; level<pGrid.size(); level++) {
            if(level==pGrid.size()-1) {
                dataOnLevels.emplace_back(std::move(data));
            } else {
                dataOnLevels.emplace_back(pGrid[level]->size3,0.0);
            }
        }
        size_t levelmax = dataOnLevels.size()-1;
        if(levelmax>0) {
            fft(dataOnLevels.back(), dataOnLevels.back(), -1);
            for(int level=levelmax-1; level>=0; --level) {

                pGrid[level]->addFieldFromDifferentGrid(*pGrid[levelmax], dataOnLevels[levelmax],
                                                        dataOnLevels[level]);


                fft(dataOnLevels[level], dataOnLevels[level], 1);
                pGrid[level]->applyFilter(pFilters->getFilterForDensityOnLevel(level),
                                          dataOnLevels[level]);
            }
            fft(dataOnLevels.back(), dataOnLevels.back(), 1);
            pGrid[levelmax]->applyFilter(pFilters->getFilterForDensityOnLevel(levelmax),
                                         dataOnLevels[levelmax]);
        }

        return dataOnLevels;

    }

    void applySpectrumToField(std::vector<std::complex<T>> &field,
                                const std::vector<T> &spectrum,
                                const Filter<T> &filter,
                                const Grid<T> &grid ) {
        T white_noise_norm = sqrt(T(grid.size3));
	
        #pragma omp parallel for
        for(size_t i=0; i<grid.size3; i++) {
            T k = grid.getAbsK(i);
            if(enforceExactPowerSpectrum) {
                T existing_norm = abs(field[i]);
                field[i] *= sqrt(spectrum[i] * filter(k)) * white_noise_norm/existing_norm;
            } else {
                field[i] *= sqrt(spectrum[i] * filter(k));
            }

        }
    }

    void applyCovarianceToWhiteNoiseField() {

        if(getNumLevels()>1)
            splitLevel0(); // TODO: this should probably be managed elsewhere

        for(size_t i=0; i<getNumLevels(); ++i) {
            applySpectrumToField(pGrid[i]->getFieldFourier(), C0s[i],
                                 pFilters->getFilterForCovarianceOnLevel(i),
                                 *pGrid[i]);
        }

    }

    void zeroLevel(int level) {
        if(level==-1) {
            std::fill(pField_k_0_high.begin(), pField_k_0_high.end(), 0);
        } else {
            auto &field = pGrid[level]->getField();
            std::fill(field.begin(), field.end(), 0);
        }
    }

    void forEachLevel(std::function<void(Grid<T> &)> newLevelCallback) {
        for(size_t level=0; level<nLevels; level++) {
            newLevelCallback(*pGrid[level]);
        }
    }

    void forEachCellOfEachLevel(
            std::function<void(size_t)> levelCallback,
            std::function<void(size_t, size_t, size_t, std::vector<std::complex<T>> &)> cellCallback,
            bool kspace=true)
    {


        for(size_t level =0; level < nLevels; level++) {

            levelCallback(level);
            size_t level_base = cumu_Ns[level];
            auto & field = kspace?pGrid[level]->getFieldFourier():pGrid[level]->getField();

            #pragma omp parallel for
            for(size_t i=0; i<Ns[level]; i++) {
                size_t i_all_levels = level_base + i;
                cellCallback(level, i, i_all_levels, field);
            }
        }

    }

    void forEachCellOfEachLevel(
            std::function<void(size_t, size_t, size_t, std::vector<std::complex<T>> &)> cellCallback,
            bool kspace=true) {

        forEachCellOfEachLevel([](size_t i){},cellCallback);
    }



    std::complex<T> accumulateOverEachCellOfEachLevel(
            std::function<void(size_t)> newLevelCallback,
            std::function<std::complex<T>(size_t, size_t, size_t)> getCellContribution)
    {

        T res_real(0), res_imag(0);


        for(size_t level =0; level < nLevels; level++) {

            newLevelCallback(level);
            size_t level_base = cumu_Ns[level];

            #pragma omp parallel for reduction(+:res_real,res_imag)
            for(size_t i=0; i<Ns[level]; i++) {
                size_t i_all_levels = level_base + i;
                std::complex<T> res = getCellContribution(level, i, i_all_levels);

                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }

        return std::complex<T>(res_real, res_imag);
    }

    std::complex<T> accumulateOverEachCellOfEachLevel(
            std::function<std::complex<T>(size_t, size_t, size_t, T)> getCellContribution) {
        T weight;
        auto newLevelCallback = [&weight,this](size_t comp) {
            weight = weights[comp];
        };
        auto getCellContributionWrapper = [&weight, &getCellContribution](size_t component, size_t i, size_t cumu_i) {
            return getCellContribution(component, i, cumu_i, weight);
        };
        return accumulateOverEachCellOfEachLevel(newLevelCallback, getCellContributionWrapper);
    }



    virtual std::complex<T> v1_dot_y(const std::vector<std::complex<T>> &v1, bool fourier=true){

        T weight;
        const std::vector<std::complex<T>> *pField;

        auto newLevelCallback = [&weight,&pField,this,fourier](size_t comp) {
            weight = weights[comp];
            if (fourier)
                pField = &pGrid[comp]->getFieldFourier();
            else
                pField = &pGrid[comp]->getFieldReal();
        };

        auto getCellContribution = [&weight, &pField,&v1,this](size_t comp, size_t i, size_t cumu_i) {
            return conj(v1[cumu_i])*(*pField)[i]*weight;
        };

        return accumulateOverEachCellOfEachLevel(newLevelCallback, getCellContribution);
    }

    virtual std::complex<T> v1_dot_v2(const std::vector<std::complex<T>> &v1,
                                      const std::vector<std::complex<T>> &v2){

        auto getCellContribution = [this,&v1,&v2](size_t comp, size_t i, size_t cumu_i, T weight) {
            return conj(v1[cumu_i])*v2[cumu_i]*weight;
        };

        return accumulateOverEachCellOfEachLevel(getCellContribution);

    }

    virtual std::complex<T> v1_cov_v2(const std::vector<std::complex<T>> &v1,
                                      const std::vector<std::complex<T>> &v2){

        auto getCellContribution = [this,&v1,&v2](size_t comp, size_t i, size_t cumu_i, T weight) {
            return conj(v1[cumu_i])*cov(v2,comp,i)*weight;
        };

        return accumulateOverEachCellOfEachLevel(getCellContribution);

    }

    virtual std::complex<T> v1_cov_v2_with_pseudo_crossterm(const std::vector<std::complex<T>> &v1,
                                                            const std::vector<std::complex<T>> &v2){

        T weight;
        Filter<T> *pFilt;
        auto newLevelCallback = [&weight,&pFilt,this](size_t comp) {
            weight = weights[comp];
            pFilt = &pFilters->getFilterForDensityOnLevel(comp);
        };
        auto getCellContribution = [&weight,&v1,&v2,&pFilt,this](size_t component, size_t i, size_t cumu_i) {
            T k=std::sqrt(pGrid[component]->getKSquared(i));
            return conj(v1[cumu_i])*cov(v2,component,i)*weight/(*pFilt)(k);
        };

        return accumulateOverEachCellOfEachLevel(newLevelCallback,getCellContribution);

    }

    virtual std::complex<T> v_cov_v(const std::vector<std::complex<T>> &v){
        return v1_cov_v2(v,v);
    }

    virtual std::complex<T> v_cov_v_with_pseudo_crossterm(const std::vector<std::complex<T>> &v){
        return v1_cov_v2_with_pseudo_crossterm(v,v);
    }


    virtual void get_realization(std::complex<T> *r) {
        // Get the realization and store it in the specified array
        size_t j=0;
        for(size_t comp=0; comp< nLevels; comp++) {
            auto field = pGrid[comp]->getFieldFourier();

            for(size_t i=0; i<Ns[comp]; i++) {
                r[j]=field[i];
                j++;
            }
        }
        assert(j==Ntot);
    }

    virtual void get_mean(std::complex<T> *r) {
        // Get the mean of all realizations and store it in the specified array
        for(size_t i=0; i<Ntot; i++) {
            if(i%1000==0) progress("get mean x", float(i)/Ntot);
            r[i]=0;
        }
        end_progress();
    }

    T get_field_chi2() {

        T chi2=0;

        for(size_t i=0; i<getNumLevels(); ++i) {

            auto & field = pGrid[i]->getFieldFourier();
            const auto & spectrum = C0s[i];
            T norm = T(pGrid[i]->size3);

            for(size_t j=0; j<field.size(); ++j) {
	      if (spectrum[j]!=0)
                chi2+=pow(abs(field[j]),2.0)/(spectrum[j]*norm);
            }
        }

        return chi2;

    }

};
#endif
