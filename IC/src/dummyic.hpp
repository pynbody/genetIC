#ifndef _DUMMYIC_HPP_INCLUDED
#define _DUMMYIC_HPP_INCLUDED

#include "ic.hpp"

template<typename MyFloat>
class DummyIC  : public IC<MyFloat> {
protected:
  IC<MyFloat> *pUnderlying;
public:
  DummyIC(IC<MyFloat> *pUnderlying) : IC<MyFloat>(pUnderlying->interpreter), pUnderlying(pUnderlying) {

  }

  void initGrid(unsigned int level) override {

      if(this->n[level]<0 || this->boxlen[level]<0)
          return;

      this->nPartLevel[level] = ((long)this->n[level]*this->n[level])*this->n[level];
      this->dx[level] = this->boxlen[level]/this->n[level];

      if(this->pGrid.size()!=level)
          throw std::runtime_error("Trying to re-initialize a grid level");

      if(pUnderlying->n[level]!=this->n[level])
        throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid n)");

      if(pUnderlying->dx[level]!=this->dx[level])
        throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid dx)");

      if(pUnderlying->x_off[level]!=this->x_off[level] || pUnderlying->y_off[level]!=this->y_off[level] || pUnderlying->z_off[level]!=this->z_off[level])
        throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid origin)");

      this->pGrid.push_back(this->pUnderlying->pGrid[level]);

      this->initMapper();

  }

  void drawRandom() override {

  }

  void zeroLevel(int level) override {

  }

  void interpolateIntoLevel(int level) override {

  }

  void splitLevel0() override { }


  void recombineLevel0() override { }

  void applyPowerSpecForLevel(int level, bool high_k) override { }

  void applyPowerSpec() override { }

  void dumpGrid(int level) override { }

  void dumpPS(int level) override { }

  void zeldovichForLevel(int level) override { }

  void zeldovich() override { }

  void write() override { }

  void prepare() override { }

  void constrain(string name, string type, float value) override { }

  void done() override { }
};

#endif
