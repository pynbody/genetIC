#ifndef _DUMMYIC_HPP_INCLUDED
#define _DUMMYIC_HPP_INCLUDED

#include "ic.hpp"

template<typename MyFloat>
class DummyICGenerator : public ICGenerator<MyFloat> {
protected:
  ICGenerator<MyFloat> *pUnderlying;
public:
  DummyICGenerator(ICGenerator<MyFloat> *pUnderlying) : ICGenerator<MyFloat>(pUnderlying->interpreter), pUnderlying(pUnderlying) {

  }

  void zeroLevel(int level) override {

  }

  void applyPowerSpec() override { }

  void dumpGrid(int level) override { }

  void dumpPS(int level) override { }


  void write() override { }

  void constrain(string name, string type, float value) override { }

  void done() override { }
};

#endif
