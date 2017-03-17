#ifndef _SIGNALING_HPP
#define _SIGNALING_HPP

#include "src/tools/boost/signals2.hpp"

class Signaling {
public:

  typedef boost::signals2::signal<void()> signal_t;
  typedef boost::signals2::scoped_connection connection_t;

  connection_t connect(const std::function<void()> &subscriber) {
    return connection_t(signal.connect(subscriber));
  }


protected:
  void changed() {
    signal();
  }

private:
  signal_t signal;


};

#endif