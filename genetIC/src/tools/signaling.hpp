#ifndef _SIGNALING_HPP
#define _SIGNALING_HPP

#include "boost/signals2.hpp"

namespace tools {
  class Signaling {
  public:

    typedef boost::signals2::signal<void()> signal_t;
    typedef boost::signals2::scoped_connection connection_t;

    connection_t connect(const std::function<void()> &subscriber) {
      return connection_t(signal.connect(subscriber));
    }

    Signaling(const Signaling& /* copy*/) {

    }

    Signaling () {

    }

  protected:
    void changed() {
      signal();
    }

  private:
    signal_t signal;


  };
}
#endif