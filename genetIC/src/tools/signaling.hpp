#ifndef _SIGNALING_HPP
#define _SIGNALING_HPP

#include "boost/signals2.hpp"

namespace tools {
/*! \class Signaling
    \brief Class to handle sending signals between classes, mostly to perform updates.

    Specifically, this is used by multi-level fields to signal the multi-level context
    to update itself when this becomes necessary.
*/
  class Signaling {
  public:

    typedef boost::signals2::signal<void()> signal_t;
    typedef boost::signals2::scoped_connection connection_t;

    //! Sets up a connection with the specified function to be called when signalling
    connection_t connect(const std::function<void()> &subscriber) const {
      return connection_t(signal.connect(subscriber));
    }

    //! Copy constructor
    Signaling(const Signaling & /* copy*/) {

    }

    //! Default constructor
    Signaling() {

    }

  protected:
    //! Calls the signal
    void changed() {
      signal();
    }

  private:
    mutable signal_t signal; //!< Signal to send when called.


  };
}
#endif
