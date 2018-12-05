#ifndef _PROGRESS_HPP
#define _PROGRESS_HPP

#include <mutex>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>
#include <condition_variable>

/*!
    \namespace tools::progress
    \brief Provides a tool to create a progress bar

 */
namespace tools {
  namespace progress {
    /*! \class ProgressBar
        \brief Class detailing output of a progress bar to the terminal.
    */
    class ProgressBar {
    protected:
      std::string title; //!< Title of the progress bar
      std::ostream &stream; //!< stream to output to.
      std::thread thread; //!< Current thread being used to operate the progress bar
      std::chrono::time_point<std::chrono::system_clock> start; //!< Time at which the progress bar started

      std::condition_variable cv; //!< Condition variable, used to wait between updates, and notify of progress bar termination
      std::mutex mutex;

      float progress; //!< current progress
      int width; //!< Doesn't appear to be used
      bool terminated; //!< If true, the progress bar will terminate instead of updating
      int updateTimeout; //!< Time to wait between updates before timing out

      size_t nOpsTotal; //!< Total number of operations that have to be performed.
      size_t nOpsCurrent; //!< Total number of operations that have been performed.

      //! Update the display if possible
      void update();

      //! Updates the progress bar at the specified intervals, until it receives the signal to terminate
      void runUpdateLoop();

      //! Signals the progress bar to terminate
      void terminate();

      //! Returns an estimate of the time remaining to complete the operation
      void remainingTime();

    public:
      //! Constructor with a given title and number of operations. Starts the progress bar as soon as it is created
      ProgressBar(const std::string &&title, size_t nOps = 0);

      //! Destructor - terminates progress bar if still in progress.
      virtual ~ProgressBar();

      //! Updates the progress variable
      void setProgress(float p);

      //! Updates the progress variable
      void tick();


      //! Returns the size of the current terminal display
      size_t terminalWidth();

      //! Returns true if the terminal display is wider than 0
      bool isInteractive();

      //! Clears the current line - used to update the display
      void clearLine();
    };
  }
}
#endif
