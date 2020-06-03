#include <mutex>
#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>

#include "progress.hpp"

namespace tools {
  namespace progress {

    //! Returns the size of the current terminal display
    size_t ProgressBar::terminalWidth() {
      struct winsize size;
      if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) < 0)
        return 0;
      else
        return size_t(size.ws_col);
    }

    //! Returns true if the terminal display is wider than 0
    bool ProgressBar::isInteractive() {
      return terminalWidth() > 0;
    }

    //! Returns an estimate of the time remaining to complete the operation
    void ProgressBar::remainingTime() {
      auto end = std::chrono::system_clock::now();
      double elapsed = ((end - start)).count() / std::chrono::system_clock::duration::period::den;
      double toEnd = elapsed * (1.0 - progress) / progress;

      if (progress > 0.01 && elapsed > 1.0 && toEnd > 1.0) {
        stream << int(toEnd) << "s remaining";
      }
    }

    //! Clears the current line - used to update the display
    void ProgressBar::clearLine() {
      size_t width = terminalWidth();
      stream << "\r";
      for (size_t i = 0; i < width; ++i)
        stream << " ";
      stream << "\r";
    }

    //! Update the display if possible
    void ProgressBar::update() {
      if (!isInteractive())
        return;

      clearLine();

      stream << title << " ";
      stream.precision(3);
      stream.width(4);
      stream << (progress * 100);
      stream << "% |";

      int remaining_width = int(terminalWidth());
      remaining_width -= title.size();
      remaining_width -= 30;

      if (remaining_width < 0) return;

      int bar_chars = remaining_width * progress;
      if (bar_chars > remaining_width)
        bar_chars = remaining_width;

      int i;
      for (i = 0; i < bar_chars - 1; i++) {
        stream << "-";
      }
      if (progress < 1.0)
        stream << ">";
      else
        stream << "-";
      i++;

      for (; i < remaining_width; i++) {
        stream << " ";
      }

      stream << "| ";

      remainingTime();

    }

    //! Updates the progress bar at the specified intervals, until it receives the signal to terminate
    void ProgressBar::runUpdateLoop() {
      std::unique_lock<std::mutex> lock(mutex);
      start = std::chrono::system_clock::now();
      while (!terminated) {
        cv.wait_for(lock, std::chrono::milliseconds(200));
        update();
      }
      update();
      stream << std::endl;
    }


    //! Signals the progress bar to terminate
    void ProgressBar::terminate() {
      setProgress(1.0);
      terminated = true;
      cv.notify_all();
      thread.join();
    }

    //! Updates the progress variable
    void ProgressBar::tick() {
      std::unique_lock<std::mutex> lock(mutex);
      nOpsCurrent += 1;
      setProgress(float(nOpsCurrent) / nOpsTotal);
    }

    //! Constructor with a given title and number of operations. Starts the progress bar as soon as it is created
    ProgressBar::ProgressBar(const std::string &&title, size_t nOps) : title(title), stream(std::cerr),
                                                                       nOpsTotal(nOps) {
      width = 80;
      updateTimeout = 1000;
      progress = 0.0;
      nOpsCurrent = 0;
      terminated = false;
      thread = std::thread([this]() { this->runUpdateLoop(); });
    }

    //! Destructor - terminates progress bar if still in progress.
    ProgressBar::~ProgressBar() {
      terminate();
    }

    //! Updates the progress variable
    void ProgressBar::setProgress(float p) {
      if (p > 1.0)
        p = 1.0;
      progress = p;
    }
  }
}

