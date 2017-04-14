#include <mutex>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>
#include <sys/ioctl.h>
#include <zconf.h>

#include "progress.hpp"

namespace tools {
  namespace progress {

    size_t ProgressBar::terminalWidth() {
      struct winsize size;
      if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) < 0)
        return 0;
      else
        return size_t(size.ws_col);
    }

    bool ProgressBar::isInteractive() {
      return terminalWidth() > 0;
    }

    void ProgressBar::remainingTime() {
      auto end = std::chrono::system_clock::now();
      double elapsed = ((end - start)).count() / std::chrono::system_clock::duration::period::den;
      double toEnd = elapsed * (1.0 - progress) / progress;

      if (progress > 0.01 && elapsed > 1.0 && toEnd > 1.0) {
        stream << int(toEnd) << "s remaining";
      }
    }

    void ProgressBar::clearLine() {
      size_t width = terminalWidth();
      stream << "\r";
      for (size_t i = 0; i < width; ++i)
        stream << " ";
      stream << "\r";
    }

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


    void ProgressBar::terminate() {
      setProgress(1.0);
      terminated = true;
      cv.notify_all();
      thread.join();
    }

    void ProgressBar::tick() {
      std::unique_lock<std::mutex> lock(mutex);
      nOpsCurrent += 1;
      setProgress(float(nOpsCurrent) / nOpsTotal);
    }

    ProgressBar::ProgressBar(const std::string &&title, size_t nOps) : title(title), stream(std::cerr),
                                                                       nOpsTotal(nOps) {
      width = 80;
      updateTimeout = 1000;
      progress = 0.0;
      nOpsCurrent = 0;
      terminated = false;
      thread = std::thread([this]() { this->runUpdateLoop(); });
    }

    ProgressBar::~ProgressBar() {
      terminate();
    }

    void ProgressBar::setProgress(float p) {
      if (p > 1.0)
        p = 1.0;
      progress = p;
    }

  }
}