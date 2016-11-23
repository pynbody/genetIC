#include <mutex>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>

#include "progress.hpp"

namespace progress {

  void ProgressBar::remainingTime() {
    auto end = std::chrono::system_clock::now();
    double elapsed = ((end - start)).count() / std::chrono::system_clock::duration::period::den;
    double toEnd = elapsed * (1.0 - progress) / progress;

    if (progress > 0.01 && elapsed > 1.0 && toEnd > 1.0) {
      stream << int(toEnd) << "s remaining       ";
    } else {
      stream << "                  ";
    }
  }

  void ProgressBar::update() {
    stream << title << " ";
    stream.precision(3);
    stream.width(4);
    stream << (progress * 100);
    stream << "% |";

    int remaining_width = width - title.size() - 2;

    int bar_chars = remaining_width * progress;
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

    stream << "\r";
  }

  void ProgressBar::runUpdateLoop() {
    std::unique_lock<std::mutex> lock(mutex);
    start = std::chrono::system_clock::now();
    while (!terminated) {
      cv.wait_for(lock, std::chrono::milliseconds(100));
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
    nOpsCurrent+=1;
    setProgress(float(nOpsCurrent)/nOpsTotal);
  }

  ProgressBar::ProgressBar(const std::string &&title, size_t nOps) : title(title), stream(std::cerr), nOpsTotal(nOps) {
    width = 80;
    updateTimeout = 1000;
    progress = 0.0;
    terminated = false;
    thread = std::thread([this]() { this->runUpdateLoop(); });
  }

  ProgressBar::~ProgressBar() {
    terminate();
  }

  void ProgressBar::setProgress(float p) {
    progress = p;
  }

}

/*
int main() {
  ProgressBar pb("Hello, world");
  float prog = 0.0;
  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  while(prog<1.0) {
    prog+=0.001;
    pb.setProgress(prog);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

}
*/
