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
        class ProgressBar {
        protected:
            std::string title;
            std::ostream &stream;
            std::thread thread;
            std::chrono::time_point<std::chrono::system_clock> start;

            std::condition_variable cv;
            std::mutex mutex;

            float progress;
            int width;
            bool terminated;
            int updateTimeout;

            size_t nOpsTotal;
            size_t nOpsCurrent;


            void update();

            void runUpdateLoop();

            void terminate();

            void remainingTime();

        public:
            ProgressBar(const std::string &&title, size_t nOps = 0);

            virtual ~ProgressBar();

            void setProgress(float p);

            void tick();

            size_t terminalWidth();

            bool isInteractive();

            void clearLine();
        };
    }
}
#endif
