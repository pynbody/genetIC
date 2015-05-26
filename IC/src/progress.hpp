#ifndef _PROGRESS_HPP_INCLUDED
#define _PROGRESS_HPP_INCLUDED

extern void progress(const std::string & message, size_t i, size_t len);
extern void progress(const std::string & message, float progress);
extern void end_progress();

#endif
