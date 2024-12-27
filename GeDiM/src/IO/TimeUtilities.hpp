#ifndef __TimeUtilities_H
#define __TimeUtilities_H

#include <iostream>
#include <limits>
#include <vector>

#include "chrono"

namespace Gedim
{

template <class DT = std::chrono::milliseconds, class ClockT = std::chrono::steady_clock> class Timer
{
    using timep_t = decltype(ClockT::now());

    timep_t _start = ClockT::now();
    timep_t _end = {};

  public:
    void tick()
    {
        _end = timep_t{};
        _start = ClockT::now();
    }

    void tock()
    {
        _end = ClockT::now();
    }

    template <class duration_t = DT> duration_t duration() const
    {
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time");
        return std::chrono::duration_cast<duration_t>(_end - _start);
    }
};

} // namespace Gedim

#endif
