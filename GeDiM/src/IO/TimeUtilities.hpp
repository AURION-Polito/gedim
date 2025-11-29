// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __TimeUtilities_H
#define __TimeUtilities_H

#include "IOUtilities.hpp"
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
        Gedim::Output::Assert(_end != timep_t{} && "Timer must toc before reading the time");
        return std::chrono::duration_cast<duration_t>(_end - _start);
    }
};

} // namespace Gedim

#endif
