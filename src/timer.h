//
//  timer.h
//
//  Created by Temesgen H. Dadi on 31/08/15.
//
//
#include <chrono>
#include <ctime>

#ifndef SEQAN_TIMER_H
#define SEQAN_TIMER_H

template<typename TimeT = std::chrono::seconds,
typename ClockT=std::chrono::high_resolution_clock,
typename DurationT=double>
class Timer
{
private:
    std::chrono::time_point<ClockT> _start, _end, _lapStart, _lapEnd;
public:
    Timer();
    DurationT lap();
    DurationT elapsed();
};

template<typename TimeT, typename ClockT, typename DurationT>
Timer<TimeT, ClockT, DurationT>::Timer()
{
    _start = _end = _lapStart = _lapEnd = ClockT::now();
}
template<typename TimeT, typename ClockT, typename DurationT>
DurationT Timer<TimeT, ClockT, DurationT>::lap()
{
    _lapEnd = ClockT::now();
    auto delta = std::chrono::duration_cast<TimeT>(_lapEnd-_lapStart);
    _lapStart = _lapEnd = ClockT::now();
    return delta.count();
}

template<typename TimeT, typename ClockT, typename DurationT>
DurationT Timer<TimeT, ClockT, DurationT>::elapsed()
{
    _end = ClockT::now();
    auto delta = std::chrono::duration_cast<TimeT>(_end-_start);
    return delta.count();
}

#endif //SEQAN_TIMER_H
