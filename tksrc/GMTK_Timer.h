/*-
 * GMTK_Timer.h
 *     Polling based timer class
 *
 * Written by Chris Bartels 
 *
 * Copyright (c) 2003, < fill in later >
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle, and Jeff Bilmes make no representations about
 * the suitability of this software for any purpose.  It is provided
 * "as is" without express or implied warranty.
 *
 */

#ifndef GMTK_TIMER_H
#define GMTK_TIMER_H

#include <string>
#include <time.h>

class TimerClass {

public:

    TimerClass(); 
    TimerClass( time_t seconds ); 

    void   DisableTimer();
    void   EnableTimer();
    bool   Expired();
    void   Reset( time_t seconds );
    time_t SecondsLeft();
    time_t parseTimeString( std::string timeString );

private:

    bool   enabled;
    time_t expiration_time;
};

#endif

