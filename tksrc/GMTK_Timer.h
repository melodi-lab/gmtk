/*-
 * GMTK_Timer.h
 *     Polling based timer class
 *
 * Written by Chris Bartels 
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
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

