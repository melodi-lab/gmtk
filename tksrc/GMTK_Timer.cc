/*-
 * GMTK_Timer.cc
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


#include <string>
#include <vector>

#include <ctype.h> 

#include "debug.h"
#include "error.h"

#include "GMTK_Signals.h"
#include "GMTK_Timer.h"

using namespace std;

/*-
 *-----------------------------------------------------------------------
 * TimerClass::TimerClass
 *   default constructor
 *
 * Preconditions:
 *   none
 *
 * Postconditions:
 *   timer is disabled 
 *               
 * Side Effects:
 *   none
 *   
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
TimerClass::
TimerClass()
{
  Reset(0);
  DisableTimer();
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::TimerClass
 *   constructor initializing the timer with an integer
 *         
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   timer is enabled 
 *               
 * Side Effects:
 *   none
 *   
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
TimerClass::
TimerClass(
  time_t seconds
  )
{
  Reset(seconds);
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::DisableTimer
 *
 * Preconditions:
 *   timer is disabled
 *   
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *-----------------------------------------------------------------------
 */
void TimerClass::
DisableTimer()
{
  enabled = false; 
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::EnableTimer
 *      
 * Preconditions:
 *   timer is enabled
 *   
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void TimerClass::
EnableTimer()
{
  enabled = true; 
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::Expired
 *      
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   none   
 *               
 * Side Effects:
 *   none
 *
 * Results:
 *   true if timer has expired, false if it has not 
 *-----------------------------------------------------------------------
 */
bool TimerClass::
Expired()
{
  return( SecondsLeft() < 1 ); 
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::Reset
 *      
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   timer is reset to a specified number of seconds 
 *               
 * Side Effects:
 *   none
 *   
 * Results:
 *   none
 *-----------------------------------------------------------------------
 */
void TimerClass::
Reset(
  time_t seconds
  )
{
  time_t crrnt_time;

  crrnt_time = time(NULL); 
  expiration_time = crrnt_time + seconds; 
  EnableTimer();
}

/*-
 *-----------------------------------------------------------------------
 * TimerClass::SecondsLeft
 *      
 * Preconditions:
 *   none
 *   
 * Postconditions:
 *   none
 *               
 * Side Effects:
 *   none
 *   
 * Results:
 *   The number of seconds left on the timer.  Seconds will be negative
 *   if timer has expired.  If the timer is disabled it will return 1,
 *   if a terminate signal has been received it will return -1. 
 *-----------------------------------------------------------------------
 */
time_t TimerClass::
SecondsLeft()
{
  time_t crrnt_time;
  time_t seconds_left;
  
  if (TerminateSignalReceived()) {
    seconds_left = -1;
  }
  else if (! enabled) {
    seconds_left = 1; 
  }
  else {
    crrnt_time = time(NULL);
    seconds_left = ( expiration_time - crrnt_time ); 
  }

  return( seconds_left ); 
}


/*-
 *-----------------------------------------------------------------------
 * TimerClass::ParseTimeString 
 *   Parses a time string and converts the expression to seconds.  The
 *   strings it recognizes are seconds, minutes, hours, days, and weeks.
 *   Partial words are accepted and  whitespace and non-alphanumeric 
 *   characters are ignored.
 *
 * Preconditions:
 *   none
 * 
 * Postconditions:
 *   none
 *
 * Side Effects:
 *   none
 *
 * Results:
 *   The number of seconds represented by the given string  
 *-----------------------------------------------------------------------
 */
time_t TimerClass::
parseTimeString(
  string timeString
  )
{
  vector<string> values;
  vector<string> types;
  string         new_string;
  unsigned       crrnt_string;
  time_t         total_seconds;
  bool           reading_digits = true;

  //////////////////////////////////////////////////////////////////////
  // Parse the values and labels from the string into vectors values 
  // and types
  //////////////////////////////////////////////////////////////////////
  crrnt_string = 0; 
  values.push_back(new_string);
  types.push_back(new_string);

  for(unsigned i=0; i<timeString.length(); i++) {

    timeString[i] = tolower(timeString[i]);

    if ((reading_digits) && isalpha(timeString[i])) {
      reading_digits = false;
      types[crrnt_string].push_back(timeString[i]);
    }
    else if ((reading_digits) && isdigit(timeString[i])) {
      values[crrnt_string].push_back(timeString[i]);
    }
    else if ((! reading_digits) && isalpha(timeString[i])) {
      types[crrnt_string].push_back(timeString[i]);
    }
    else if ((! reading_digits) && isdigit(timeString[i])) {
      reading_digits = true;
      values.push_back(new_string);
      types.push_back(new_string);
      crrnt_string++;
      values[crrnt_string].push_back(timeString[i]);
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Translate the time labels into seconds 
  //////////////////////////////////////////////////////////////////////
  total_seconds = 0;

  for(unsigned j = 0; j<values.size(); j++) {

    if ( types[j].compare( 0, types[j].size(), "seconds", 
                           0, types[j].size()) == 0 ) {
      total_seconds += atoi(values[j].c_str());
    }  
    else if ( types[j].compare( 0, types[j].size(), "minutes", 
                                0, types[j].size()) == 0 ) {
      total_seconds += 60*atoi(values[j].c_str());
    }  
    else if ( types[j].compare( 0, types[j].size(), "hours", 
                                0, types[j].size()) == 0 ) {
      total_seconds += 60*60*atoi(values[j].c_str());
    }  
    else if ( types[j].compare( 0, types[j].size(), "days", 
                                0, types[j].size()) == 0 ) {
      total_seconds += 24*60*60*atoi(values[j].c_str());
    }  
    else if ( types[j].compare( 0, types[j].size(), "weeks", 
                                0, types[j].size()) == 0 ) {
      total_seconds += 7*24*60*60*atoi(values[j].c_str());
    } 
    else {
       error("ERROR: %s is not a valid time string", types[j].c_str()); 
    }  
  }

  return(total_seconds); 
}

