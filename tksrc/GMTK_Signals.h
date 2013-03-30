/*-
 * GMTK_Signals.h
 *     Handling of Unix signals
 *
 * Written by Chris Bartels & Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 *
 *
 */

#ifndef GMTK_SIGNALS_H 
#define GMTK_SIGNALS_H 

void InstallSignalHandlers();
void InstallSignalHandlersTime();
bool TerminateSignalReceived();

#endif

