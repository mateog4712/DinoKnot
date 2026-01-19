/***************************************************************************
                          constants.h  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define NONE            'N'         // no structure
#define HAIRP           'H'         // closes a hairpin loop
#define INTER           'I'         // closes an internal loop
#define MULTI           'M'         // closes a regular multi-loop
#define outer           'o'         // closes a "multiloop" with no inner crossing base pair, just regions on the left and right

#define M_WM            'B'         // closes a regular partial multi-loop
#define M_WMv            'v'         // closes a regular partial multi-loop
#define M_WMp            'p'         // closes a regular partial multi-loop
#define M_VMp            'm'         // Closes a partial multiloop across two RNA
#define M_VMc            'c'         // region for a multiloop across two RNA

#define FREE            'W'         // this base is free to be paired or unpaired to other base
#define LOOP            'V'         // closes a loop

#define P_WMB			'R'
#define P_VP			'D'
#define P_WI			'G'
#define P_BE			'J'
#define P_WIP			'L'
#define P_WMBP			'T'
#define P_WMBW       'X'
#define P_VPL        'Y'
#define P_VPR        'Z'

#endif

