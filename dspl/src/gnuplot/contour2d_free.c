/*
* Copyright (c) 2015-2024 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of DSPL.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"

void DSPL_API contour2d_free(contour2d_t* c)
{
    int nl;
	if(c->lines)
	{
		for (nl = 0; nl < c->nlines; nl++)
		{
			if(c->lines[nl].points)
				free(c->lines[nl].points);
		}
		free(c->lines);
		c->lines = NULL;
	}
	c->nlines = 0;
	c->level = 0.0;
}