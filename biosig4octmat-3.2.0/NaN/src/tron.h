/*

$Id$
Copyright (c) 2007-2009 The LIBLINEAR Project.
Copyright (c) 2010 Alois Schloegl <alois.schloegl@gmail.com>
This function is part of the NaN-toolbox
http://pub.ist.ac.at/~schloegl/matlab/NaN/

This code was extracted from liblinear-1.51 in Jan 2010 and 
modified for the use with Octave 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

*/

#ifndef _TRON_H
#define _TRON_H

class function
{
public:
	virtual double fun(double *w) = 0 ;
	virtual void grad(double *w, double *g) = 0 ;
	virtual void Hv(double *s, double *Hs) = 0 ;

	virtual int get_nr_variable(void) = 0 ;
	virtual ~function(void){}
};

class TRON
{
public:
	TRON(const function *fun_obj, double eps = 0.1, int max_iter = 1000);
	~TRON();

	void tron(double *w);
	void set_print_string(void (*i_print) (const char *buf));

private:
	int trcg(double delta, double *g, double *s, double *r);
	double norm_inf(int n, double *x);

	double eps;
	int max_iter;
	function *fun_obj;
	void info(const char *fmt,...);
	void (*tron_print_string)(const char *buf);
};
#endif
