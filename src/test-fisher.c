#include <stdlib.h> 
#include <math.h> 
#include <string.h> 
#include "subread.h" 
#include "HelperFunctions.h"



double factorial_float(int a)
{

        double ret = 0;
        while(a)
                ret += log(a--);
        return ret;
}



double fisherSub(int a, int b, int c, int d)
{
        double ret = factorial_float(a+b) + factorial_float(c+d) + factorial_float(a+c) + factorial_float(b+d) ;
        ret -= factorial_float(a) + factorial_float(b) + factorial_float(c) + factorial_float(d) + factorial_float(a+b+c+d);
        return pow(2.71828183, ret);
}




/**
 * See HELP string or run with no arguments for usage.
 * <p>
 * The code used to calculate a Fisher p-value comes originally from a
 * <a href="http://infofarm.affrc.go.jp/~kadasowa/fishertest.htm">JavaScript program</a>
 * by T. Kadosawa (kadosawa@niaes.affrc.go.jp).
 * Retrieved from http://www.users.zetnet.co.uk/hopwood/tools/StatTests.java on 3/Jul/2012
 *
 * @author David Hopwood
 * @date   2000/04/23
 */

double fisher_exact_test(int a, int b, int c, int d)
{

	if (a * d > b * c) {
	    a = a + b; b = a - b; a = a - b;
	    c = c + d; d = c - d; c = c - d;
	}

	if (a > d) { a = a + d; d = a - d; a = a - d; }
	if (b > c) { b = b + c; c = b - c; b = b - c; }

	double p_sum = 0.0;

	double p = fisherSub(a, b, c, d);
	while (a >= 0) {
	    p_sum += p;
	    if (a == 0) break;
	    --a; ++b; ++c; --d;
	    p = fisherSub(a, b, c, d);
	}

	return p_sum;
}


long double fastfact(int x){
	return logl(x)*x - x + 0.5 * logl(2*M_PI* x) + 1/(12.*x) - 1./(360.* x*x*x) +  1./(1260.* x*x*x*x*x) -  1./(1680.*x*x*x*x*x*x*x);// + (x>60?0:(1./(1188.*x*x*x*x*x*x*x*x*x ) ));
}

main(){
	unsigned int 	a = 10  , c = 11,
			b = 11 ,  d = 5000;

	double fisher, fisher_old;
	fisher = fast_fisher_test_one_side(a,b,c,d, NULL, 0);
	fisher_old = fisher_exact_test(a,b,c,d);
	printf("Log fisher = %.7f ; Old fisher = %.7f\n", log(fisher),log(fisher_old));

	long double x1 =  1E-19L + 1E-20L;
	long double x2 =  1L - expl(logl(0.5L) + logl(2.0L));
	printf("New Vals: x1=%LG, x2=%LG\n", x1,x2);
}

