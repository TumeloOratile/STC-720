proc iml;
********************HWMA Xbar chart*******************;
/**** Steady-State when tau is different from 1******/;
tau =2;
*Number of simulations;
sim = 20;
rlvec = j(sim,1,.);

*Size of the Phase II sample;
n = 100;

*lambda is the smothing constant;
lambda = 0.25; 
l1 = 1-lambda;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0.05;
 sigma0 = 2;

*Optimal parameter;
L=2.5;

*Shift: for the IC case, delta0=0 and delta1=0-for the OOC case, vary delta1 from 0.25 to 3 with an increment of 0.25;
delta0 = 0.00;		
delta1 = 0.00;		
stdev = sigma0;

mean0  = delta0 * (stdev);
mean1  = delta1 * (stdev);

df=5;
v=1;
wp=2;

do k=1 to sim; 
	count=0;
	signal=0;
	count2=0;
	vec={};
	
	zi_1 = mu0;
	vec={};

	do i=1 to 10000000000 until (signal=1); 
	    yi = j(n,1,.);
	    
	    * Generating observations from the Normal distribution;
	    yi0 = j(n,1,.);
		yi1 = j(n,1,.);
		
	    * Generating a Phase II sample;
	    * Generating observations from the Normal distribution;
	    yi0= randnormal(n,mean0,Sigma0);	*generates IC obs;
		yi1= randnormal(n,mean1,Sigma0);	*generates OOC obs;
		
		if (i<tau) then 
		yi = yi0; 
		
	    else yi = yi1;
	    
	    count=count+1;
		
	    *Control limits;
	    if i=1 then do; 
	    	lcl= mu0-L*sqrt(lambda**2*sigma0**2/n);
	        ucl= mu0+L*sqrt(lambda**2*sigma0**2/n);
			goto test;
		end;
		
	    lcl= mu0-L*sqrt(lambda**2*stdev**2/n+(l1)**2*sigma0**2/(n*(i-1)));
	    ucl= mu0+L*sqrt(lambda**2*stdev**2/n+(l1)**2*sigma0**2/(n*(i-1)));
	
		*Obtaining the charting (or Plotting) statistics;
	    test:
	    x_bar = sum(yi)/(n);
	    vec=vec//x_bar;
	
	   *HWMA Plotting statistic;
		zi = lambda*x_bar + (1-lambda)*zi_1;          
	
	    zi_1 = mean(vec);
		
	    *Comparing the plotting statistics to the control limits;
	    if (((zi>=ucl)|(zi<=lcl))&(i>tau)) then signal=1;
		else signal=0;
	
	    rlvec[k,1]=count-tau;
	end;
end;

results=rlvec;
print sim n lambda delta1 L;

create runlength_normal from results[colname={arl}];
append from results;
quit;

data  runlength_normal;
set work.runlength_normal;

proc univariate data= runlength_normal noprint;
var arl;		   
histogram;
inset mean std p5 q1 median q3 p95 / format = 10.2;
run;



proc iml;
********************SSS-HWMA Xbar chart*******************;
/**** Steady-State when tau is different from 1******/;
tau = 2;
*Number of simulations;
sim = 20;
rlvec = j(sim,1,.);
trendvec = j(sim,1,.);

*Size of the Phase II sample;
n = 100;

*lambda is the smothing constant;
lambda = 0.25; 
l1 = 1-lambda;

*Process parameters: Case K i.e. distribution parameters are known;
 mu0 = 0.05;
 sigma0 = 2;

*Optimal parameter;
L_1=2.8;
L_2=2.5;

*synthetic chart parameters;
trend_length = 3;	*no. of points required to confirm trend;
cl = mu0;
trend_threshold = cl;		*threshold for significant upward or downward shift;

*Shift: for the IC case, delta0=0 and delta1=0-for the OOC case, vary delta1 from 0.25 to 3 with an increment of 0.25;
delta0 = 0.00;		
delta1 = 0.5;		
stdev = sigma0;

mean0  = delta0 * (stdev);
mean1  = delta1 * (stdev);

df=5;
v=1;
wp=2;



do k=1 to sim; 
	count=0;
	signal=0;
	trend_up = 0;	*counter for upward trend;
	trend_down = 0;		*counter for downward trend;
	confirmed = 0;
	vec={};	
	zi_1 = mu0;
	do i=1 to 10000000000 until (signal=1 | confirmed = 1); 
	    yi = j(n,1,.);
	    
	    * Generating observations from the Normal distribution;
	    yi0 = j(n,1,.);
		yi1 = j(n,1,.);
		
	    * Generating a Phase II sample;
	    * Generating observations from the Normal distribution;
	    yi0= randnormal(n,mean0,Sigma0);	*generates IC obs;
		yi1= randnormal(n,mean1,Sigma0);	*generates OOC obs;
		
		if (i<tau) then yi = yi0; 
	    else yi = yi1;
	    
	    count=count+1;	*counts the number of obs before first shift;
		
	    *Control limits;
	    if i=1 then do; 
	    	lcl= mu0-L_1*sqrt(lambda**2*sigma0**2/n);
	        ucl= mu0+L_2*sqrt(lambda**2*sigma0**2/n);
			goto test;
		end;
		
	    lcl= mu0-L_1*sqrt(lambda**2*stdev**2/n+(l1)**2*sigma0**2/(n*(i-1)));
	    ucl= mu0+L_2*sqrt(lambda**2*stdev**2/n+(l1)**2*sigma0**2/(n*(i-1)));
	
		*Obtaining the charting (or Plotting) statistics;
	    test:
	    x_bar = sum(yi)/(n);
	    vec=vec//x_bar;
	
	   *SSS-HWMA Plotting statistic;
		zi = lambda*x_bar + (1-lambda)*zi_1;          
	    zi_1 = mean(vec);
		
	    *Comparing the plotting statistics to the control limits;
	    if (((zi>=ucl)|(zi<=lcl))&(i>tau))  then do;
		    signal=1;
			
			*monitoring trends in the syynthetic chart;
			if zi >= trend_threshold then trend_up = trend_up + 1;
			
			if zi >= trend_threshold then trend_down = trend_down + 1;
			
			if(trend_up >= trend_length | trend_down >= trend_length) then do;
				confirmed = 1;		*confirm the alarm after synthetic chart validation;
	    	end;
		end;
		else do;
			trend_up = 0;     * Reset trend counters if signal is not confirmed;
			trend_down = 0;
			signal = 0;
			confirmed = 0;
		end;
		rlvec[k,1] = count-tau;
		trendvec[k,1] = confirmed;
	end;
end;

results=rlvec;
print results trendvec, sim n lambda delta1 L_1 L_2;

create runlength_normal2 from results[colname={arl}];
append from results;
quit;

proc univariate data= runlength_normal2 noprint;
var arl;		   
histogram;
inset mean std p5 q1 median q3 p95 / format = 10.2;
run;




