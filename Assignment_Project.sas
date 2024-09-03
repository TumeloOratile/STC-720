*HWMA control chart;
proc iml;
m = 20;
n = 100;
x = j(n, m, .);

*case k;
lambda = 0.25;
xbar_0 = 0.05;		*randomly chosen;
s2_0 = 2;
xbar_vec = j(1, m, .);

*Simulate a multivariate normal sample X_jk;
call randseed(123);
do j = 1 to m; *loop through m samples;
	do k = 1  to n;	*loop through n obs.;
		x[k,j] = rand("normal", xbar_0, s2_0);
	end;
	*get sample means;
	xbar_vec[j] = x[,j][+]/n;
end;
print x, xbar_vec;

*Calculate HWMA statistics and control limits;
xbar_bar = xbar_0;
result = {};
do j =  1 to m;
	if j = 1 then do;
		hwma = lambda * xbar_vec[j] + (1-lambda) * xbar_bar;
		L = 2.5;
		ucl = xbar_0 + L*sqrt(lambda**2 * (s2_0/n));
		lcl = xbar_0 - L*sqrt(lambda**2 * (s2_0/n));
		cl = xbar_0;
	end;
	else do;
		xbar_bar = (xbar_bar + xbar_vec[j-1])/(j-1);
		hwma = lambda * xbar_vec[j] + (1-lambda) * xbar_bar;
		L = 2.5;
		ucl = xbar_0 + L*sqrt(lambda**2 * (s2_0/n) + (1-lambda)**2 * s2_0/(n*(j-1)));
		lcl = xbar_0 - L*sqrt(lambda**2 * (s2_0/n) + (1-lambda)**2 * s2_0/(n*(j-1)));
		cl = xbar_0;
	end;
	
	result = result // (j || xbar_vec[j] || xbar_bar || hwma || lcl || ucl || cl || L);
end;
nms = {"t" "jth_xbar" "xbar_bar" "hwma" "lcl" "ucl" "cl" "L"};
print result[colname = nms];

create _data from result [colname = nms];
append from result;
quit;

*SSS-HWMA control chart;
proc iml;
m = 20;
n = 100;
x = j(n, m, .);

*case k;
lambda = 0.25;
xbar_0 = 0.05;		*randomly chosen;
s2_0 = 2;
xbar_vec = j(1, m, .);

*Simulate a multivariate normal sample X_jk;
call randseed(123);
do j = 1 to m; *loop through m samples;
	do k = 1  to n;	*loop through n obs.;
		x[k,j] = rand("normal", xbar_0, s2_0);
	end;
	*get sample means;
	xbar_vec[j] = x[,j][+]/n;
end;
print x, xbar_vec;

*Calculate HWMA statistics and control limits;
xbar_bar = xbar_0;
result = {};
do j =  1 to m;
	if j = 1 then do;
		sss_hwma = lambda * xbar_vec[j] + (1-lambda) * xbar_bar;
		L_1 = 2.5;
		L_2 = 2.7;
		ucl = xbar_0 + L_1*sqrt(lambda**2 * (s2_0/n));
		lcl = xbar_0 - L_2*sqrt(lambda**2 * (s2_0/n));
		cl = xbar_0;
	end;
	else do;
		xbar_bar = (xbar_bar + xbar_vec[j-1])/(j-1);
		sss_hwma = lambda * xbar_vec[j] + (1-lambda) * xbar_bar;
		L_1 = 2.5;
		L_2 = 2.7;
		ucl = xbar_0 + L_1*sqrt(lambda**2 * (s2_0/n) + (1-lambda)**2 * s2_0/(n*(j-1)));
		lcl = xbar_0 - L_2*sqrt(lambda**2 * (s2_0/n) + (1-lambda)**2 * s2_0/(n*(j-1)));
		cl = xbar_0;
	end;
	
	result = result // (j || xbar_vec[j] || xbar_bar || sss_hwma || lcl || ucl || cl || L_1 || L_2);
end;
nms = {"t" "jth_xbar" "xbar_bar" "sss_hwma" "lcl" "ucl" "cl" "L_1" "L_2"};
print result[colname = nms];

create _data2 from result [colname = nms];
append from result;
quit;

ods graphics on/ height=400px width=600px;
proc sgplot data = work._data;
series x = t y = hwma/ lineattrs = (color = grey thickness = 2);
scatter x = t y = hwma/ markerattrs = (color = sandybrown size = 8 symbol = circlefilled);
series x=t y=ucl / lineattrs=(color=tomato pattern=dash thickness=2);
series x=t y=lcl / lineattrs=(color=forestgreen pattern=dash thickness=2);
refline cl/ lineattrs=(color = black pattern = dash thickness = 1);
xaxis label = "time (t)" values=(1 to 20 by 1);
yaxis label = "HWMA statistic";
title "HWMA control chart";
run;
ods graphics off;

ods graphics on/ height=400px width=600px;
proc sgplot data = work._data2;
series x = t y = sss_hwma/ lineattrs = (color = grey thickness = 2);
scatter x = t y = sss_hwma/ markerattrs = (color = sandybrown size = 8 symbol = circlefilled);
series x=t y=ucl / lineattrs=(color=skyblue pattern=dash thickness=2);
series x=t y=lcl / lineattrs=(color=violet pattern=dash thickness=2);
refline cl/ lineattrs=(color = black pattern = dash thickness = 1);
xaxis label = "time (t)" values=(1 to 20 by 1);
yaxis label = "SSS-HWMA statistic";
title "SSS-HWMA control chart";
run;
ods graphics off;



