proc delete data = _all_; run;

%include "C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 4\Scripts\decimal-nlmixed.sas";

proc import out 	 = diabetes replace dbms = tab
			datafile = 'C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 4\Scripts\diabetes.txt';
run;

data diabetes;
	set diabetes;
	if(left = right)  then di = 1;
	if(left ^= right) then di = 0; 
	if(left = 0) then delete;
run;

proc sort data = diabetes;
	by gender;
run;

/*Weibull*/
proc lifereg data = diabetes;
	by gender;
	model (left, right) = / dist = weibull;
run;

proc nlmixed data = diabetes tech = quanew maxiter = 10000;
	parms 	mu = 9.35, beta = 2.6;
	bounds 	mu > 0, beta > 0;

	by gender;	
	ui		=	right;
	li		=	left;

	if(di = 0) then do;
	Sli	  =	exp(-(li / mu)**beta);
	Sui	  =	exp(-(ui / mu)**beta);
	llike = log(Sli - Sui);
	end;

	if(di = 1) then do;
	llike =	log(beta) - beta * log(mu) + beta * log(ui) - (ui / mu)**beta;
	end;

	model ui ~ general(llike);
	estimate 'lambda' 1 / mu; 
	ods output ParameterEstimates=emvs1(drop = tvalue Probt alpha df gradient);
	ods output AdditionalEstimates=lambda1(drop = tvalue Probt alpha df gradient rename = (label = parameter));
run;

/*Marshall-Olkin Extended Weibull*/
proc nlmixed data = diabetes tech = quanew maxiter = 10000 df=999999;
	parms 	mu = 9.35, beta = 2.6, alpha = 0.5;
	bounds 	mu > 0, beta > 0, alpha >0, alpha < 1;

	by gender;	
	ui		=	right;
	li		=	left;
	
	if(di = 0) then do;
	Sli	  =	alpha * exp(-(li / mu)**beta) / (1 - (1 - alpha) * exp(-(li / mu)**beta));
	Sui	  =	alpha * exp(-(ui / mu)**beta) / (1 - (1 - alpha) * exp(-(ui / mu)**beta));
	llike = log(Sli - Sui);
	end;

	if(di = 1) then do;
	llike =	log(alpha) + log(beta) - beta * log(mu) + beta * log(ui) - (ui / mu)**beta
		  - 2 * log(1 - (1 - alpha) * exp(-(ui / mu)**beta));
	end;
	model ui ~ general(llike);
	estimate 'lambda' 1 / mu; 
	ods output  ParameterEstimates=emvs2(drop = tvalue Probt alpha df gradient);
	ods output  AdditionalEstimates=lambda2(drop = tvalue Probt alpha df gradient rename = (label = parameter));
run;

/*Exponentiated Weibull*/
proc nlmixed data = diabetes tech = quanew maxiter = 10000 df=999999;
	parms 	mu = 9, beta = 2, alpha = 2;
	bounds 	mu > 0, beta > 0, alpha >0;

	by gender;	
	ui		=	right;
	li		=	left;
	
	if(di = 0) then do;
		Sli	  =	1 - (1 - exp(-(li / mu)**beta))**alpha;
		Sui	  =	1 - (1 - exp(-(ui / mu)**beta))**alpha;
		llike = log(Sli - Sui);
	end;

	if(di = 1) then do;
		llike =	log(alpha) + log(beta) - beta * log(mu) + beta * log(ui)
			  - (ui / mu)**beta + (alpha - 1) * log(1 - exp(-(ui / mu)**beta));
	end;
	model ui ~ general(llike);
*	estimate 'lambda' 1 / mu; 
	ods output  ParameterEstimates=emvs3(drop = tvalue Probt alpha df gradient);
*	ods output  AdditionalEstimates=lambda2(drop = tvalue Probt alpha df gradient rename = (label = parameter));
run;

/*Odd Weibull*/
proc nlmixed data = diabetes tech = quanew maxiter = 10000 df=999999;
	parms 	mu = 9, beta = 2, alpha = 2;
	bounds 	mu > 0, beta > 0, alpha >0;

	by gender;	
	ui		=	right;
	li		=	left;
	
	if(di = 0) then do;
		Sli	  	=	1 / (1 + (exp(li / mu)**beta - 1)**alpha);
		Sui	  	=	1 / (1 + (exp(ui / mu)**beta - 1)**alpha);
		llike 	= log(Sli - Sui);
	end;

	if(di = 1) then do;
		llike 	=	;
	end;
	model ui ~ general(llike);
*	estimate 'lambda' 1 / mu; 
	ods output  ParameterEstimates=emvs3(drop = tvalue Probt alpha df gradient);
*	ods output  AdditionalEstimates=lambda2(drop = tvalue Probt alpha df gradient rename = (label = parameter));
run;

data emvs1;
	set emvs1 lambda1;
run;
proc sort data = emvs1;
	by gender;
run;
data emvs2;
	set emvs2 lambda2;
run;
proc sort data = emvs2;
	by gender;
run;

proc export data = emvs1 outfile = "emvs1.txt" replace dbms = tab; run;
proc export data = emvs2 outfile = "emvs2.txt" replace dbms = tab; run;

/*Weibull*/
proc mcmc data = diabetes outpost = post1 nmc = 1000000 seed = 666 thin = 100 nbi = 1500 diagnostics = all dic;
	parms 	mu = 30 beta = 4;

	prior	mu   	~ 	igamma(shape = 0.0001, scale = 0.0001); 
	prior	beta   	~ 	igamma(shape = 0.0001, scale = 0.0001); 

	by gender;	
	ui		=	right;
	li		=	left;

	if(di = 0) then do;
		Sli	  =	exp(-(li / mu)**beta);
		Sui	  =	exp(-(ui / mu)**beta);
		llike = log(Sli - Sui);
	end;

	if(di = 1) then do;
		llike =	log(beta) - beta * log(mu) + beta * log(ui) - (ui / mu)**beta;
	end;

	model ui ~ general(llike);
	ods output PostSumInt = bayes1(drop = N);
	ods output DIC=dic1(where=(Criterion="DIC (smaller is better)"));
run;


/*Marshall-Olkin Extended Weibull*/
proc mcmc data = diabetes outpost = post2 nmc = 1000000 seed = 666 thin = 100 nbi = 1500 diagnostics = all dic;
	parms 	mu = 30 beta = 4 alpha = 0.5;

	prior	mu   	~ 	uniform(0, 1000); *igamma(shape = 0.0001, scale = 0.0001); 
	prior	beta   	~ 	uniform(0, 1000); *igamma(shape = 0.0001, scale = 0.0001); 
	prior	alpha	~	uniform(0, 1); 

	by gender;	
	ui		=	right;
	li		=	left;

	if(di = 0) then do;
		Sli	  =	alpha * exp(-(li / mu)**beta) / (1 - (1 - alpha) * exp(-(li / mu)**beta));
		Sui	  =	alpha * exp(-(ui / mu)**beta) / (1 - (1 - alpha) * exp(-(ui / mu)**beta));
		llike = log(Sli - Sui);
	end;

	if(di = 1) then do;
		llike =	log(alpha) + log(beta) - beta * log(mu) + beta * log(ui) - (ui / mu)**beta
			  - 2 * log(1 - (1 - alpha) * exp(-(ui / mu)**beta));
	end;

	model ui ~ general(llike);
	ods output PostSumInt = bayes2(drop = N);
	ods output DIC=dic2(where=(Criterion="DIC (smaller is better)"));
run;

proc print data = dic1;
	format _numeric_ comma10.4;
run; 
proc print data = dic2;
	format _numeric_ comma10.4;
run; 

proc export data = bayes1 outfile = "bayes1.txt" replace dbms = tab; run;
proc export data = post1  outfile = "post1.txt" replace dbms = tab; run;
proc export data = bayes2 outfile = "bayes2.txt" replace dbms = tab; run;
proc export data = post2  outfile = "post2.txt" replace dbms = tab; run;


