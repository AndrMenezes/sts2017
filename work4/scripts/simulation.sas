proc delete data = _all_; run;

%let lambda	= 	1.0;
%let beta	=	2.0;
%let alpha	= 	0.1;
%let n		=	500;


data simulation;
	call streaminit(6666);
	do i = 1 to &n;
		if(i <= 100) then n1 = 1; else n1 = 0;
		if(i <= 300) then n2 = 1; else n2 = 0;
		if(i <= 500) then n3 = 1; else n3 = 0;
		ui 	= 	rand('UNIFORM');
		ti	=	1 / &lambda * (-log(ui / (&alpha + ui * (1 - &alpha))))**(1 / &beta);
		output;
	end;
run;

proc export data = simulation outfile = 'C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 4\Scripts\sim-MOEW.txt' replace dbms = tab; run;

* Estimativas via Máxima Verossimilhança;
%macro mles(dados =, n =);
	proc nlmixed data = &dados tech = quanew update = bfgs df = 99999;
		parms 	lambda = 1.0, beta = 2.0, alpha = 0.1;
		bounds 	lambda > 0, beta > 0, alpha > 0;

		where &n = 1;

		llike	 =	log(alpha) + log(beta) + beta * log(lambda) + (beta - 1) * log(ti)
					- (lambda * ti)**beta - 2 * log(1 - (1 - alpha) * exp(-(lambda * ti)**beta));
		model ti ~ general(llike);
		ods output  ParameterEstimates=emvs(drop = tvalue Probt alpha df gradient);
	run;
	proc export data = emvs outfile = "emvs&n..txt" replace dbms = tab; run;
%mend;

%mles(dados=simulation, n=n1);
%mles(dados=simulation, n=n2);
%mles(dados=simulation, n=n3);

* Estimativas via Bayesiana;
proc mcmc data = simulation outpost = post nmc = 100000 seed = 666 thin = 20 nbi = 100 diagnostics = all;
	parms 	lambda = 1.0 beta = 2.0 alpha = 0.1;

	where   n1 = 1;
	prior	lambda	~ 	gamma(shape = 0.0001, iscale = 0.0001); 
	prior	beta   	~ 	gamma(shape = 0.0001, iscale = 0.0001); 
	prior	alpha	~	uniform(0, 1);  

	llike	=	log(alpha) + log(beta) + beta * log(lambda) + (beta - 1) * log(ti)
				- (lambda * ti)**beta - 2 * log(1 - (1 - alpha) * exp(-(lambda * ti)**beta));
	model ti ~ general(llike);
	*ods output PostSumInt = bayes(drop = N);
run;
%macro bayes(dados=, n=);
	proc mcmc data = &dados outpost = post nmc = 100000 seed = 666 thin = 20 nbi = 100 diagnostics = all;
		parms 	lambda = 1.0 beta = 2.0 alpha = 0.1;

		where &n = 1;

		prior	lambda	~ 	gamma(shape = 0.0001, iscale = 0.0001); 
		prior	beta   	~ 	gamma(shape = 0.0001, iscale = 0.0001); 
		prior	alpha	~	uniform(0, 1);  

		llike	=	log(alpha) + log(beta) + beta * log(lambda) + (beta - 1) * log(ti)
					- (lambda * ti)**beta - 2 * log(1 - (1 - alpha) * exp(-(lambda * ti)**beta));

		model ti ~ general(llike);
		ods output PostSumInt = bayes(drop = N);
	run;
	proc export data = bayes outfile = "bayes&n..txt" replace dbms = tab; run;
	proc export data = post  outfile = "post&n..txt" replace dbms = tab; run;
%mend;

%bayes(dados=simulation, n=n1);
%bayes(dados=simulation, n=n2);
%bayes(dados=simulation, n=n3);
