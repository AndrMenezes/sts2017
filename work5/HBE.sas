proc delete data = _all_; run;

libname "C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 5\R";
proc import out = dados replace dbms = csv
			datafile = "exemplo-HBE.csv";
run;

proc import out = boot replace dbms = csv
			datafile = "boot-samples-HBE.csv";
run;

/*Estimativas de máxima verossimilhança para a amostra observada*/
proc nlmixed data = dados tech = quanew update = bfgs df = 999999;
	parms mu1 = 5, mu2 = 10, r = 8;
	bounds mu1 > 0, mu2 > 0, r > 0;

	Sxy	=	exp(-((x / mu1) ** r + (y / mu2) ** r) ** (1 / r));
	fxy	=	(x * y) ** (r - 1) / ((mu1 * mu2)**r) * ((x / mu1) ** r + (y / mu2) ** r)**(1 / r - 2) *
			(r - 1 + ((x / mu1) ** r + (y / mu2) ** r)**(1 / r)) * Sxy;

	llike	=	log(fxy);

	model x ~ general(llike);

	estimate "R" 	mu2 ** r / (mu1**r + mu2**r);
	estimate "rho" 	gamma(1 / r)**2 / (r * gamma(2 / r)) - 1;
run;

/*Estimativas de bayesianas para a amostra observada*/
proc mcmc data = dados outpost = post nmc = 1000000 seed = 666 thin = 100 nbi = 1500 
		  diagnostics = all dic monitor = (mu1 mu2 r Rxy rho);
	parms 	mu1 = 2 mu2 = 3 r = 2;

	prior	mu1	~ 	gamma(shape = 0.0001, iscale = 0.0001); 
	prior	mu2 ~ 	gamma(shape = 0.0001, iscale = 0.0001); 
	prior 	r	~ 	gamma(shape = 0.0001, iscale = 0.0001);

	Sxy	=	exp(-((x / mu1) ** r + (y / mu2) ** r) ** (1 / r));
	fxy	=	(x * y) ** (r - 1) / ((mu1 * mu2)**r) * ((x / mu1) ** r + (y / mu2) ** r)**(1 / r - 2) *
			(r - 1 + ((x / mu1) ** r + (y / mu2) ** r)**(1 / r)) * Sxy;

	ll	=	log(fxy);

	model general(ll);

	Rxy = mu2 ** r / (mu1**r + mu2**r);
	rho = gamma(1 / r)**2 / (r * gamma(2 / r)) - 1;
run;

/* 
data boot;
	set boot;
	if(B ^= 1 and B^= 2) then delete;
run; 
*/

options nonotes;
ods results off;
/*Estimativas de máxima verossimilhança para as amostras bootstrap*/
proc nlmixed data = boot tech = quanew update = bfgs df = 9999999;
	parms mu1 = 5.2379, mu2 = 10.5160, r = 7.5571;
	bounds mu1 > 0, mu2 > 0, r > 0;

	by	B;

	bounds mu1 > 0, mu2 > 0, r > 0;

	Sxy	=	exp(-((x / mu1) ** r + (y / mu2) ** r) ** (1 / r));
	fxy	=	(x * y) ** (r - 1) / ((mu1 * mu2)**r) * ((x / mu1) ** r + (y / mu2) ** r)**(1 / r - 2) *
			(r - 1 + ((x / mu1) ** r + (y / mu2) ** r)**(1 / r)) * Sxy;

	llike	=	log(fxy);

	model x ~ general(llike);

	estimate "rho" 	gamma(1 / r)**2 / (r * gamma(2 / r)) - 1;
	estimate "R" 	mu2 ** r / (mu1**r + mu2**r);

	ods output ParameterEstimates  = parmsest(drop = df tValue Probt Alpha Lower Upper Gradient);
	ods output AdditionalEstimates = rest(drop = df tValue Probt Alpha Lower Upper );
run;
options notes;
ods results on;

proc export data = parmsest replace outfile = "parms-HBE.csv" dbms = csv replace; run;
proc export data = rest     replace outfile = "Rrho-HBE.csv"  dbms = csv replace; run;
