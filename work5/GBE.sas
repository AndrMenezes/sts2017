proc delete data = _all_; run;

proc import out = dados replace dbms = csv
			datafile = "C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 5\R\exemplo2.csv";
run;

proc import out = boot replace dbms = csv
			datafile = "C:\Users\User\Dropbox\4° Série\Tópicos Especiais em Estatística\Trabalhos\Trabalho 5\R\boot-samples.csv";
run;

/*Estimativas de máxima verossimilhança para a amostra observada*/
 
proc nlmixed data = dados;
	parms lambda1 = 2, lambda2 = 3, theta = 0.8;
	bounds lambda1 > 0, lambda2 > 0, theta > 0 , theta < 1;

	Sxy = 	exp(- (lambda1 * x + lambda2 * y + theta * lambda1 * lambda2 * x * y));
	fxy	=	((1 - theta) * lambda1 * lambda2 + theta * lambda1 ** 2 * lambda2 * x + theta * lambda1 * lambda2 ** 2 * y + 
            theta ** 2 * lambda1 ** 2 * lambda2 ** 2 * x * y) * Sxy;

	ll	=	log(fxy);

	model x ~ general(ll);
	
	p	=	(lambda1 + lambda2) / (sqrt(2 * theta * lambda1 * lambda2));
	aux =	cdf ("normal", p, 0, 1);
	pi	=	constant("pi");
	z	=	1 / theta;
	u	=	0.65;
  	A   = 	log((0.56146 / z + u) * (1 + z));
 	B   = 	z ** 4 * exp(7.7 * z) * (2 + z) ** (3.7);
  	Ei  =	(A ** (-7.7) + B) ** (-0.13);
	estimate "R" 1 / 2 + sqrt(pi) * (lambda1 - lambda2) / (2 * sqrt(theta * lambda1 * lambda2)) *
							exp((lambda1 + lambda2)**2 / (4 * theta * lambda1 * lambda2)) * (1 - aux); 
	estimate "rho" z * exp(z) * Ei - 1;

	ods output ParameterEstimates  = parmsest(drop = df tValue Probt Alpha Lower Upper Gradient);
	*ods output AdditionalEstimates = rest(drop = df tValue Probt Alpha Lower Upper );
run;

/*Estimativas de bayesianas para a amostra observada*/

proc mcmc data = dados outpost = post nmc = 1000000 seed = 666 thin = 100 nbi = 1500 
		  diagnostics = all dic monitor = (lambda1 lambda2 theta R rho);
	parms 	lambda1 = 2 lambda2 = 3 theta = 0.8;

	prior	lambda1	~ 	igamma(shape = 0.0001, scale = 0.0001); 
	prior	lambda2 ~ 	igamma(shape = 0.0001, scale = 0.0001); 
	prior 	theta	~	beta(1, 1);

	Sxy = 	exp(- (lambda1 * x + lambda2 * y + theta * lambda1 * lambda2 * x * y));
	fxy	=	((1 - theta) * lambda1 * lambda2 + theta * lambda1 ** 2 * lambda2 * x + theta * lambda1 * lambda2 ** 2 * y + 
            theta ** 2 * lambda1 ** 2 * lambda2 ** 2 * x * y) * Sxy;

	ll	=	log(fxy);

	p	=	(lambda1 + lambda2) / (sqrt(2 * theta * lambda1 * lambda2));
	aux =	cdf ("normal", p, 0, 1);
	pi	=	constant("pi");
	R	= 	1 / 2 + sqrt(pi) * (lambda1 - lambda2) / (2 * sqrt(theta * lambda1 * lambda2)) *
							exp((lambda1 + lambda2)**2 / (4 * theta * lambda1 * lambda2)) * (1 - aux); 
	z	=	1 / theta;
	u	=	0.65;
  	A   = 	log((0.56146 / z + u) * (1 + z));
 	B   = 	z ** 4 * exp(7.7 * z) * (2 + z) ** (3.7);
  	Ei  =	(A ** (-7.7) + B) ** (-0.13);
	rho	=	z * exp(z) * Ei - 1;
	model general(ll);
	ods output PostSumInt = bayes(drop = N);
run;

/*
data boot;
	set boot;
	if(B ^= 1 and B^=2 and B^=10 ) then delete;
run;
*/

options nonotes;
ods results off;
/*Estimativas de máxima verossimilhança para as amostras bootstrap*/
proc nlmixed data = boot tech = quanew update = bfgs df = 9999999;
	parms lambda1 = 2.2785, lambda2 = 3.1237, theta = 0.8254;
	bounds lambda1 > 0, lambda2 > 0, theta > 0 , theta < 1;

	by	B;

	Sxy = 	exp(- (lambda1 * x + lambda2 * y + theta * lambda1 * lambda2 * x * y));
	fxy	=	((1 - theta) * lambda1 * lambda2 + theta * lambda1 ** 2 * lambda2 * x + theta * lambda1 * lambda2 ** 2 * y + 
            theta ** 2 * lambda1 ** 2 * lambda2 ** 2 * x * y) * Sxy;

	ll	=	log(fxy);

	model x ~ general(ll);
	
	p	=	(lambda1 + lambda2) / (sqrt(2 * theta * lambda1 * lambda2));
	aux =	cdf ("normal", p, 0, 1);
	pi	=	constant("pi");
	z	=	1 / theta;
	u	=	0.65;
  	A   = 	log((0.56146 / z + u) * (1 + z));
 	K   = 	z ** 4 * exp(7.7 * z) * (2 + z) ** (3.7);
  	Ei  =	(A ** (-7.7) + K) ** (-0.13);
	estimate "R" 1 / 2 + sqrt(pi) * (lambda1 - lambda2) / (2 * sqrt(theta * lambda1 * lambda2)) *
							exp((lambda1 + lambda2)**2 / (4 * theta * lambda1 * lambda2)) * (1 - aux); 
	estimate "rho" z * exp(z) * Ei - 1;

	ods output ParameterEstimates  = parmsest(drop = df tValue Probt Alpha Lower Upper Gradient);
	ods output AdditionalEstimates = rest(drop = df tValue Probt Alpha Lower Upper );
run;
options notes;
ods results on;

proc export data = parmsest replace outfile = "parms.csv" dbms = csv replace; run;
proc export data = rest     replace outfile = "Rrho.csv"  dbms = csv replace; run;
