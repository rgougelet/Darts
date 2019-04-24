%% Inits
sig_length = 100;
srate = 512;
n = srate*sig_length;
t = linspace(0,sig_length,n);
nyq = srate/2;

%% Test signals
const = zeros(1,n);
sine = nmlz(sin(1:n));
two_sine = nmlz(sin(t)+sin(5*t));
sinec = nmlz(sin(1:n)./(1:n));
ssine = nmlz(0.1*sin(1:n));
chir = nmlz(chirp(t,3,40,14,'quadratic'));
saw = nmlz(sawtooth(1:n));
unif = nmlz(rand(1,n));
sunif = 0.5*nmlz(rand(1,n));
gau = nmlz(randn(1,n));
sgau = 0.5*nmlz(randn(1,n));
bimod = nmlz([randn(1,n/2)-10, randn(1,n/2)+10]);
trimod = nmlz([randn(1,round(n/3))-10,randn(1,round(n/3)), randn(1,round(n/3))+10]);
squ = nmlz(square(2*pi*0.25*t));
diri = nmlz(diric(1:n,10));
voltco = nmlz(vco(sawtooth(2*pi*1:n,0.75),[0.1 0.4]*srate,srate));
onoff = nmlz([ones(1,n/2),zeros(1,n/2)]);
lepto = nmlz(pearsrnd(0,1,0,100,1,n));
noise_sine = nmlz(sin(1:n)+randn(1,n));
noise_lepto = nmlz(lepto+randn(1,n));
noise_squ = nmlz(square(2*pi*0.25*t)+randn(1,n));

%% Entropy
disp('Entropy...........')
disp(['  Constant ', num2str(entropy(				const)) ] )
disp(['  Sine ', num2str(entropy(				sine)) ] )
disp(['  2Sine ', num2str(entropy(				two_sine)) ] )
disp(['  Sinc ', num2str(entropy(				sinec)) ] )
disp(['  Scaled sine ', num2str(entropy(	ssine)) ] )
disp(['  Chirp ', num2str(entropy(				chir)) ] )
disp(['  Sawtooth ', num2str(entropy(		saw	)) ] )
disp(['  Uniform ', num2str(entropy(			unif)) ] )
disp(['  Scaled Uniform ', num2str(entropy(			sunif)) ] )
disp(['  Gaussian ', num2str(entropy(		gau)) ] )
disp(['  Scaled Gaussian ', num2str(entropy(		sgau)) ] )
disp(['  Bimodal ', num2str(entropy(			bimod)) ] )
disp(['  Trimodal ', num2str(entropy(			trimod)) ] )
disp(['  Square ', num2str(entropy(			squ)) ] )
disp(['  Dirichlet ', num2str(entropy(		diri)) ] )
disp(['  VCO ', num2str(entropy(		voltco)) ] )
disp(['  On/Off ', num2str(entropy(		onoff)) ] )
disp(['  Lepto ', num2str(entropy(		lepto)) ] )

%% Aggregate entropy
disp('Aggregate entropy...........')
disp(['  Sine ', num2str(entropy(				hist(			sine					,nyq)))])
disp(['  Scaled sine ', num2str(entropy(	hist(			ssine						,nyq)))])
disp(['  Sawtooth ', num2str(entropy(		hist(			saw					,nyq)))])
disp(['  Uniform ', num2str(entropy(			hist(			unif							,nyq)))])
disp(['  Gaussian ', num2str(entropy(		hist(			gau						,nyq)))])
disp(['  Bimodal ', num2str(entropy(			hist(			bimod	,nyq)))])
disp(['  Trimodal ', num2str(entropy(		hist(			trimod	,nyq)))])
disp(['  Square ', num2str(entropy(			hist(			squ, nyq)))])
disp(['  Dirichlet ', num2str(entropy(		hist(			diri, nyq)))])
disp(['  VCO ', num2str(entropy(					hist(			voltco, nyq)))])
disp(['  On/Off ', num2str(entropy(			hist(			onoff, nyq)))])

%% Spectral entropy
disp('Spectral entropy...........')
disp(['  Constant ', num2str(entropy(				pspectrum(			sine					,srate)))])
disp(['  Sine ', num2str(entropy(				pspectrum(			sine					,srate)))])
disp(['  2Sine ', num2str(entropy(				pspectrum(			two_sine					,srate)))])
disp(['  Sinc ', num2str(entropy(				pspectrum(			sinec					,srate)))])
disp(['  Scaled sine ', num2str(entropy(	pspectrum(			ssine						,srate)))])
disp(['  Chirp ', num2str(entropy(				pspectrum(			chir					,srate)))])
disp(['  Sawtooth ', num2str(entropy(		pspectrum(			saw					,srate)))])
disp(['  Uniform ', num2str(entropy(			pspectrum(			unif							,srate)))])
disp(['  Scaled Uniform ', num2str(entropy(			pspectrum(			sunif							,srate)))])
disp(['  Gaussian ', num2str(entropy(		pspectrum(			gau						,srate)))])
disp(['  Scaled Gaussian ', num2str(entropy(		pspectrum(			sgau						,srate)))])
disp(['  Bimodal ', num2str(entropy(			pspectrum(			bimod	,srate)))])
disp(['  Trimodal ', num2str(entropy(		pspectrum(			trimod	,srate)))])
disp(['  Square ', num2str(entropy(			pspectrum(			squ, srate)))])
disp(['  Dirichlet ', num2str(entropy(		pspectrum(			diri, srate)))])
disp(['  VCO ', num2str(entropy(					pspectrum(			voltco, srate)))])
disp(['  On/Off ', num2str(entropy(			pspectrum(			onoff, srate)))])
disp(['  Lepto ', num2str(entropy(			pspectrum(			lepto, srate)))])
disp(['  Noise+Sine ', num2str(entropy(			pspectrum(			noise_sine, srate)))])
disp(['  Noise+Square ', num2str(entropy(			pspectrum(			noise_squ, srate)))])

%% Persistence spectral entropy
disp('Persistence spectral entropy...........')
disp(['  Sine ', num2str(entropy(				pspectrum(			sine					,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1, 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Scaled sine ', num2str(entropy(	pspectrum(			ssine						,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1, 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Sawtooth ', num2str(entropy(		pspectrum(			saw					,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Uniform ', num2str(entropy(			pspectrum(			unif							,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Gaussian ', num2str(entropy(		pspectrum(			gau						,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Bimodal ', num2str(entropy(			pspectrum(			bimod	,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Trimodal ', num2str(entropy(		pspectrum(			trimod	,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Square ', num2str(entropy(			pspectrum(			squ ,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  Dirichlet ', num2str(entropy(		pspectrum(			diri ,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  VCO ', num2str(entropy(					pspectrum(			voltco ,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])
disp(['  On/Off ', num2str(entropy(			pspectrum(			onoff ,srate,'persistence', 'FrequencyLimits', [0 60], 'FrequencyResolution',0.1)))])