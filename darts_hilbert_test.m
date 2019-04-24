%% test instantaneous hilbert measures
	srate = 512;
	t = 0:(1/srate):5;
	freqTS = abs(interp1(linspace(t(1),t(end),10),.5+rand(1,10),t,'spline'));
	centfreq = mean(freqTS);
	k = (centfreq/srate)*2*pi/centfreq;
	sig = sin(2*pi.*centfreq.*t + k*cumsum(freqTS-centfreq));
	amp_shift = [1:-(1/length(t)):0];
	sig = sig.*amp_shift(1:length(sig));
	hilsig = hilbert(sig);
	
	% amp
	figure;
	subplot(311)
	amp = abs(hilsig);
	plot(t,sig); hold on
	plot(t,amp)
	
	% phase
	subplot(312)
	phases = angle(hilsig);
	plot(t,sig); hold on
	plot(t,phases)
	phase_c = cos(phases);
	phase_s = sin(phases);
	
	% freq
	subplot(313)
	freq = srate*diff(unwrap(phases))/(2*pi);
	plot(t(1:end-1),freq)
	hold on
	plot(t,freqTS,'r')
	set(gca,'ylim',[0 2])