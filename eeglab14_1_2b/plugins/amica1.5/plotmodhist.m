function plotmodhist(x,mod,ch,srate,Fmax)
% function plotmodhist(data,mod,chanlocs,srate,Fmax)
%
% Plot the topoplot, density, log density, and spectra for the components
% in the AMICA model mod, with density calculated frm data x. Use chanlocs
% for channel locations in topoplot, and srate for spectrum axis labeling
% up to frequency Fmax (default srate/2).
%
% Commands:
%           <enter> - advance to next component
%           0       - go back to previous component
%           n       - go to component n in this model
%           m000    - go to first comonent of model m
%


if max(size(mod.LLt)) > 1
    x = x(:,find(mod.LLt(1,:)~=0));
    mod.LLt = mod.LLt(:,find(mod.LLt(1,:)~=0));
end

[nx,N] = size(x);

c = mod.c;
W = mod.W;
%A = mod.A;
S = mod.S;
gm = mod.mod_prob;
alpha = mod.alpha;
mu = mod.mu;
sbeta = mod.sbeta;
rho = mod.rho;
mn = mod.data_mean;

pdftype = 1; %mod.pdftype;

m = size(alpha,1);
n = size(c,1);
M = length(gm);

if M == 1
    LLt = ones(1,N);
end
nbins = min(100,sqrt(N));
psdwindow = srate;
psdoverlap = round(srate/2);
nfft = srate;
if nargin < 5
    Fmax = round(srate/2);
end
fftfactor = 0.7;


s = zeros(n,N);


B = 1000;
num_blocks = floor(N/B);

% sphere the data
for k = 1:num_blocks
    xstrt = (k-1)*B + 1;
    if k < num_blocks
        xstp = k*B;
    else
        xstp = N;
    end
    for i = 1:nx
        x(i,xstrt:xstp) = x(i,xstrt:xstp) - mn(i);
        %x(i,xstrt:xstp) = s(1,xstrt:xstp) - mn(i);
    end    
    x(1:n,xstrt:xstp) = S(1:n,:) * x(:,xstrt:xstp);
end

f = figure;

[mval,mind] = max(mod.LLt);

for h = 1:M
    A(:,:,h) = pinv(W(:,:,h)*S(1:n,:));
end

h = 1;
while h <= M

    if M > 1
        inds = find(mind == h);
        %vsum = sum(mod.v,2);
        if isempty(inds)
            h = h + 1;
            continue;
        end         
    end
    
    for k = 1:num_blocks
        xstrt = (k-1)*B + 1;
        if k < num_blocks
            xstp = k*B;
        else
            xstp = N;
        end
        s(:,xstrt:xstp) = W(:,:,h) * x(1:n,xstrt:xstp);
        for i = 1:n
            s(i,xstrt:xstp) = s(i,xstrt:xstp) - c(i,h);
        end
    end
    
    i = 1;
    while i <= n
        if M > 1
            smn = mean(s(i,inds));           
            [hst,b] = hist(s(i,inds)-smn,nbins);
            %umax = max(s(i,inds));
            %umin = min(s(i,inds));
            %deltau = (umax-umin)/nbins;
            %b = (umin+deltau/2):deltau:(umax-deltau/2);
            %ds = min(s(i,:)-umin,umax-umin);
            %ds = max(0,ds);
            %u(i,:) = 1 + round((nbins - 1) * ds / (umax - umin));
            %for k = 1:nbins
            %    hst(k) = sum(mod.v(h,:).*(u(i,:)==k))/vsum(h);
            %end
            %%%%%hst  = diff([0 find(diff(sort(u(i,:)))) N]);

        else
            [hst,b] = hist(s(i,:),nbins);
        end
        db = abs(b(2)-b(1));
        hst = hst / (sum(hst)*db);

        p = zeros(1,length(b));

        clf;
        subplot(2,2,1);
        topoplot(A(:,i,h),ch);
        if M > 1
            title(['model ' int2str(h) ' component ' int2str(i) ' -- ' num2str(length(inds)*100/N) '% of points']);        
        else
            title(['Component ' int2str(i)]);
        end
        
        subplot(2,2,2);  plot(0);
        axis([b(1)-2*db b(end)+2*db 0 max(hst)+0.1]); hold on;
        if M > 1
            title(['model ' int2str(h) ' component ' int2str(i) ' -- ' num2str(length(inds)*100/N) '% of points']);
        else
            title(['Component ' int2str(i)]);
        end

        for j = 1:m
            y = (b-mu(j,i,h)) * sbeta(j,i,h);
            switch pdftype
                case 2
                    d = alpha(j,i,h) * sbeta(j,i,h) * exp(-0.5*y.^2) / sqrt(2*pi);
                case 1
                    d = alpha(j,i,h) * sbeta(j,i,h) * exp(-abs(y).^rho(j,i,h)) / (2*gamma(1+1/rho(j,i,h)));
                case 0
                    d = alpha(j,i,h) * sbeta(j,i,h) * 0.25 ./ cosh(y/2).^2;
            end
            plot(b,d,'g');
            p = p + d;
        end
        
        plot(b,hst,'r'); plot(b,p);
        

        fftmaxind = round((nfft/2+1)*Fmax/(srate/2));
        %linds = length(inds);
        if M > 1 && length(inds) < N - psdwindow
            [sp1,f1] = pwelch(s(i,inds) * norm(A(:,i,h)),psdwindow,psdoverlap,nfft,srate);
            [sp2,f2] = pwelch(s(i,setdiff(1:N,inds)) * norm(A(:,i,h)),psdwindow,psdoverlap,nfft,srate);
            subplot(2,2,3); plot(0); plot(f1(1:fftmaxind),10*log(sp1(1:fftmaxind)),'b'); hold on;
            plot(f2(1:fftmaxind),10*log(sp2(1:fftmaxind)),'r');

            [sp3,f3] = pwelch(s(i,:) * norm(A(:,i,h)),psdwindow,psdoverlap,nfft,srate);
            plot(f3(1:fftmaxind),10*log(sp3(1:fftmaxind)),'g');
            title('Spectrum of model and non-model time points');
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            axis([0 Fmax -Inf Inf]);
        else
            subplot(2,2,3);
            [sp,f] = pwelch(s(i,:) * norm(A(:,i,h)),psdwindow,psdoverlap,nfft,srate);
            plot(log10(f(1:fftmaxind)),10*log(sp(1:fftmaxind)));
            title('Spectrum of component activations');
            xlabel('log10 Frequency (Hz)'); ylabel('Power (dB)');

        end        
        
        subplot(2,2,4); hold on;
        axis([b(1)-2*db b(end)+2*db -Inf Inf]);
        plot(b,log(hst),'r'); plot(b,log(p));
        title('Log histogram and component density model');

        sum(p)*db
        
        inp = input('');
        if inp == 0
            if i > 1
                i = i - 1;
            else
                if h > 1
                    h = h - 1;
                    i = n;
                end
            end                    
        elseif (inp >= 1 & inp <= n)
            i = inp
        elseif (inp >= 1000)
            h = round(inp / 1000);
            break
        else
            if i < n
                i = i + 1;
            else
                h = h + 1;
                break
            end
        end
    end
end