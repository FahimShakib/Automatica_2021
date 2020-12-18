function [u,freq] = idinput_custom(N,type,band,levels,nosine)
%IDINPUT Generates input signals for identification.
%   U = IDINPUT(N,TYPE,BAND,LEVELS)
%
%   U: The generated input signal. A column vector or a N-by-nu matrix.
%   N: The length of the input.
%   N = [N1 Nu] gives a N1-by-Nu input matrix (Nu input channels).
%   N = [P Nu M] gives a M*P-by-Nu input, periodic with period P
%       and with M periods.
%   Default values are Nu = 1 and M = 1.
%
%   TYPE: One of the following:
%         'RGS': Generates a Random, Gaussian Signal.
%         'RBS': Generates a Random, Binary Signal.
%         'PRBS': Generates a Pseudo-random, Binary Signal.
%         'SINE': Generates a sum-of-sinusoid signal.
%         Default: TYPE = 'RBS'.
%   BAND: A 1 by 2 row vector that defines the frequency band for the
%         input's frequency contents.
%         For the 'RS', 'RBS' and 'SINE' cases BAND = [LFR,HFR], where
%         LFR and HFR are the lower and upper limits of the passband,
%         expressed in fractions of the Nyquist frequency (thus always
%         numbers between 0 and 1).
%         For the 'PRBS' case BAND = [0,B], where B is such that the
%         signal is constant over intervals of length 1/B (the Clock Period).
%         Default: BAND =[0 1].
%   LEVELS = [MI, MA]: A 1 by 2 row vector, defining the input levels.
%         For 'RBS', 'PRBS', and 'SINE', the levels are adjusted so
%         that the input signal always is between MI and MA.
%         For the 'RGS' case, MI is the signal's mean value minus one
%         standard deviation and MA is the signal's mean plus one standard
%         deviation.
%         Default LEVELS = [-1 1]. For 'RGS', this means a signal with zero
%         mean and 1 standard deviation.
%
%   In the 'PRBS' case, if M > 1, the length of the data sequence and the
%   period is adjusted so that always an integer number of maximum length
%   PRBS periods are obtained. If M = 1 the period is chosen so that it
%   becomes longer than P = N. In the multiinput case the signals are
%   maximally shifted. This means that P/Nu is an upper bound for the model
%   orders that can be used to identify systems excited by such a signal.
%
%   In the 'SINE' case, the sinusoids are chosen from the frequency grid
%   freq = 2*pi*[1:Grid_Skip:fix(P/2)]/P intersected with pi*[BAND(1)
%   BAND(2)]. (for Grid_Skip see below.) For multi-input signals, the
%   different inputs use different frequencies from this grid. An integer
%   number of full periods is always delivered. The selected frequencies
%   are obtained as [U,FREQS] = IDINPUT(....), where row ku of FREQS
%   contains the frequencies of input number ku. The resulting signal
%   is affected by a 5th input argument SINEDATA:
%
%   U = IDINPUT(N,TYPE,BAND,LEVELS,SINEDATA)
%
%   where:
%   SINEDATA = [No_of_Sinusoids, No_of_Trials, Grid_Skip],
%   meaning that No_of_Sinusoids are equally spread over the indicated
%   BAND, trying No_of_Trials different, random, relative phases,
%   until the lowest amplitude signal is found.
%   Default SINEDATA = [10,10,1];
%
%   See also IDMODEL/SIM, IDDATA.

%   L. Ljung 3-3-95
%   Copyright 1986-2017 The MathWorks, Inc.

freq = [];
narginchk(1,5)

if nargin < 5
   nosine=[];
end

if nargin < 4
   levels = [];
end

if nargin < 3
   band = [];
end

if nargin < 2
   type = [];
elseif isstring(type)
   type = char(type);
end

if ~idpack.isNonnegIntMatrix(N) || ~isvector(N)
   ctrlMsgUtils.error('Ident:dataprocess:idinput1')
end

if size(N,2)==3
   P = N(1); nu = N(2); M = N(3);
elseif size(N,2)==2
   P = N(1); nu =N(2); M = 1;
elseif size(N,2)==1
   nu = 1; P = N; M = 1;
else
   ctrlMsgUtils.error('Ident:dataprocess:idinput1')
end

if isempty(nosine),nosine=[10,10,1];end
nosine = nosine(:).';
if length(nosine)==1
   nosine = [nosine,10,1];
elseif length(nosine)==2
   nosine = [nosine,1];
end

if isempty(levels),levels=[-1,1];end
if isempty(band),band=[0 1];end
if isempty(type),type='rbs';end

if band(2)<band(1) && ~strcmpi(type,'prbs')
   ctrlMsgUtils.error('Ident:dataprocess:idinput2')
end

if levels(2)<levels(1)
   ctrlMsgUtils.error('Ident:dataprocess:idinput3')
end

if nosine(1)<1 || nosine(2)<1
   ctrlMsgUtils.error('Ident:dataprocess:idinput4')
end

if ~ischar(type)
   ctrlMsgUtils.error('Ident:dataprocess:idinput5')
end

if strcmpi(type,'rs') || strcmpi(type,'rgs')
   u=randn(5*P,nu);
   if ~all(band==[0 1])
      u = idfilt(u,8,band);
   end
   u = u(2*P+1:end-2*P,:); % to take out transients
   %u = u - ones(P,1)*mean(u)+(levels(2)+levels(1))/2;
   for ku = 1:nu
      u(:,ku) = u(:,ku)-mean(u(:,ku));%Now it is zero mean.
      % The standdard deviation shall be (lev(2)-lev(1))/2
      % and the mean shall be (lev(2)+lev(1))/2
      u(:,ku)=(levels(2)+levels(1))/2 + u(:,ku)*(levels(2)-levels(1))/2/...
         sqrt(u(:,ku)'*u(:,ku)/length(u(:,ku)));
      %u(:,ku) = u(:,ku)/norm(u(:,ku))*sqrt(P)*(levels(2)-levels(1))/2;
   end
   
elseif strcmpi(type,'rbs')
   u=randn(5*P,nu);
   if ~all(band==[0 1]),u = idfilt(u,8,band);end
   u = sign(u(2*P+1:end-2*P,:)); % to take out transients
   u = (levels(2)-levels(1))*(u+1)/2+levels(1);
   
elseif strcmpi(type,'prbs')
   u = controllib.internal.util.genPRBS(P,nu,M,band,levels);
   
elseif strcmpi(type,'sine')
   odd = nosine(3);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%   !! HERE I MADE SOME CHANGES !!   %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   test_tol = 1e-8;
   possfreq = 2*pi*(0:odd:fix(P/2)-1)/P;
   possfreq = possfreq(band(1)*pi - possfreq <=  test_tol & ...
                       band(2)*pi - possfreq >= -test_tol);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%           END OF CHANGES           %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   nl = length(possfreq);
   if nl<nu*nosine(1)
      ctrlMsgUtils.error('Ident:dataprocess:idinput9',nosine(1))
   end
   fnr = choose(nl,nu*nosine(1));
   numtrial = nosine(2);
   for ku = 1:nu
      freqs = possfreq(fnr(ku:nu:end));
      freq(ku,:) = freqs;
      for kn = 1:numtrial
         ut = zeros(P,1);
         for ks = 1:nosine(1)
            ut = ut+cos((0:P-1)'*freqs(ks)+rand(1,1)*2*pi);
         end
         mm = max(ut); mn = min(ut);
         if kn==1, u1 = ut; bestamp = mm-mn; end
         if mm-mn<bestamp, u1 = ut; bestamp = mm-mn;end
      end
      mn = min(u1); mm = max(u1); mm = max(mm,-mn); mn = -mm;
      u(:,ku) = (levels(2)-levels(1))*(u1-mn)/(mm-mn)+levels(1);
   end
else
   ctrlMsgUtils.error('Ident:dataprocess:idinput5')
end
if M>1
   uu = u;
   for km = 2:M
      u = [uu;u];
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fnr = choose(nr,nl)
%Choose nl values from integers 1:nr as evenly spread as possible
nnr = 1:nr;
s = 1;
fnr = zeros(1,nl);
while nl>0
   nr = length(nnr);
   if nl == 3
      fnr(s) = nnr(1);
      fnr(s+1)= nnr(ceil(nr/2));
      fnr(s+2) = nnr(end);
      nl = 0;
   elseif nl>1
      k = floor(nr/(nl-1));
      fnr(s)=nnr(1);
      fnr(s+1)=nnr(end);
      s=s+2;
      nnr=nnr(k+1:end-k);
      nl = nl-2;
   else
      k=ceil(nr/2);
      fnr(s) = nnr(k);
      s=s+1;
      nl = nl -1;
   end
end
fnr = sort(fnr);

%&&&&&&&&&&&&
%  u = iddata([],u);
%  if M>1
%      set(u,'Period',P*ones(nu,1))
%  end
%
