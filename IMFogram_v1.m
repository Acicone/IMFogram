function varargout = IMFogram_v1(IMF,fs,winOverPerc,IFint,NIF,winLen,verbose,scale,fig)

% IMFogram_v1(IMF,fs,winOverPerc,IFint,NIF,winLen,fig)
%
% Produces the time-frequency plots of the Absolute values of the IMFs Amplitudes
% produced by (Fast) Iterative Filtering
%
%  Inputs
%
%  IMFs        = signal to be extended. Each row is a different IMF.
%                This matrix can include the trend which will not be
%                considered in the calculations
%  fs          = sampling frequency of the original signal
%  winOverPerc = Percentage value of overlapping of the sliding windows
%                Value between 0 and 100. Default value (50)
%  IFint       = Lower and Upper bound of the frequency interval.
%                Default value = the average frequencies of the first and
%                last IMFs
%  NIF         = Number of rows of the IMFogram, i.e. number of
%                Instantaneous Frequencies values to be plotted in the
%                interval IFint. Default value (2000)
%  winLen      = Number of the sample points of the time window on which
%                computing the aewrage Frequency and Amplitude.
%                Default value = length of the signal / 200
%  verbose     = boolean which allows to specify if we want the algorithm
%                to give use statistics about the IMFs and their extrema
%                Default value set to (true)
%  scale       = frequency scale. Options
%                'linear'
%                'log' log along the frequency axis
%                'loglog' log along both the frequency and the amplitude
%                         values axes
%                'mel'
%                Default value set to ('linear')
%  fig         = handling of the figure to be generated.
%                If set to 0 no figure will be displayed.
%
%  Outputs
%
%  A       = NxM matrix containing the Amplitudes values. N = NIF and
%            represent the number of Instantaneous Frequencies values.
%            M is the number of time windows on which we devide the
%            original signal
%  IF      = Instantaneous Frequencies values for each row of A
%  IMF_Num = NxM matrix containing the IMF number(s) corresponding to each
%            row and column of A
%  IMF_iF  = Instantaneous Frequencies function associated with each IMF
%  IMF_iA  = Instantaneous Amplitudes function associated with each IMF
%  figA    = handle of the output figure to be generated
%  X       = NxM matrix containing the time
%  Y       = NxM matrix containing the frequency
%
%  Five options: None required -> IMFogram_v1 outputs the plot
%                4 required -> IMFogram_v1 outputs X,Y,A and IF
%                5 required -> IMFogram_v1 outputs X,Y,A, IF, the plot and its
%                handle
%                6 required -> IMFogram_v1 outputs X,Y,A, IF, IMF_Num, the plot
%                and its handle
%                8 required -> IMFogram_v1 outputs X, Y, A, IF, IMF_Num, IMF_iF,
%                IMF_iA, the plot and its handle
%
%   See also IF2mel, hz2mel, mel2hz

% we use overlapping windows


if nargin<1 help IMFogram_v1; return; end

if nargin<1 help IMFogram_v1; return; end
if nargin<3 || isempty(winOverPerc), winOverPerc = 50; end % the windows will overlap 50%
if nargin<7 || isempty(verbose), verbose = true; end
if nargin<8 || isempty(scale), scale = 'linear'; end
if nargin < 9, fig=[]; end
figA=0;

[M0,N0]=size(IMF);

if nargin<5 || isempty(NIF), NIF=2000; end

zerocrossings=cell(1,M0);
IMF_iA=cell(1,M0);
IMF_iF=cell(1,M0);

min_iF=zeros(1,M0);
max_iF=zeros(1,M0);
% function for the zero crossing identification
zci = @(v) find(diff(sign(v)));

for i = 1:M0
    zerocrossings{i} = zci(IMF(i,:));
    
    
    
    if length(zerocrossings{i})<2
        IMF_iA{i}=zeros(1,N0);
        IMF_iF{i}=zeros(1,N0);
    else
        temp_val = fs./(2*diff(zerocrossings{i}));
        max_iF(i)=max(temp_val);
        min_iF(i)=min(temp_val);
        temp_IMF_maxval=zeros(1,length(zerocrossings{i})-1);
        temp_IMF_maxval_pos=zeros(1,length(zerocrossings{i})-1);
        for j=1:length(zerocrossings{i})-1
            [temp_IMF_maxval(j),temp_pos]=max(abs(IMF(i,zerocrossings{i}(j):zerocrossings{i}(j+1))));
            temp_IMF_maxval_pos(j)=temp_pos+zerocrossings{i}(j)-1;
        end
        % we prepare to comute the iAmplitude curve
        % we remove all repeated entries from temp_IMF_maxval_pos and the
        % corresponding entries in temp_IMF_maxval
        [IMF_maxval_pos,IMF_maxval_pos_pos]=unique([1 temp_IMF_maxval_pos N0],'first');
        IMF_maxval=[temp_IMF_maxval(1) temp_IMF_maxval temp_IMF_maxval(end)];
        IMF_maxval=IMF_maxval(IMF_maxval_pos_pos);
        IMF_iA{i}=interp1(IMF_maxval_pos,IMF_maxval,1:N0,'linear');
        if zerocrossings{i}(1)==1 && zerocrossings{i}(end)==N0
            IMF_iF{i}=interp1([zerocrossings{i}],[temp_val temp_val(end)],1:N0,'linear');
        elseif zerocrossings{i}(1)~=1 && zerocrossings{i}(end)~=N0
            IMF_iF{i}=interp1([1 zerocrossings{i} N0],[temp_val(1) temp_val temp_val(end) temp_val(end)],1:N0,'linear');
        elseif zerocrossings{i}(1)~=1 && zerocrossings{i}(end)==N0
            IMF_iF{i}=interp1([1 zerocrossings{i}],[temp_val(1) temp_val temp_val(end)],1:N0,'linear');
        elseif zerocrossings{i}(1)==1 && zerocrossings{i}(end)~=N0
            IMF_iF{i}=interp1([zerocrossings{i} N0],[temp_val temp_val(end) temp_val(end)],1:N0,'linear');
        end
    end
end

if nargin < 4 || isempty(IFint)
    % we compute the highest average frequency
    IFint(2)=max(max_iF);
    % and the lowest average frequency
    IFint(1)=min(min_iF);
end
if nargin < 6 || isempty(winLen)
    winLen=floor(N0/200); % we define the smallest time window on which we compute the average amplitudes and frequencies
end
DeltaWin = floor(winLen*(1-winOverPerc/100)); % Nonoverlapping Window length
if DeltaWin==0
    error(' DeltaWin = 0! Try increasing ''winLen'' or decreasing ''winOverPerc''')
end
Nwin=floor((N0-(winLen-DeltaWin))/DeltaWin);



if strcmp(scale,'mel')
    IF = IF2mel(IFint(1), IFint(2), NIF, 0);
elseif strcmp(scale,'linear')
    IF=linspace(IFint(1),IFint(2),NIF); % Instantaneous Frequency values used to produce the IMFogram
elseif strcmp(scale,'log') || strcmp(scale,'loglog')
    if IFint(1)>0
        IF=logspace(log10(IFint(1)),log10(IFint(2)),NIF); % Instantaneous Frequency values used to produce the IMFogram
    else
        IF=logspace(log10(IFint(1)+eps),log10(IFint(2)),NIF);
    end
end

N_IF=length(IF);

A=zeros(N_IF,Nwin); % we prepare the Matrix used to plot the IMFogram
if nargout>3
    IMF_Num=cell(N_IF,Nwin); % we prepare the Matrix containing the IMF numbers
end
for ii=1:M0 % we scan all IMFs containing more than 1 extrema
    for jj=1:Nwin % we study the IMFogram over windows in time
        temp_val=sum(IMF_iF{ii}(1+(jj-1)*DeltaWin:winLen+(jj-1)*DeltaWin))/winLen;
        if temp_val<IFint(2)+(IF(2)-IF(1)) && temp_val>IFint(1)-(IF(2)-IF(1))
            [v,pos]=min(abs(temp_val-IF)); % with mink we can spread the outcome to nearby frequencies to make more readable the plot % ceil(5*N_IF/1000) %originally: ceil(1.5*N_IF/100)
            
            if nargout>3
                IMF_Num{pos,jj}=[IMF_Num{pos,jj} ii];
            end
            A(pos,jj)=A(pos,jj)+sum(IMF_iA{ii}(1+(jj-1)*DeltaWin:winLen+(jj-1)*DeltaWin))/winLen; % we consider the average amplitude to mimic what is done in the spectrogram where we compute the amplitude of the stationary sin or cos
        end
    end
end

if nargout<1
    if nargin < 9 || isempty(fig)
        figure
    else
        figure(fig)
    end
    [X,Y] = meshgrid((1:Nwin)*DeltaWin/fs,IF);
    surf(X,Y,A,'edgecolor','none')
    %     N_colors=1000;
    %     colormap([1 1 1;
    %         ones(N_colors,1) (1./1000.^(1/N_colors:1/N_colors:1))' zeros(N_colors,1)])
    shading interp
    oldcolormap=colormap('parula');
    colormap([1 1 1; flipud(oldcolormap) ]);
    title('IMFogram')
    view(0,90)
    %yticklabels(xlables)
    
elseif nargout==4
    [X,Y] = meshgrid((1:Nwin)*DeltaWin/fs,IF);
    varargout{1}=X;
    varargout{2}=Y;
    varargout{3}=A;
    varargout{4}=IF;
elseif nargout==5
    [X,Y] = meshgrid((1:Nwin)*DeltaWin/fs,IF);
    if isempty(fig) || not(fig==0)
        if nargin < 9 || isempty(fig)
            figA=figure;
        else
            figA=figure(fig);
        end
        surf(X,Y,A,'edgecolor','none')
        %     N_colors=1000;
        %     colormap([1 1 1;
        %         ones(N_colors,1) (1./1000.^(1/N_colors:1/N_colors:1))' zeros(N_colors,1)])
        shading interp
        oldcolormap=colormap('parula');
        colormap([1 1 1; flipud(oldcolormap) ]);
        title('IMFogram')
        view(0,90)
        %yticklabels(xlables)
    end
    varargout{1}=X;
    varargout{2}=Y;
    varargout{3}=A;
    varargout{4}=IF;
    varargout{5}=figA;
elseif nargout==6
    [X,Y] = meshgrid((1:Nwin)*DeltaWin/fs,IF);
    if isempty(fig) || not(fig==0)
        if nargin < 9 || isempty(fig)
            figA=figure;
        else
            figA=figure(fig);
        end
        surf(X,Y,A,'edgecolor','none')
        %     N_colors=1000;
        %     colormap([1 1 1;
        %         ones(N_colors,1) (1./1000.^(1/N_colors:1/N_colors:1))' zeros(N_colors,1)])
        shading interp
        oldcolormap=colormap('parula');
        colormap([1 1 1; flipud(oldcolormap) ]);
        title('IMFogram')
        view(0,90)
        %yticklabels(xlables)
    end
    varargout{1}=X;
    varargout{2}=Y;
    varargout{3}=A;
    varargout{4}=IF;
    varargout{5}=IMF_Num;
    varargout{6}=figA;
elseif nargout==8
    [X,Y] = meshgrid((1:Nwin)*DeltaWin/fs,IF);
    if isempty(fig) || not(fig==0)
        if nargin < 9 || isempty(fig)
            figA=figure;
        else
            figA=figure(fig);
        end
        
        surf(X,Y,A,'edgecolor','none')
        %     N_colors=1000;
        %     colormap([1 1 1;
        %         ones(N_colors,1) (1./1000.^(1/N_colors:1/N_colors:1))' zeros(N_colors,1)])
        shading interp
        oldcolormap=colormap('parula');
        colormap([1 1 1; flipud(oldcolormap) ]);
        title('IMFogram')
        view(0,90)
        %yticklabels(xlables)
    end
    varargout{1}=X;
    varargout{2}=Y;
    varargout{3}=A;
    varargout{4}=IF;
    varargout{5}=IMF_Num;
    varargout{6}=IMF_iF;
    varargout{7}=IMF_iA;
    varargout{8}=figA;
end

end

%% Auxiliar functions

function melIF = IF2mel(min_iF, max_iF, nfilts, htkmel)



minmel = hz2mel(min_iF, htkmel);
maxmel = hz2mel(max_iF, htkmel);
melIF = mel2hz(minmel+[0:(nfilts-1)]/(nfilts-1)*(maxmel-minmel), htkmel);

end

function f = mel2hz(z, htk)
%   f = mel2hz(z, htk)
%   Convert 'mel scale' frequencies into Hz
%   Optional htk = 1 means use the HTK formula
%   else use the formula from Slaney's mfcc.m
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
    htk = 0;
end

if htk == 1
    f = 700*(10.^(z/2595)-1);
else
    
    f_0 = 0; % 133.33333;
    f_sp = 200/3; % 66.66667;
    brkfrq = 1000;
    brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
    logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)
    
    linpts = (z < brkpt);
    
    f = 0*z;
    
    % fill in parts separately
    f(linpts) = f_0 + f_sp*z(linpts);
    f(~linpts) = brkfrq*exp(log(logstep)*(z(~linpts)-brkpt));
    
end

end

function z = hz2mel(f,htk)
%  z = hz2mel(f,htk)
%  Convert frequencies f (in Hz) to mel 'scale'.
%  Optional htk = 1 uses the mel axis defined in the HTKBook
%  otherwise use Slaney's formula
% 2005-04-19 dpwe@ee.columbia.edu

if nargin < 2
    htk = 0;
end

if htk == 1
    z = 2595 * log10(1+f/700);
else
    % Mel fn to match Slaney's Auditory Toolbox mfcc.m
    
    f_0 = 0; % 133.33333;
    f_sp = 200/3; % 66.66667;
    brkfrq = 1000;
    brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
    logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)
    
    linpts = (f < brkfrq);
    
    z = 0*f;
    
    % fill in parts separately
    z(linpts) = (f(linpts) - f_0)/f_sp;
    z(~linpts) = brkpt+(log(f(~linpts)/brkfrq))./log(logstep);
    
end

end
