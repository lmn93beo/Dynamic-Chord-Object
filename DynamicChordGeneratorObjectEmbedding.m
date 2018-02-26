%%
% DRC mark 2
if exist('Ptotal')
    Ptotal = P*0;
end
includeObject = 1; % toggle for whether there is an object or not;
for i = 1%:20
duration = 1;

sampleRate = 100000;
toneLength = 0.025;
window = 0.025;
range = 1000*(4:0.125:46);% possibly make logspace
maxTones = 15;
toneSegSamps = (1/sampleRate):(1/sampleRate):(toneLength);
toneL = length(toneSegSamps);

soundfield = zeros(1,sampleRate*duration);

currNum = 0;
ObjectseedF = (randi(length(range),1,1));%probably need to make this draw from an octave space distribution (e.g. with logspace function)
maxObjectStep = 5;% maximum size of step for the object related tones 


for n = 1:(-1+duration/window)
    numTones(n) = randi(maxTones);%-currNum);
    tonePos = randi(length(range),numTones(n),1); %probably need to make this draw from an octave space distribution (e.g. with logspace function)
    intensities = random('Normal',1,0.1,numTones(n)+1,1); %round(length(range)/2)+randi(round(length(range)/2),numTones(n),1);
    intensities = 0.9*intensities./sum(intensities);
    winPoints = window*sampleRate;
    divVal = 0.9;
    for m = 1:(numTones(n))
        soundfield(((n-1)*winPoints+1):((n-1)*winPoints+toneL)) = soundfield(((n-1)*winPoints+1):((n-1)*winPoints+toneL))+intensities(m).*sin(range(tonePos(m))*toneSegSamps*2*pi);
        range(tonePos(m));
    end
    if n > 1
        currNum = numTones(n)-1;
    end
    
    % embed object
    if includeObject == 1
        if n == 1
            currObjectPoint = ObjectseedF;
        else
            if currObjectPoint+maxObjectStep < length(range)
                currObjectPoint = (currObjectPoint+round(randi(maxObjectStep)-maxObjectStep/2));
            elseif currObjectPoint-maxObjectStep > 1% to keep the object from going out of range
                currObjectPoint = (currObjectPoint-round(randi(round(maxObjectStep/2))));
            else
                currObjectPoint = (currObjectPoint+round(randi(round(maxObjectStep/2))));
            end
        end
        currObjectFreq = range(currObjectPoint);
        soundfield(((n-1)*winPoints+1):((n-1)*winPoints+toneL)) = soundfield(((n-1)*winPoints+1):((n-1)*winPoints+toneL))+mean(intensities).*sin(currObjectFreq*toneSegSamps*2*pi);
    end
    
end

soundfieldHold = soundfield;
for u = 1 % I think addition of non-tonal background noise is not necessary for this task
    soundfield = soundfieldHold;
    if u == 1
        tag = 'none'
        lowValCut = 0;
        highValCut = 0;
    elseif u == 3
        tag = 'low'
        lowValCut = 4000;
        highValCut = 25000;
    elseif u == 2
        tag = 'high'
        lowValCut = 25000;
        highValCut = 46000;
    end
if u > 1
    noise = wgn(sampleRate*duration,1,0);

    hf = design(fdesign.bandpass('N,F3dB1,F3dB2',10,lowValCut,highValCut,sampleRate));

    noise = filter(hf,noise)/2;
    sound(noise,100000)

    soundfield = soundfield+noise';
else
    noise = wgn(sampleRate*duration,1,0);
    noise = noise/200
    soundfield = soundfield+noise';
end



soundfield = 0.9*soundfield/max([abs(min(soundfield)) max(soundfield)]); %re-range for sound export
figure
%subplot(2,1,1);spectrogram(soundfield,500,400,[],sampleRate,'yaxis')
spectrogram(soundfield,500,400,[],sampleRate,'yaxis')
[S,F,T,P] = spectrogram(soundfield,500,400,[],sampleRate,'yaxis');

soundStart = soundfield(1:(length(soundfield)-round(winPoints/2)))
soundOffset = soundfield((round(winPoints/2)+1):end)

% [Cxy,W] = mscohere(soundStart',soundOffset',window*sampleRate/2,window*sampleRate/4,[],sampleRate)
% CxySmo = nanconv(Cxy,gausswin(20*maxObjectStep)./sum(gausswin(20*maxObjectStep)))
% subplot(2,1,2);plot(CxySmo,W);ylim([min(range) max(range)])
if ~exist('Ptotal')
    Ptotal = P*0;
end
Ptotal = Ptotal+P;
eval(['DRC' num2str(n) tag ' = soundfield;'])

%audiowrite(['dynamicRandChord' num2str(i) tag '.wav'],soundfield,100000,'BitsPerSample',16)
end
end
%imagesc(Ptotal)
%coherence assessment