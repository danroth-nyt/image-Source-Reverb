%Reverb calculator using image/source modelling
%parameters
Lx = 40; %length of room, m
Ly = 38; %width of room, m
Lz = 23; %height of room, m
a = 15; %listener x location, m
b = 21; %listener y location, m
c = 3; %listener z location, m
p = 30; %source x location, m
q = 22; %source y location, m
r = 12; %source z location, m
alpha = 0.6; % room reflection coefficient, 0-1
cs = 343; %speed of sound m/s
Fs = 48000; %sample rate
N = 10; %max image source order
dE = 0.2; %distance between ears

%create Impulse Response(IR) vector
simtime = (1/cs) * sqrt((N*Lx+p-a)^2 + (N*Ly+q-b)^2 + (N*Lz+r-c)^2);
IRvector = zeros(round(simtime * Fs),2);

%listener positions
angle = atan((p-a)/(q-b));
aRight = a - ((dE/2)*cos(angle));
aLeft = a + ((dE/2)*cos(angle));
bRight = b - ((dE/2)*cos(angle));
bLeft = b + ((dE/2)*cos(angle));

%calculate IR coefficients
for d = -N : N %x dimension
    if mod(N,2) ~= 0
        Aleft = (d+1)*Lx - p - aLeft;
    else
        Aleft = d*Lx + p - aLeft;
    end
    if mod(N,2) ~= 0
        Aright = (d+1)*Lx - p - aRight;
    else
        Aright = d*Lx + p - aRight;
    end
    [dwall1, dwall2] = wallHits(d);
    for e = -N : N %y dimension
        if mod(N,2) ~= 0
            Bleft = (e+1)*Ly - q - bLeft;
        else
            Bleft = e*Ly + q - bLeft;
        end
        if mod(N,2) ~= 0
            Bright = (e+1)*Ly - q - bRight;
        else
            Bright = e*Ly + q - bRight;
        end
        [ewall1, ewall2] = wallHits(e);
            for f = -N : N %z dimension
                if mod(N,2) ~= 0
                    C = (f+1)*Lz - r - c;
                else
                    C = f*Lz + r - c;
                end
                [fwall1, fwall2] = wallHits(f);
                w = abs(dwall1 + dwall2) + abs(ewall1 + ewall2) + abs(fwall1 + fwall2);
                Lleft = sqrt(Aleft^2 + Bleft^2 + C^2);
                tLeft = Lleft/cs;
                gLeft = (1/Lleft) * alpha^w;
                Lright = sqrt(Aright^2 + Bright^2 + C^2);
                tRight = Lright/cs;
                gRight = (1/Lright) * alpha^w;
                IRvector(round(tLeft*Fs),1) = IRvector(round(tLeft*Fs)) + gLeft;
                IRvector(round(tRight*Fs),2) = IRvector(round(tRight*Fs)) + gRight;
            end
    end
end

%HF absorption, 6th Order Butterworth 8kHz Cutoff
[B1, A1] = butter(6, 8000/(Fs/2), 'low');
IRvector = filter(B1,A1,IRvector);

%load and process audio
[dry, Fs] = audioread('dry.wav');
wet(:,1) = conv(dry,IRvector(:,1));
wet(:,2) = conv(dry,IRvector(:,2));
wet = wet / max(abs(wet));

%play audio
player1 = audioplayer(dry,Fs);
playblocking(player1);
player2 = audioplayer(wet,Fs);
play(player2);

