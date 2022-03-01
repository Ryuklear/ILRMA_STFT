function [mix,audio,Param] = audio_rc(audio,Param,STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����̓ǂݍ��݂ƃC���p���X�����̏�ݍ��݂��s���܂��D
%STFT�͂����ł͍s���܂���D�����܂ō���(�ϑ�)�M���̍쐬�D
%
%[Inputs]
%       audio: �e�����Ɋւ���������\����
%       Param: main�Œ�`�����e�p�����[�^�̍\����
%        STFT: STFT�Ɋւ���p�����[�^�̍\����
%[Outputs]
%         mix:��ݍ��݌�̍����M�� (mic x justsamples)
%       audio:�����������\����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = Param.ns;
M = Param.mic;
fsResample = Param.fsResample;

%Input data for all sound source 
for i = 1:N
   [audio(i).sig, fs] = audioread(audio(i).name);
end

%Resample
for i = 1:N
   audio(i).resample = resample(audio(i).sig, fsResample, fs); 
end

%just_sample = fix((size(audio(1).resample,1)-fftSize)/shiftSize)*shiftSize+fftSize;
Param.sig_length = fix((fsResample*Param.seconds-STFT.fftSize)/STFT.shiftSize)*STFT.shiftSize+STFT.fftSize;

%Ajustment of resample length  
%Changing monoral signal
for i = 1:N
   audio(i).resample = audio(i).resample(1:Param.sig_length,1)';
end

mix = zeros(M,size(audio(1).resample,2));
for n = 1:N
audio(n).music = zeros(M,size(audio(1).resample,2));
end

%Convolution
for n = 1:N
  for m = 1:M
      music = conv(audio(n).resample,audio(n).impulse(m,:));
      audio(n).music(m,:) = music(1:length(audio(n).resample));
  end
end

%Make mixed signal
for m = 1:M
    for n = 1:N
      mix(m,:) = mix(m,:) + audio(n).music(m,:);
    end
end

end
%%%%%%%%%%EOF%%%%%%%%%%