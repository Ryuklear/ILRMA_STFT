function[B,r] = rank1_learning(filename,Param,STFT)
%�y��̒P�������s��Ƃ��Ē��o
r = 1;         %�����N�P 
fs = Param.fsResample;
L = STFT.fftSize;
eta = STFT.shiftSize;

[music1,Fs] = audioread(filename);
% T = size(music1,1);  Fs���̃T���v�����i���Ԃ�T/Fs[s]�j

music1 = resample(music1,fs,Fs);
T1 = size(music1,1);   %music1�̑��T���v����
%music1��STFT���傤�ǂ̒���
T2 = fix((T1 - L + eta)/eta) * eta + L;

music1 = music1(1:T2);

switch STFT.type
    case 1
    [X_comp,~] = STFT(music1',L,eta,STFT.w_name); %���f���l�����܂�
    case 2
     X_comp = STFT_basic(music1',STFT.analysis_window, eta, STFT.frag);
end
X_mono = abs(X_comp).^2;   %��X�̌����ł̓p���[�X�y�N�g��!

I = size(X_mono,1); %�s��̗�̒����i�c�j
J = size(X_mono,2); %�s�̒����i���j

%SVD���{��
[U,~,~] = svd(X_mono); 

%���s��B�̐݌v
B = zeros(size(X_mono,1),r);
%B��1��ڂ�U��1���
B(:,1) = abs(U(:,1));

%B�̕��̒l��0��
B(B<0)=0;
%norm(B,'fro')^2

    B(:,r) = B(:,r)/norm(B(:,r));
end