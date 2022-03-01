function[B,r] = rankr_learning(filename,Param,STFT)
%�y��̒P�������s��Ƃ��Ē��o(�قژ_���ʂ�) 

fs = Param.fsResample;
L = STFT.fftSize;
eta = STFT.shiftSize;

[music1,Fs] = audioread(filename);
music1 = resample(music1,fs,Fs);
T1 = size(music1,1);  
T2 = fix((T1 - L + eta)/eta) * eta + L;
music1 = music1(1:T2);

switch STFT.type
    case 1
    [X_comp,~] = STFT(music1',L,eta,STFT.w_name); 
    case 2
     X_comp = STFT_basic(music1',STFT.analysis_window, eta, STFT.frag);   
end
X_mono = abs(X_comp).^2;   %��X�̌����ł̓p���[�X�y�N�g��!
J = size(X_mono,2); %�s�̒����i���j

%SVD���{��
[U,S,V] = svd(X_mono); 
%diag()�ł͑Ίp�������擾�i�����ł͓��ْl�擾�j
s = diag(S);
%�����N�������߂邽�߂̌덷�͈�(10%�ȉ�)
epsilon = norm(X_mono,'fro')/10;

%�����N����(B�x�N�g����)������
for r = 1:(length(s)-1)
    %���{�ڂŌ덷�͈͂𖞂�����
    %[.]�͊e�v�f���Ƃ̈Ӗ�
    if sqrt(sum(s(r+1:end).^2)) <= epsilon   %�t���x���j�E�X�m����?
        break
    end
end

%���s��B�̐݌v
B = zeros(size(X_mono,1),r);
%B��1��ڂ�U��1���
B(:,1) = abs(U(:,1));

%X_mono�̃����Nr�ߎ�
X_mono_r = U(:,1:r)*S(1:r,1:r)*(V(:,1:r)');

%B��2��ڈȍ~��X_mono_r����擾����
for k = 2:r
    PB_x = zeros(size(X_mono,1),J);
    Y = zeros(1,J);
    
   for i = 1:J
       PB_x(:,i) = B(:,1:k-1)*((B(:,1:k-1)'*B(:,1:k-1))\(B(:,1:k-1)'*X_mono_r(:,i)));
       Y(1,i) = abs(pi/2 - acos(PB_x(:,i)'*X_mono_r(:,i)/(norm(PB_x(:,i))*norm(X_mono_r(:,i)))));
   end
    [~, I] = min(Y);
    B(:,k)= X_mono_r(:,I)/norm(X_mono_r(:,I));
end

%B�̕��̒l��0��
B(B<0)=0;
%norm(B,'fro')^2

for k=2:r
    B(:,k) = B(:,k)/norm(B(:,k));
end