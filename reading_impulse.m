function [audio] = reading_impulse(audio,impulse_name,Param,fftSize)
%RWCP���������f�[�^�x�[�X�@E2A�C���p���X�����f�[�^
%�f�[�^�@48kHz 32bit(float) raw�t�@�C��(�o�C�i���f�[�^)
%E2A�͎c������0.3s�̎c���ώ��Ŏ��^���ꂽ�f�[�^
%�������~�}�C�N���̃C���p���X�������g�p
Fs = 48000;

for n = 1:Param.ns
for m = 1:Param.mic
   fileID = fopen(impulse_name((n-1)*Param.mic+m,:));    
   h = fread(fileID,'float');
   h = resample(h,Param.fsResample,Fs);
   fclose(fileID);
   audio(n).impulse(m,:) = h(1:fftSize);
end
end

%�ǂݍ��񂾃f�[�^�͍\����audio.impulse�Ƃ��ĕԂ����

end
%%%%%EOF%%%%%