function[Y,W,A,T,V] = ILRMA1_initial(X_tilde,T_0,V_0,Param)
%ILRMA1_initial�Ƃ����֐����쐬�A���͂�3�ϐ�(X_tilde,T_0,V_0)
%�o�͂�5�ϐ�(Y,W,A,T,V)��Ԃ�
%I = freqention bin, J = time frames
%input

N = Param.ns;

%X_tilde�̃T�C�Y��I,J,N�Ɋi�[����B
%for�����񂷂���
[I,J,M] = size(X_tilde);
X = zeros(J,M,I);
for i=1:I
    for j=1:J
        for m=1:M
            %�]�u�H���Ȃɂ����Ă���H
            %j,n,i�̏��ɒ���
            X(j,m,i)=X_tilde(i,j,m);
        end
    end
end
    
%Initialize W(demixing matrix) with identity matrix or complex-valued random values 
% and T and V with nonnegative random values(1);
%�P�ʍs��(eye)���쐬
%W = eye(N,M);
W = Param.W_opt;
A = eye(M,N);
%initialize W
for i = 1:I
    %W(:,:,i) = W(:,:,1); %initialize W
    A(:,:,i) = A(:,:,1);
end

%A = eye(N,N);
%initialize A
% for i = 1:I
%     A(:,:,i) = A(:,:,1) ; %initialize W
% end

%initialize T & V
T = T_0;%rand(I,L,N);
V = V_0;%rand(L,J,N);

Y=zeros(I,J,N);
X_tilde2 = zeros(M,J,I);

for i=1:I
    for n=1:N
        for m = 1:M
            X_tilde2(m,:,i) = X_tilde(i,:,m);
        end
        Y(i,:,n) = W(n,:,i) * X_tilde2(:,:,i);
    end
end

P = Y.*conj(Y);% P=zeros(I,J,2);

R = zeros(I,J,N);
for n=1:N
    R(:,:,n) = T(:,:,n) * V(:,:,n);
end

U = zeros(M,M,I,N);

fprintf('ILRMA Iteration:   ');
for KKK=1:Param.it
    fprintf('\b\b\b\b%4d', KKK);
    for n=1:N
        %Update of basis matrix
        if Param.supervised == false
        T(:,:,n) = max( T(:,:,n) .* sqrt( ( (P(:,:,n) .* (R(:,:,n).^(-2))) * (V(:,:,n)') ) ./ ( (R(:,:,n).^(-1)) * (V(:,:,n)') ) ) , eps);
        R(:,:,n) = T(:,:,n) * V(:,:,n);
        end
        V(:,:,n) = max( V(:,:,n) .* sqrt( ( (T(:,:,n)') * (P(:,:,n) .* (R(:,:,n).^(-2))) ) ./ ( (T(:,:,n)') * (R(:,:,n).^(-1)) ) ) , eps);
        R(:,:,n) = T(:,:,n) * V(:,:,n);
    end
        
    for i=1:I
        %A(:,:,i) = inv(W(:,:,i)); 
        %W_update = zeros(N,N);
        for n=1:N
            U(:,:,i,n) = (1/J) * ( (X(:,:,i)') * ( X(:,:,i) .* ((R(i,:,n).^(-1))' * ones(1,M)) ) ).'; 
%             E = eye(N);
%             en = E(:,n); 
%             wni = (W(:,:,i) * U(:,:,i,n)) \ en;
            wni = U(:,:,i,n) \ A(:,n,i);
            %W_update(n,:) = wni' / sqrt(wni' * U(:,:,i,n) * wni);
            W(n,:,i) = wni' / sqrt(wni' * U(:,:,i,n) * wni);
        end
        %W(:,:,i) = W_update;
    end
        
    for i=1:I
        for n=1:N
            Y(i,:,n) = W(n,:,i) * X_tilde2(:,:,i);
        end
    end
    
    P = Y.*conj(Y);
    
    for i=1:I
        A(:,:,i) = pinv(W(:,:,i));
    end    
    
    PA = A.*conj(A);
    for n = 1:N
        lambda = sqrt( sum(sum(squeeze(PA(:,n,:)))) /I);
        A(:,n,:) = A(:,n,:)/lambda;
        W(n,:,:) = W(n,:,:)*lambda;
        P(:,:,n) = P(:,:,n)*(lambda^2);
        R(:,:,n) = R(:,:,n)*(lambda^2);
        if Param.supervised == false
        T(:,:,n) = T(:,:,n)*(lambda^2);
        end
        V(:,:,n) = V(:,:,n)*(lambda^2);
    end
   
%     for n = 1:N
%         lambda = sqrt( sum(sum(P(:,:,n))) / (I*J) );
%         
%         for i = 1:I
%             W(n,:,i) = W(n,:,i)/lambda;
%         end
%         P(:,:,n) = P(:,:,n)/(lambda^2);
%         R(:,:,n) = R(:,:,n)/(lambda^2);
%         T(:,:,n) = T(:,:,n)/(lambda^2);
%     end
end

for i=1:I
    for n=1:N
        Y(i,:,n) = W(n,:,i) * X_tilde2(:,:,i);
    end
end
fprintf(' ILRMA1_initial done.\n');
end
