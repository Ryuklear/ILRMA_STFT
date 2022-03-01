function[Y,A] = ILRMA1_Sparse(X_tilde,T_0,V_0,Param,Sparse,L)
%I = freqention bin, J = time frames
%input
%X_tilde(I,J,M):observed multichannel omplex-valued signals I*J*M
%T_0:basis matrix for the nth sourse
%V_0:activation matrix for the nth source

N = Param.ns;

[I,J,M] = size(X_tilde);
X=zeros(J,M,I);
for i=1:I
    for j=1:J
        for m=1:M
            X(j,m,i)=X_tilde(i,j,m);
        end
    end
end
    
%Initialize W(demixing matrix) with identity matrix or complex-valued random values 
% and T and V with nonnegative random values(1);
%W = eye(N,M);
W = Param.W_opt;
A = eye(M,N);
%initialize W
for i = 1:I
    %W(:,:,i) = W(:,:,1);
    A(:,:,i) = A(:,:,1);
end

%initialize T & V
T = T_0;%rand(I,Kn,N);
V = V_0;%rand(Kn,J,N);

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

W_hat = zeros(N,M,I);
A_hat = zeros(M,N,I);
U = zeros(M,M,I,N);
U_hat = zeros(M,M,I,N);
h_all = zeros(2*I-2,N,M);
h_sparse = zeros(Sparse.T_im,N,M);
p = zeros(size(h_all,1),N,M);

%%%%%%%%%%�C���p���X�����̃p���[�����ɔ����d�݂̒���%%%%%%%%%
% weight_impulse = zeros(1,L);
% weight_impulse(1:L/2) = -log10(1-exp(-432./(1:L/2)));
% weight_impulse(L/2+1:L) = weight_impulse(L/2)*ones(1,L-L/2);

weight_impulse = zeros(1,8192);
weight_impulse(1:2048) = -log10(1-exp(-216./(1:2048)));
weight_impulse(2049:8192) = weight_impulse(2048)*ones(1,8192-2048);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%�X�V�A���S���Y��%%%%%%%%%%
fprintf('ILRMA+SP Iteration:   ');
for KKK=1:Param.it
    fprintf('\b\b\b\b%4d', KKK);
    for n=1:N
        %Update of basis matrix
        if Param.supervised == false
        T(:,:,n) = max( T(:,:,n) .* sqrt( ( (P(:,:,n) .* (R(:,:,n).^(-2))) * (V(:,:,n)') ) ./ ( (R(:,:,n).^(-1)) * (V(:,:,n)') ) ) , eps);
        R(:,:,n) = T(:,:,n) * V(:,:,n);
        end
        V(:,:,n) = max( V(:,:,n) .* sqrt( ( (T(:,:,n)') * (P(:,:,n) .* (R(:,:,n).^(-2))) ) ./ ( (T(:,:,n)') * (R(:,:,n).^(-1)) ) ) , eps);
       if Param.supervised
           %�w�K�ϊ��Œl������������0��
         if Param.R_x(n) ~= max(Param.R_x)
           V(Param.R_x(n)+1:max(Param.R_x),:,n) = 0;
        end
       end
       
        R(:,:,n) = T(:,:,n) * V(:,:,n);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%% ILRMA3_W(kitahara oda) %%%%%%%%%%%%%%%%%%%%%%%%%%       
    for i=1:I
        for n=1:N
            % U,U_hat�����
            U(:,:,i,n) = (1/J) * ( (X(:,:,i)') * ( X(:,:,i) .* ((R(i,:,n).^(-1))' * ones(1,M)) ) ).';
            %norm(U(:,:,i,n))
            U_hat(:,:,i,n) = U(:,:,i,n) + Sparse.parameter*eye(M);
           
            % U_in,U_in_hat(wni,win_hat)�����(32)(33)
            wni = U_hat(:,:,i,n) \ A(:,n,i); %(32) wni = u
            wni_hat = Sparse.parameter*(U_hat(:,:,i,n) \ (W_hat(n,:,i)'));  %(33) wni_hat = u_hat *W_hat

            %h,h_hat�����(34)(35)
            h    = wni' *U_hat(:,:,i,n) *wni;  %(34)
            h_hat= wni' *U_hat(:,:,i,n) *wni_hat; %(35)

            % W�����(36-1)(36-2)
            if h_hat == 0
                W(n,:,i) = (wni'/(sqrt(h))) + wni_hat';   %(36-1)
            else
                W(n,:,i) =  (conj(h_hat)/(2*h)) *(-1 +sqrt(1+4*h/((abs(h_hat))^2)))* (wni') + wni_hat';  %(36-2)
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        lambda = sqrt( (sum(sum(squeeze(PA(:,n,1))))+ sum(sum(squeeze(PA(:,n,I)))) + 2*sum(sum(squeeze(PA(:,n,2:I-1))))) /L);
        A(:,n,:) = A(:,n,:)/lambda;
        W(n,:,:) = W(n,:,:)*lambda;
        P(:,:,n) = P(:,:,n)*(lambda^2);
        R(:,:,n) = R(:,:,n)*(lambda^2);
        if Param.supervised == false
        T(:,:,n) = T(:,:,n)*(lambda^2);
        end
        V(:,:,n) = V(:,:,n)*(lambda^2);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%% ILRMA3_H(kitahara oda) %%%%%%%%%%%%%%%%%%%%%%%%%%    
a = zeros(2*I-2,N,M);

for m = 1:M
   for n = 1:N
       a(1:L/2+1,n,m) = reshape(A(m,n,:),[L/2+1,1]);
       a(L/2+2:L,n,m) = conj(flip(a(2:L/2,n,m)));
   end
end

for m = 1:M
    for n = 1:N
        h_all(:,n,m) = ifft(a(:,n,m));
        h_all(:,n,m) = weighted_threshhold(h_all(:,n,m),weight_impulse'*Sparse.nu/L);
        H = h_all(:,n,m);
       %H(abs(H)<sqrt(weight_impulse'*Sparse.nu/L)) = 0;
        H(abs(H)<sqrt(weight_impulse'*3/2*Sparse.nu/L)) = 0;
        h_all(:,n,m) = H;
        h_sparse(:,n,m) = h_all(1:Sparse.T_im,n,m);
    end
end

for n = 1:N
    Sum = zeros(size(h_sparse,1),1);
    for m = 1:M
       Sum = Sum + h_sparse(:,n,m).^2;
    end
    lambda_h = sqrt(sum(Sum));
    
    for m = 1:M
        h_sparse(:,n,m) = h_sparse(:,n,m)/lambda_h;
    end
end
   
for m = 1:M
    for n = 1:N
       p(:,n,m) = fft([h_sparse(:,n,m);zeros(L-Sparse.T_im,1)]);
       A_hat(m,n,:) = p(1:I,n,m);
    end
end

    for i=1:I
        W_hat(:,:,i) = pinv(A_hat(:,:,i));
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if drawConv
%     cost(it+1,1) = costFunction_local(P,R,W,I,J);
% end
for i=1:I
    for n=1:N
        Y(i,:,n) = W(n,:,i) * X_tilde2(:,:,i);
    end
end

 fprintf(' ILRMA1_Sparse done.\n');
% if drawConv
%     figure;
%     plot( (0:it), cost );
%     set(gca,'FontName','Times','FontSize',16);
%     xlabel('Iteration','FontName','Arial','FontSize',16);
%     ylabel('Value of cost function','FontName','Arial','FontSize',16);
% end
end
