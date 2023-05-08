clc;clear all;
twpi = 2*pi;
derad = pi/180;
j = sqrt(-1);
M1 = 8;%x轴阵元数
M2 = 8;%y轴阵元数
kelm = M1+M2-1;%实际阵元数
c = 1500;%预计声速（m/s）
f0 = 1e4;%载波频率（Hz）
fs = 2*f0;%采样频率（Hz）
lambda0 = c/f0;
dd = lambda0/2;%阵元间距
d1 = 0:dd:(M1-1)*dd;
d2 = 0:dd:(M2-1)*dd;
theta = [30 60];
iwave = length(theta);%信源数

err = 0;
for ff=1:10
snr = 20;
snap = 200;%采样数
t = [0:snap-1]/fs;%采样时长
c_real = 1500;%实际声速（m/s）
LAMDA = c_real/f0;%实际波长（m）
A1 = exp(j*twpi*d1'*sin(theta*derad)/LAMDA);%x轴方向向量
A2 = exp(j*twpi*d2'*cos(theta*derad)/LAMDA);%y轴方向向量；
A = [A1;A2(2:end,:)];
S = randn(iwave,snap).*(ones(iwave,1)*exp(j*twpi*(f0*t)));%独立信源

%%%%%接收数据%%%%%
Z = awgn(A*S,snr,'measured');
X = Z(1:M1,:);
Y = Z([1,M1+1:end],:);
Rx = X*X' / snap;
Ry = Y*Y' / snap;
Rxy = X*Y'/ snap;
J = [zeros(M1*M2-1,1),eye(M1*M2-1)];
G = J*kron(Ry.',Rx)*J'/snap;
Q = G^(-1/2);%协方差矩阵的误差的二阶统计特性
kapper = 1e-4;
beta = chi2inv(1-kapper,M1*M2-1);
tic
%%%%%利用CVX工具箱求解凸优化问题%%%%%
cvx_solver sdpt3
cvx_begin sdp quiet
    variable Tu_x(M1,M1) complex hermitian toeplitz semidefinite
    variable Tu_y(M2,M2) complex hermitian toeplitz semidefinite
    variable Rc(M1,M2) complex
    minimize( 1 /(2*sqrt(M1*M2)) * (trace(Tu_x) + trace(Tu_y)) )	%目标函数
    subject to
    [Tu_y,Rc';Rc,Tu_x] == hermitian_semidefinite(M1 + M2);
    square_pos( norm( Q*J*reshape(Rxy - Rc,M1*M2,1), 2) ) <= beta;
cvx_end
toc

Rxy1 = Rc(:,1:M2-1);
Rxy2 = Rc(:,2:M2);
R = [Rxy1;Rxy2];
[U,Sigma,V] = svd(R);
[~,ind] = sort(diag(Sigma),'descend');
U = U(:,ind);
Us = U(:,1:iwave);
U1 = [Us(1:M1-1,:);Us(M1+1:2*M1-1,:)];
U2 = [Us(2:M1,:);Us(M1+2:2*M1,:)];
[T,FAI] = eig(pinv(U1)*U2);
B = Us*T;
doa = zeros(1,iwave);
for ii=1:iwave
   doa(ii) = acot( phase(B(1,ii)*conj(B(1+M1,ii)))/phase(FAI(ii,ii)) )/derad; 
end
doa = sort(doa);
err = err + sum(abs(doa-theta));
end