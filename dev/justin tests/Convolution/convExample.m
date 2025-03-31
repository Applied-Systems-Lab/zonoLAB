clear, clc, close all;
%%
load("convtest.mat")
%%
A = testPoint;
B = randn(3,7);

Cfull = conv2(A,B,"full");
Csame = conv2(A,B,"same");
Cvalid = conv2(A,B,"valid");

%%
convLength = size(B,2);
convHeight = size(B,1);
Empty = zeros(size(A,1)+2*(convHeight-1),size(A,2)+2*(convLength-1));
k = 1;
for i = 1: size(Empty,1)-size(B,1)+1
    for j = 1: size(Empty,2)-size(B,2)+1
        temp = Empty;
        temp(i:i+(convHeight-1),j:j+(convLength-1))= rot90(B,2);

        convShape(k,:) = reshape(temp',[size(temp,1)*size(temp,2) 1]);
        k = k + 1;
    end
end

Ashape = Empty;
Ashape(size(B,1):size(Ashape,1)-size(B,1)+1,size(B,2):size(Ashape,2)-size(B,2)+1) = A;
Ashape = reshape(Ashape',[size(Ashape,1)*size(Ashape,2) 1]);

convOut = convShape*Ashape;
convSquare = reshape(convShape*Ashape,[size(A,2)+size(B,2)-1 size(A,1)+size(B,1)-1])';

startH = floor(convHeight/2) + 1;
startL = floor(convLength/2) + 1;

TESTER = convSquare(startH:startH+size(testPoint,1)-1,startL:startL+size(testPoint,2)-1);
%%Comparison
norm(Cfull-convSquare,2)

% conv(A,B) = conv(B,A)

% One matrix becomes toeplitz Any convolution toeplitz matrix the size of
% If A is an n x m matrix and B is a p x q matrix
% a single matrix can be constructed to calculate convolution with the dimensions
% n+p-1 by p

% zfill = zeros(size(A,2)+size(B,1)-1,size(B,1));
% 
% z = zeros(1,size(B,1)-1);
% a1toe = toeplitz([a1 z],[a1(1,1) z]);
% a2toe = toeplitz([a2 z],[a2(1,1) z]);
% a3toe = toeplitz([a3 z],[a3(1,1) z]);
% 
% M = [   a1toe zfill zfill zfill;
%         a2toe a1toe zfill zfill;
%         a3toe a2toe a1toe zfill;
%         zfill a3toe a2toe a1toe;
%         zfill zfill a3toe a2toe;
%         zfill zfill zfill a3toe;
%         ];
% 
% % Matrix gets Stacked
% Brow = [b1; b2; b3; b4];
% 
% matrixConvRow = M*Brow;
% 
% matrixConv = reshape(matrixConvRow, 6,6);
% 
% matrixConv' - Cfull
% 
% 
% %%
% clear; clc;
% A = [1,2;3,4;5,6;7,8];
% B = [1,2,3,4,5;6,7,8,9,10;11,12,13,14,15];
% 
% convVal = conv2(A,B);
% 
% [n,m] = size(A);
% [p,q] = size(B);
% 
% for i = 1:n
%     Arow{i} = A(i,:);
% end
% 
% for i = 1:p
%     Brow{i} = B(i,:);
% end
% % size of initial toplitez
% % m+q-1 by q
% topli = zeros(m+q-1,q);
% 
% 
% 
% 
% % size of single matrix toplitez