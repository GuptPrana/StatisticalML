function bp = DwShell(A, varargin)
%DWSHELL returns the boundary points of the input matrix. This function is
%based on the combination of rotations around z axis and imaginary axis (y
%axis). The user can specify two angular resolution, since the algorithm is
%quite brutal O(n^2), no more than 200 samples in each rotating direction
%is suggested. The rotations are uniformly distributed from 0 to 2*pi.
% dwshell(A) returns 50x50 (xy,xz) boundary points of DW(A)
% dwshell(A,100) returns 100x100 (xy,xz) boundary points of DW(A)
% dwshell(A,100,150) returns 100x150 (xy,xz) boundary points of DW(A)

%% input processing
switch nargin
    case 1
        n1 = 50; n2 = n1;
    case 2
        n1 = varargin{1}; n2 = n1;
    case 3
        n1 = varargin{1}; n2 = varargin{2};
    otherwise
        error("Wrong number of inputs.")
end

%% rotation
n = size(A,1); % the size of the square matrix
phi = linspace(0,2*pi,n1); % rotation angles around z axis
theta = linspace(0,2*pi,n2); % rotationi angles around y axis
bp = zeros(3,n1*n2);
dwTup = [(A+A')/2;(A-A')/(2*1i);A'*A]; % tuple matrices of DW shell
for i = 1:n2
    TXZ = [cos(theta(i)),0,-sin(theta(i));...
            0,1,0;...
            sin(theta(i)),0,cos(theta(i))];
    mRotatedXZ = kron(TXZ, eye(n))*dwTup;
    for j = 1:n1
        TXY = [cos(phi(j)), -sin(phi(j)), 0;...
            sin(phi(j)), cos(phi(j)), 0;...
            0, 0, 1];
        mRotated = kron(TXY,eye(n))*mRotatedXZ;
        
        HRotated = mRotated(1:n,:); % rotated first Hermitian block
        [hEigVec,hEig] = eig(HRotated);
        [d,ind] = sort(diag(hEig),'descend');
        hEig = diag(d);
        hEigVec = hEigVec(:,ind); % sort the eigenvalues and the associated eigenvectors
        xB = hEigVec(:,1); % choose the eigenvector corresponding to the maximum eigenvalue
        bp(:,n1*(i-1)+j) = real(kron(eye(3),xB')*dwTup*xB);
%         pRotated = kron(eye(3),xB')*mRotated*xB; % three tuples
%         % rotate the boundary point back
%         bp(:,n1*(i-1)+j) = real(TXZ\(TXY\pRotated)); % bp = boundary points
    end
end

