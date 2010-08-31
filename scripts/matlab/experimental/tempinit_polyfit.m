function beta = tempinit_polyfit(order, uv, X)
%
% TEMPINIT_POLYFIT: fit polynomial to scattered parameterized data
%
%     beta = tempinit_polyfit(order, uv, X)
%
%     computes coefficients of polynomial of given order that gives
%     best fit to data X. 
%

F=zeros(length(uv), order*(order+1)/2);
u = uv(:,1);
v = uv(:,2);

for j=0:order
    for k=0:j
        F(:,1+j*(j+1)/2+k)=(u.^k) .* (v.^(j-k));        
    end
end

beta=F\X;

% ALTERNATIVE CODE

% Fit Fourier harmonics to the data (does that make any sense?)
%nh=12;
%F=zeros(length(u), nh^2);
%ubar = 0.8 * u + 0.1;
%vbar = 0.8 * v + 0.1;
%for j=0:nh-1
%    for k=0:nh-1
%        jk=1+j*nh+k;
%        F(:,jk) = cos(pi * (2 * j + 1) * ubar / 2) .* cos(pi * (2 * k + 1) * vbar / 2);
%    end
%end
%
%beta=F\X; Xfit = F*beta; beta_R=F\R;
%scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
%scatter3(Xfit(:,1),Xfit(:,2),Xfit(:,3),'r.'); hold off; axis image;