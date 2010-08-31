function X = tempinit_polyinterp(order, beta, uv)
%
% TEMPINIT_POLYINTERP: interpolate using polynomial fit
%
%     X = tempinit_polyinterp(order, beta, uv)
%
%     Computes X(uv) based on polynomial coefficients beta
%
u = uv(:,1);
v = uv(:,2);

% Interpolate the x,y,z of the nodes
Fn=zeros(length(u), order*(order+1)/2);
for j=0:order
    for k=0:j
        Fn(:,1+j*(j+1)/2+k)=(u.^k) .* (v.^(j-k));
    end
end

X = Fn * beta;