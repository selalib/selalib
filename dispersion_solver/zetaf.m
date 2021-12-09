function w=zetaf(z)
%--------------------------------------------------
%               2015-04-20 Tsing Wong 
%               2015-05-12 update
% compute the faddeeva function:
%       w(z)=exp(-z^2)erfc(-iz);
% Reference: 
%        [1] M.Zaghloul&A. Ali. "Algorithm 916:
%            computing the Faddeyeva and Voigt func-
%            tions". ACM Trans. Math. Soft. 38 (2), 
%            15 (2011). Preprint available at arXiv:
%            1106.0151
% --------------------------------------------------
x=real(z);y=imag(z);
a=0.5;% a-convergence parameter
x1=x(:);
n=20+max(max(ceil(abs(x1/a))));% n- truncation parameter
%n=20+ceil(abs(x/a));
s1=0;s2=0;s3=0;s4=0;s5=0;
for k=1:n
    %sk=1./(a^2*k^2+y.*y).*exp(-(a^2*k^2+x.*x));
    s1=s1+1./(a^2*k^2+y.*y).*exp(-(a^2*k^2+x.*x));
    s2=s2+1./(a^2*k^2+y.*y).*exp(-(a*k+x).^2);
    s3=s3+1./(a^2*k^2+y.*y).*exp(-(a*k-x).^2);
    s4=s4+a*k./(a^2*k^2+y.*y).*exp(-(a*k+x).^2);
    s5=s5+a*k./(a^2*k^2+y.*y).*exp(-(a*k-x).^2);
end
sin1=sin(x.*y)./(x.*y);
sin1(x.*y==0)=1;% there will be a bug (1 or -1)
sin2=sin(2*x.*y)./(2*x.*y);
sin2(x.*y==0)=1;% there will be a bug (1 or -1)
% caculate the real part of w(z)
Rw1=exp(-x.*x).*erfcx(y).*cos(2*x.*y);
Rw2=2*a*x.*sin(x.*y).*exp(-x.*x).*sin1/pi;
Rw3=(2*a/pi)*(-y.*cos(2*x.*y).*s1+y.*s2/2+y.*s3/2);
Rw=Rw1+Rw2+Rw3;
% caculate the imaginary part of w(z)
Iw1=-exp(-x.*x).*erfcx(y).*sin(2*x.*y);
Iw2=(2*a*x/pi).*exp(-x.*x).*sin2;
Iw3=(2*a/pi)*(y.*sin(2*x.*y).*s1-s4/2+s5/2);
Iw=Iw1+Iw2+Iw3;
% caculate  w(z)
w=Rw+1i*Iw;
% the relationship between the plasma dispersion function 
% and the faddeeva function : zeta=i*sqrt(pi)w(z);
w=1i*sqrt(pi)*w;
end