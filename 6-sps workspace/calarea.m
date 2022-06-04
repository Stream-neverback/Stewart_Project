function [ out ] = calarea( in,flg )
%calculate area for arcsets
out = 1/2 * sum( ( in(2,:).*in(3,:).*(-cos(in(5,:))+cos(in(4,:)))...
    +in(1,:).*in(3,:).*(sin(in(5,:))-sin(in(4,:)))...
    +in(3,:).^2.*(in(5,:)-in(4,:)) ).*flg );

end