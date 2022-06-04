function Pset = cirintsuv(R1,R2,u0)
if((u0>R1+R2)||(u0<abs(R2-R1)))
    Pset=nan;
    
elseif((u0==R1+R2)||(u0==abs(R2-R1)))
    Pset=[(R1^2-R2^2+u0^2)/2/u0;0];
    
else
    u=(R1^2-R2^2+u0^2)/2/u0;
    Pset=[u,u;sqrt(R1^2-u^2),-sqrt(R1^2-u^2)];
    
end 