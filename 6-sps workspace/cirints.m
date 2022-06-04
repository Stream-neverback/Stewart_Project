function Pset = cirints(C1,R1,C2,R2)
if((C1(1)==C2(1))&&(C1(2)==C2(2)))
    Pset=nan;
    
else
    U2=transuv(C1,C2,C2);
    Pset=cirintsuv(R1,R2,U2(1));
end
if(isnan(Pset))
else
    Pset=transxy(C1,C2,Pset);
end