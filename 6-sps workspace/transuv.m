function p2 = transuv(c1,c2,p1)
m=size(p1,2);
P=[p1;ones(1,m)];
if(~((c1(1)==c2(1))&&(c1(2)==c2(2))))
    t=atan2((c2(2)-c1(2)),(c2(1)-c1(1)));
    A=[cos(t),sin(t),-c1(1)*cos(t)-c1(2)*sin(t);
       -sin(t),cos(t),c1(1)*sin(t)-c1(2)*cos(t);
        0,      0,              1               ];
else
    A=[1,0,-c1(1);
       0,1,-c1(2);
       0,0,   1  ];
end
p=A*P;
p2=[p(1,:);p(2,:)];