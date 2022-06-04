function solutions = FK_Stewart(L1,L2,L3,L4,L5,L6,xsi,ysi,xmi,ymi,rollmin,rollmax,pitchmin,pitchmax,yawmin,yawmax,pxmin,pxmax,pymin,pymax,pzmin,pzmax,legmin,legmax)
    digits(15) 
    
    legmin=2200;                             %leg min/max are defined here, but for a general case they would be defined in the GUI
    legmax=2700;
    
    syms c1 c2 c3                           %set initial symbolic variables

    I=[1 0 0;                               %identity matrix
       0 1 0;
       0 0 1];                              
    C=[0 -c3 c2;                            %used to build rotation matrix(equation 4)
       c3 0 -c1;
      -c2 c1 0];                            
  
    R=(I-C)^-1*(I+C);                       %Rotation matrix(equation 3)
    D=(1+c1^2+c2^2+c3^2);                   %Denominator common to many parts of the program.

    Ta=[1 0 0 -xsi(1);                      %translation matrix from Frame G to Fram A.
        0 1 0 -ysi(1);
        0 0 1 0      ;
        0 0 0 1      ];

    Tb=[1 0 0 -xmi(1);                      %translation matrix from Fram M to Fram B.
        0 1 0 -ymi(1);  
        0 0 1 0      ;
        0 0 0 1      ];

    a1=Ta*[xsi(1);ysi(1);0;1];              
    a2=Ta*[xsi(2);ysi(2);0;1];              
    a3=Ta*[xsi(3);ysi(3);0;1];              
    a4=Ta*[xsi(4);ysi(4);0;1];              
    a5=Ta*[xsi(5);ysi(5);0;1];
    a6=Ta*[xsi(6);ysi(6);0;1]; 

    b1=Tb*[xmi(1);ymi(1);0;1];
    b2=Tb*[xmi(2);ymi(2);0;1];
    b3=Tb*[xmi(3);ymi(3);0;1];
    b4=Tb*[xmi(4);ymi(4);0;1];
    b5=Tb*[xmi(5);ymi(5);0;1];
    b6=Tb*[xmi(6);ymi(6);0;1];

    F2=-a2(1:3).'*R*b2(1:3)+(a2(1:3).'*a2(1:3)+b2(1:3).'*b2(1:3)-L2^2+L1^2)/2;    %F equations as described by equation 14
    F3=-a3(1:3).'*R*b3(1:3)+(a3(1:3).'*a3(1:3)+b3(1:3).'*b3(1:3)-L3^2+L1^2)/2;
    F4=-a4(1:3).'*R*b4(1:3)+(a4(1:3).'*a4(1:3)+b4(1:3).'*b4(1:3)-L4^2+L1^2)/2;
    F5=-a5(1:3).'*R*b5(1:3)+(a5(1:3).'*a5(1:3)+b5(1:3).'*b5(1:3)-L5^2+L1^2)/2;
    F6=-a6(1:3).'*R*b6(1:3)+(a6(1:3).'*a6(1:3)+b6(1:3).'*b6(1:3)-L6^2+L1^2)/2;

    %---------------------------------------------------------------------Phi1
    M=[-a2(1) -a2(2) b2(1) b2(2) F2;        %equation 13
       -a3(1) -a3(2) b3(1) b3(2) F3;
       -a4(1) -a4(2) b4(1) b4(2) F4;   
       -a5(1) -a5(2) b5(1) b5(2) F5;
       -a6(1) -a6(2) b6(1) b6(2) F6];

    Phi1=det(M)*D;                          %Phi1=g1 + g2*c3 + g3*c3^2 + g4*c1^2 + g5*c1*c2 + g6*c2^2=0
    
        %note: subs can do multiple substitutions in the form
        %subs(expression,{old1,old2},{new1,new2})
        %but its performance sometimes erratic.  Hence the nested routines.
        %Also, there are cleaner ways to find g1->g6 and T1->T33,  but
        %substituting in 1's and 0's provides a very fast way to find these
        %values
    
    g1=subs(subs(subs(Phi1,c1,0),c2,0),c3,0);                           %c1/c2/c3=0 to get g1
    g2=subs(subs(subs(subs(Phi1,c2,0),c1,0),c3^2,0)-g1,c3,1);           %c1/c2/c3^2=0, subtract g1 and set c3=1 to get g2
    g3=subs(simplify((subs(subs(Phi1,c1,0),c2,0)-g1)/c3)-g2,c3,1);      %c1/c2=0, subtract g1, divide by c3, subtract g2 and set c3=1 to get g3
    g4=subs(subs(Phi1,{c2,c3},{0,0})-g1,c1,1);                          %c2/c3=0, subtract c1, set c1=1 to get g4
    g5=subs(subs(subs(subs(Phi1,c2^2,0),c1^2,0),c3,0)-g1,{c1,c2},{1,1});%c2^2/c1^2/c3=0, subtract g1 and set c1/c2=1 to get g5
    g6=subs(subs(subs(Phi1,c1,0),c3,0)-g1,c2^2,1);                      %c1/c3=0, subtract g1 and set c2^2=1 to get g6

    %---------------------------------------------------------------------Phi2
    syms ct                                 %ct is a temporary symbolic reference that is used to get the T values.

    Mp=[-a2(1)     -a2(2)       b2(1)       b2(2)     F2;               %equation 18
        -a3(1)     -a3(2)       b3(1)       b3(2)     F3;
        -a4(1)     -a4(2)       b4(1)       b4(2)     F4;
        -a5(1)     -a5(2)       b5(1)       b5(2)     F5;
        (c1-c2*c3) (c2+c1*c3) (-c1-c2*c3) (-c2+c1*c3) 0 ];

    Phi2=det(Mp)*D;                         %Phi2=T1*c1 + T2*c2 + T30*c1^3 + T31*c1^2*c2 + T32*c1*c2^2 + T33*c2^3

    T1=subs(collect(subs(subs(Phi2,c1^3,0),c2,0),c1),c1,1);             %c1^3/c2=0, collect c1's and set c1=1 to get T1
    T2=subs(collect(subs(subs(Phi2,c2^3,0),c1,0),c1),c2,1);             %c2^3/c1=0, collect c1's and set c2=1 to get T2
    T30=subs(subs(subs(subs(subs(Phi2,c1^3,ct),c1,0),c2,0),ct,c1^3),c1,1);  %c1^3=ct, set c1/c2=0, set ct=c1^2 again and set c1=1 to get T30
    T31=subs(subs(subs(subs(subs(subs(subs(subs(subs(subs(expand(subs(subs(Phi2,c3^2,0))-T2*c2),c1^3,0),c3^3,0),c2^2,0),c3^2,0)),c2^3,0),c1^2,ct),c1,0),ct,1),c2,1);    %c3^2=0, subtract T3*c2, set c1^3/c3^3/c2^2/c3^2/c2^3=0, set c1^2 to ct, set c1=0, set ct=1 and c2=1 to get T32...phew    T32=subs(subs(subs(subs(subs(subs(subs(subs(Phi2,c2^2,ct),c2,0),c1^2,0),c1^3,0)-expand(T1*c1),c3^2,0),c3^3,0),ct,1),c1,1);  %c2^2=ct, set c2/c1^2/c1^3=0, subtract T1*c1, set c3^2=0, set ct/c1=1 to get T31
    T32=subs(subs(subs(subs(subs(subs(subs(subs(Phi2,c2^2,ct),c2,0),c1^2,0),c1^3,0)-expand(T1*c1),c3^2,0),c3^3,0),ct,1),c1,1);  %c2^2=ct, set c2/c1^2/c1^3=0, subtract T1*c1, set c3^2/c3^3=0, set ct/c1=1 to get T32
    T33=subs(subs(subs(subs(subs(subs(Phi2,c1,0),c3^3,0),c3^2,0),c2^3,ct),c2,0),ct,1);  %c1/c3^3/c3^2=0, set c2^3=ct, set c2=0 and set ct=1 to get T33

    %-----------------------------------------------------------------Phi3/4/5
    syms Px Py Pz Qx Qy Qz E                %symbolic variables for Px/Py/Pz/Qx/Qy/Qz/E

    P=[Px Py Pz];
    Q=[Qx Qy Qz];

    M3a=[-a2(1) -a2(2) b2(1) b2(2);         %equation 23
         -a3(1) -a3(2) b3(1) b3(2);
         -a4(1) -a4(2) b4(1) b4(2);
         -a5(1) -a5(2) b5(1) b5(2)];
    M3b=[-F2;
         -F3;
         -F4;
         -F5];

    M3=inv(M3a)*M3b;                       %solution to equation 23

    Px=M3(1);
    Py=M3(2);
    Qx=M3(3);
    Qy=M3(4);

    e9=(I-C)*P.'-(I+C)*Q.';                 %equation 9
    e10=e9(1);                              %equation 10
    e11=e9(2);                              %equation 11
    e12=e9(3);                              %equation 12

    qz10=collect(solve(e10,Qz),Pz);         %Solve equation 10 for Pz.
    qz11=collect(solve(e11,Qz),Pz);         %Solve equation 11 for Pz.
    qz12=collect(solve(e12,Qz),Pz);         %Solve equation 12 for Pz.

    pz1=factor(eval(solve(qz11-qz12,Pz)));  %Equate eq11 and eq12 through pz, so eq11-eq12=0, equation 25 
    pz2=factor(eval(solve(qz10-qz12,Pz)));  %Equate eq10 and eq12 through pz, so eq10-eq12=0, equation 26

    N1=pz1*c1*D;                            %Numerator of pz1
    D1=c1*D;                                %denominator of pz1
    N2=pz2*c2*D;                            %numerator of pz2
    D2=c2*D;                                %denominator of pz2

    pz3=simplify((N1-N2)/(D1-D2));          %pz3 is generated with N1,N2,D1,D2.

    M4a=[g4  g5  g6  0   0 ;                %equation 31 M4a*M4=M4b
         0   g4  g5  g6  0 ;
         0   0   g4  g5  g6;
         T30 T31 T32 T33 0;
         0   T30 T31 T32 T33];

    M4b=[c1^2*g1+c1^2*g2*c3+c1^2*g3*c3^2   ;%equation 31
         c1*c2*g1+c1*c2*g2*c3+c1*c2*g3*c3^2;
         c2^2*g1+c2^2*g2*c3+c2^2*g3*c3^2   ;
         T1*c1^2+T2*c1*c2                  ;
         T1*c1*c2+T2*c2^2                  ];
    M4=vpa(inv(M4a)*(-M4b));                %solution to equation 31

    [M4num(1:5),M4den(1:5)]=numden(simplify(M4(1:5)));  %break up the solution M4 into numerator and denominator
    M4e(1:5)=M4num(1:5)/E;                              %replace the denominator with E, this is to speed calculation.

    Phi=[simplify((Px^2+Py^2+pz1^2-L1^2)*D^2*c1^2)     ;     %Phi3/4/5=PX^2 + PY^2 + PZ^2 - L1^2=0
         simplify((Px^2+Py^2+pz2^2-L1^2)*D^2*c2^2)     ;     %other values were added to remove denominators
         simplify((Px^2+Py^2+pz3^2-L1^2)*D^2*(c1-c2)^2)];    
        
        %format for the algsubs function is algsubs(old=new,expression).
        %Since algsubs is a Maple command, it only accepts string input, so
        %to get the old=new string expression, I concat the old with the
        %new lower order expression via ['old', '=' char(new)].  So the
        %format of the maple command is now:
        %maple('algsubs',['old', '=' char(new)],expression,'exact'). The
        %'exact' on the end forces maple to not try and swap 'old' and
        %'new' if for some reason it cannot replace the old expression.
    for k=1:2                               %this routine substitutes the high order expressions in Phi3/4/5 for the 
        for j=1:3                           %low order expressions found in M4
            Phi(j)=maple('algsubs',['c1^4',     '=' char(M4e(1))],simplify(Phi(j)),'exact');
            Phi(j)=maple('algsubs',['c1^3*c2',  '=' char(M4e(2))],simplify(Phi(j)),'exact');
            Phi(j)=maple('algsubs',['c1^2*c2^2','=' char(M4e(3))],simplify(Phi(j)),'exact');
            Phi(j)=maple('algsubs',['c1*c2^3',  '=' char(M4e(4))],simplify(Phi(j)),'exact');
            Phi(j)=maple('algsubs',['c2^4',     '=' char(M4e(5))],simplify(Phi(j)),'exact');
        end
    end

    Phi3=simplify(subs(Phi(1),E,M4den(2))); %replace E with the actual denominator found earlier.
    Phi4=simplify(subs(Phi(2),E,M4den(2))); %I used M4den(2) because at low digit accuracy the first
    Phi5=simplify(subs(Phi(3),E,M4den(2))); %value M4den(1) will drift slightly

    E=M4den(2);                             %replace E with denominator of M4
    [W,Wd]=numden([Phi3 Phi4 Phi5]);        %since Phi3/4/5 are equal to 0, we can discard the denominators, assuming they never equal 0

    W11=subs(subs(W(1),c1,0),c2,0);         %W's are used in the final polynomial, W*3 is found last because it uses the previous W's
    W12=subs((subs(W(1),c2,0)-W11),c1,1);
    W14=subs((subs(W(1),c1,0)-W11),c2,1);
    W13=subs(subs(W(1),c1,1),c2,1)-W11-W12-W14;

    W21=subs(subs(W(2),c1,0),c2,0);
    W22=subs((subs(W(2),c2,0)-W21),c1,1);
    W24=subs((subs(W(2),c1,0)-W21),c2,1);
    W23=subs(subs(W(2),c1,1),c2,1)-W21-W22-W24;

    W31=subs(subs(W(3),c1,0),c2,0);
    W32=subs((subs(W(3),c2,0)-W31),c1,1);
    W34=subs((subs(W(3),c1,0)-W31),c2,1);
    W33=subs(subs(W(3),c1,1),c2,1)-W31-W32-W34;

    S=[(g1+g2*c3+g3*c3^2)   g4        g5          g6     ;%equation 36
       W11/E^3              W12/E^3   W13/E^3     W14/E^3;
       W21/E^3              W22/E^3   W23/E^3     W24/E^3;
       W31/E^3              W32/E^3   W33/E^3     W34/E^3];
    Sd=factor(det(S));                      %this is the final polynomial
    c3sol=solve(Sd);                        %solve for c3

    for i=size(c3sol,1):-1:1                %for each value in c3sol
        if isreal(c3sol(i))==0              %if the value isn't real
            c3sol(i)=[];                    %delete that row
        end
    end


    Sp=[S(1,2) S(1,3) S(1,4);               %Rearange matrix S to solve for c1^2,c1*c2, and c2^2
        S(2,2) S(2,3) S(2,4);
        S(3,2) S(3,3) S(3,4)];
    Spx=[-S(1,1);-S(2,1);-S(3,1)];
    Sol=inv(Sp)*Spx;                        %Sol(1)=c1^2, Sol(2)=c1*c2 and Sol(3)=c2^2

    Sm=[];                                  %create an empty matrix

    C_tol=10^-30;  %tolerance for cases when c1,c2,c3 are close to 0

    for i=1:size(c3sol,1)                                           %for each value of c3
        c1sol=sqrt(subs(Sol(1),c3,c3sol(i)));                       %substitute c3(i) into Sol(1) to get c1
        c1c2=subs(Sol(2),c3,c3sol(i));                              %substitute c3(i) into Sol(2) to get c2*
        c2sol=sqrt(subs(Sol(3),c3,c3sol(i)));                       %substitute c3(i) into Sol(3) to get c2
        if (isreal(c1sol)~=0) && (isreal(c2sol)~=0)                 %if c1 and c2 are real
            if double(c1sol) > C_tol                                %and if c1 is bigger than the tolerance value
                c2sol=c1c2/c1sol;                                   %c2=c1*c2/c1
                pxsol=subs(Px,{c1,c2,c3},{c1sol,c2sol,c3sol(i)});   %sub c1/c2/c3 into Px to get px value
                pysol=subs(Py,{c1,c2,c3},{c1sol,c2sol,c3sol(i)});   %sub c1/c2/c3 into Py to get py value
                pzsol=subs(pz1,{c1,c2,c3},{c1sol,c2sol,c3sol(i)});  %sub c1/c2/c3 into Pz to get pz value
                Sm=[Sm;c1sol c2sol c3sol(i) pxsol pysol pzsol];     %add values to solution matrix
                Sm=[Sm;-c1sol -c2sol c3sol(i) pxsol pysol -pzsol];  %add negative case values to matrix
            else
                if double(c2sol) > C_tol                            %if c2 is bigger than the tolerance
                    pxsol=subs(Px,{c1,c2,c3},{0,c2sol,c3sol(i)});   %sub c1/c2/c3 into Px to get px value
                    pysol=subs(Py,{c1,c2,c3},{0,c2sol,c3sol(i)});   %sub c1/c2/c3 into Py to get py value
                    pzsol=subs(pz2,{c1,c2,c3},{0,c2sol,c3sol(i)});  %sub c1/c2/c3 into pz2 to get pz value
                    Sm=[Sm;0 c2sol c3sol(i) pxsol pysol pzsol];     %add values to solution matrix
                    Sm=[Sm;0 -c2sol c3sol(i) pxsol pysol -pzsol];   %add negative case values to matrix
                else
                    if c3sol(i) < C_tol                             %if c3 is below tolerance
                        c3sol(i)=0;                                 %make c3=0
                    end
                    pxsol=subs(Px,{c1,c2,c3},{c1sol,0,c3sol(i)});   %sub c1/c2/c3 into Px to get px value
                    pysol=subs(Py,{c1,c2,c3},{c1sol,0,c3sol(i)});   %sub c1/c2/c3 into Py to get py value
                    pzsol=sqrt(L1^2-pxsol^2-pysol^2);               %using px and py, get pz
                    Sm=[Sm;c1sol 0 c3sol(i) pxsol pysol pzsol];     %add values to solution matrix
                end
            end
        end
    end

    for i=1:size(Sm,1)                                          %foreach value of c1
        Rsolve=subs(subs(subs(R,c1,Sm(i,1)),c2,Sm(i,2)),c3,Sm(i,3)); %substitute c1/c2/c3 into Rotation matrix
        Rsolve=double(Rsolve);                                  %convert symbolic matrix to numerical
        Sm(i,1)=vpa(-atan2(Rsolve(2,3),Rsolve(3,3))*180/pi);    %replace c1/c2/c3 with X degree in the solution matrix
        Sm(i,2)=vpa(asin(Rsolve(1,3))*180/pi);                  %replace c1/c2/c3 with Y degree in the solution matrix
        Sm(i,3)=vpa(-atan2(Rsolve(1,2),Rsolve(1,1))*180/pi);    %replace c1/c2/c3 with Z degree in the solution matrix
    end

    for i=1:size(Sm,1)              %for each solution
        trans_coords=coordtrans(Sm(i,1),Sm(i,2),Sm(i,3),Sm(i,4),Sm(i,5),Sm(i,6),xsi,ysi,xmi,ymi); %transform the coordinates from those used in the forward kinematics analysis based on the Lee抯 paper to those of the CDSL system
        Sm(i,4)=trans_coords(1);    %add the x component of the position vector rm
        Sm(i,5)=trans_coords(2);    %add the y component of the position vector rm
        Sm(i,6)=trans_coords(3);    %add the z component of the position vector rm
    end

    for i=size(Sm,1):-1:1           %foreach solution
        if rollmin  > double(Sm(i,1)) || double(Sm(i,1)) > rollmax  ||...   %if values are not within min/max
           pitchmin > double(Sm(i,2)) || double(Sm(i,2)) > pitchmax ||...
           yawmin   > double(Sm(i,3)) || double(Sm(i,3)) > yawmax   ||...
           pxmin    > double(Sm(i,4)) || double(Sm(i,4)) > pxmax    ||...
           pymin    > double(Sm(i,5)) || double(Sm(i,5)) > pymax    ||...
           pzmin    > double(Sm(i,6)) || double(Sm(i,6)) > pzmax
                Sm(i,:)=[];         %delete that solution
        end
    end

    for i=size(Sm,1):-1:1           %for each solution
        leg_length=stew_inverse(xsi,ysi,xmi,ymi,Sm(i,1),Sm(i,2),Sm(i,3),Sm(i,4),Sm(i,5),Sm(i,6));%get leg length of solutions from inverse kinematics
        if double(leg_length(1)) < legmin || double(leg_length(1)) > legmax || abs(double(leg_length(1))-L1) > 0.5 ||...    %if the legs are too short, too long, or not the same as the input
           double(leg_length(2)) < legmin || double(leg_length(2)) > legmax || abs(double(leg_length(2))-L2) > 0.5 ||...
           double(leg_length(3)) < legmin || double(leg_length(3)) > legmax || abs(double(leg_length(3))-L3) > 0.5 ||...
           double(leg_length(4)) < legmin || double(leg_length(4)) > legmax || abs(double(leg_length(4))-L4) > 0.5 ||...
           double(leg_length(5)) < legmin || double(leg_length(5)) > legmax || abs(double(leg_length(5))-L5) > 0.5 ||...
           double(leg_length(6)) < legmin || double(leg_length(6)) > legmax || abs(double(leg_length(6))-L6) > 0.5
                Sm(i,:)=[]; %delete that solution.
        end
    end


    solutions=Sm   %return the value

