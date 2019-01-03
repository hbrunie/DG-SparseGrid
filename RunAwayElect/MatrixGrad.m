function Mat = MatrixGrad(Lev,Deg,LInt,LEnd,FunCoef)
L = LEnd-LInt;
Tol_Cel_Num=2^(Lev);
h=L/tol_cel_num;
DoF = Deg * Tol_Cel_Num;

Mat = sparse(DoF,DoF);

quad_num = 10;
FluxVal = 1;

% compute the trace values
p_L = legendre(-1,deg) * 1/sqrt(h);
p_R = legendre(+1,deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val  = legendre(quad_x,deg)  * 1/sqrt(h);
Dp_val = dlegendre(quad_x,deg) * 1/sqrt(h) * 2/h;


for WorkCel=0:Tol_Cel_Num-1
    
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    xL = LInt+WorkCel*h;
    xR = x0+h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    xMid= (xR+xL)/2;
    
    IntVal=[Dp_val'*(quad_w.*FunCoef(PhyQuad).*p_val)];
    
    c = Deg*WorkCel+[1:Deg];
    
    Mat = Mat + sparse(c'*ones(1,Deg),ones(Deg,1)*c,IntVal,DoF,DoF);
 
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    TraVal = [-p_L' * FunCoef(xL) * ( p_R/2 - FluxVal/2*p_R),...
              -p_L' * FunCoef(xL) * ( p_L/2 + FluxVal/2*p_L),... % xL
               p_R' * FunCoef(xR) * ( p_R/2 - FluxVal/2*p_R),...
               p_R' * FunCoef(xR) * ( p_L/2 + FluxVal/2*p_L),... % xR
            ];

    % Adding trace value to matrix
    RowInd = [c' c' c' c']*ones(4,Deg);
    ColInd = ones(Deg,4)*[c-Deg,c,c,c+Deg];
    Mat = Mat + sparse(RowInd,ColInd,TraVal,DoF,DoF);
    
    
end

end