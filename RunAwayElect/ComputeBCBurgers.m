function bc = ComputeBCBurgers(Lev,Deg,LInt,LEnd,Fun,time,bcL,bcR,MaxC)
%function ComputeBC to compute the bc term

%-----------------------------------------------------


L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

quad_num = 10;

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
% [quad_x,quad_w] = lgwt(quad_num,-1,1);
% p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);
%
% Jacobi = h/2;

p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

WorkCel = 0;

p_val  = legendre(-1,Deg)  * 1/sqrt(h);

% alpha = FunCoef(LInt);
c = [1:Deg];
if bcL == 0 %p_val*uold(c) >= 0 %&& bcL == 0 %|| alpha-FluxVal*abs(alpha) == 1% left cell is Dirichlet boundary
%     AvgF = (Fun(LInt,time)).^2/2;
%     JumU = (Fun(LInt,time) );
%     TraceFL = AvgF + MaxC/2 * JumU;
%     IntVal = p_L'*TraceFL;
    
    IntVal = p_L'*(Fun(LInt,time)).^2/2;
    %     IntVal =  p_val'*(Fun(LInt,time)) ;
    bc(c) = - IntVal;
    
end

WorkCel = Tol_Cel_Num - 1;

c = Deg*WorkCel+[1:Deg];
if bcR == 0 %p_val*uold(c) < 0 && bcR == 0 %|| alpha + FluxVal*abs(alpha) == 1% right cell is Dirichlet boundary
%     AvgF = (  Fun(LEnd,time).^2)/2;
%     JumU = (- Fun(LEnd,time) );
%     TraceFR = AvgF + MaxC/2 * JumU;
%     IntVal = p_R'*TraceFR;
    
    IntVal = p_R'*(Fun(LEnd,time)).^2/2;
    %     IntVal =  p_val'*(Fun(LEnd,time));
    bc(c) = IntVal;
end


end