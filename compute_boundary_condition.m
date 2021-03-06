function bc = compute_boundary_condition(pde,g_func,t,Lev,Deg,xMin,xMax,Fun,LorR)
% function bc = ComputeBC(Lev,Deg,xMin,xMax,Fun,time,LorR)
% function ComputeBC to compute the bc term
% This is the evaluation for two points on 1D
% Func*v|_xMin and Func*v|_xMax
%-----------------------------------------------------

p = pde.params;

L = xMax-xMin;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
small_dx = h*1e-7; %Parameter to fix problem with Inhomogeneous BC
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

p_L = lin_legendre(-1,Deg) * 1/sqrt(h); % TODO : this happens in multiple places. Consolidate. 
p_R = lin_legendre(+1,Deg) * 1/sqrt(h);

if strcmp(LorR,'L')
    
    WorkCel = 0;
    c = [1:Deg];
    IntVal =  p_L'*Fun(xMin,p,t) ;
    if isfinite(g_func(xMin,p,t)) %make sure g_func is finite
        
        bc(c) = -g_func(xMin,p,t).*IntVal;
    
    else
        
        xLclose = xMin + small_dx; %Move away from xMin by small_dx
        bc(c) = -g_func(xLclose,p,t).*IntVal;
        
    end
else
    
    WorkCel = Tol_Cel_Num - 1;
    c = Deg*WorkCel+[1:Deg];
    IntVal =  p_R'*Fun(xMax,p,t);
    
    if isfinite(g_func(xMax,p,t)) %make sure g_func is finite
        
        bc(c) = g_func(xMax,p,t).*IntVal;
        
    else
        
        xRclose = xMax - small_dx; %Move away from xMax by small_dx
        bc(c) = g_func(xRclose,p,t).*Intval;
        
    end
    
end

end
