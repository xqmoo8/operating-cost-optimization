% f = @(x,y,z) x.^2+y.^2+z.^2-10;      % 函数表达式
% [x,y,z] = meshgrid(-10:.2:10,-10:.2:10,-10:.2:10);       % 画图范围
% v = f(x,y,z);
% h = patch(isosurface(x,y,z,v,0)); 
% isonormals(x,y,z,v,h)              
% set(h,'FaceColor','r','EdgeColor','none');
% xlabel('x');ylabel('y');zlabel('z'); 
% alpha(1)   
% grid on; view([1,1,1]); axis equal; camlight; lighting gouraud


cvx_begin
    variable Ppr nonnegative
    variable Pg nonnegative
%     variable P nonnegative
    
%     minimize ( a(1)/a(2) )
    minimize ( (4 * Pg^2 + 3 * Pg + 2) - (Ppr)^(1/3)  )
    subject to
        Pg <= 50
        Pg >= 5
        Ppr <= 5
        Ppr >= 1 
        Pg == Ppr + 5 

cvx_end