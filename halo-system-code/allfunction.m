%--------------------------------------------------------------------------
% Halo system for the ice flipping experiments
% Function registeration

% Crafted by Scott Weady
% Revised by Zihan Zhang, Courant Institute
% Updated July 2023
%--------------------------------------------------------------------------

function func_result = allfunction
    func_result.updateFig = @updateFig;
    func_result.changeY = @changeY;
    func_result.changeX = @changeX;
    func_result.changeRadius = @changeRadius;
    func_result.changeDecay = @changeDecay;
    func_result.calcI = @calcI;
end


function updateFig(M)
    global im cmap
    figure(2);
    hold on
    cla
    im = imshow(M);
end

function changeY(source,~)
    global params
    Y = get(source,'value');
    X = params.X;
    params.Y = Y;
    R = params.R;
    A = params.A;
    I = calcI(X,Y,R,A);
    updateFig(I)
end

function changeX(source,~)
    global params
    X = get(source,'value');
    params.X = X;
    Y = params.Y;
    R = params.R;
    A = params.A;
    I = calcI(X,Y,R,A);
    updateFig(I)
end

function changeRadius(source,~)
    global params
    R = get(source,'value');
    X = params.X;
    Y = params.Y;
    params.R = R;
    A = params.A;
    I = calcI(X,Y,R,A);
    updateFig(I)
end

function changeDecay(source,~)
    global params
    A = get(source,'value');
    X = params.X;
    Y = params.Y;
    R = params.R;
    params.A = A;
    I = calcI(X,Y,R,A);
    updateFig(I)
end

function I = calcI(X,Y,R,a)
    global Nx Ny
    x = linspace(-1,1,Nx);
    y = linspace(-(Ny/Nx),(Ny/Nx),Ny);
    [xx,yy] = meshgrid(x,y);
    xx = xx'; yy = yy';
    rr = sqrt((xx + X).^2 + (yy + Y).^2);
    I = 1-1./(1 + exp(-50*a*(rr-R)));
    updateFig(I);  
end



