% kappa-L method for boundary evolution
function [body_p1,body,f_m1] = kappaL(body,body_m1,f_m1,Vcm,dt,s,iks,kinv_h,filter)
    
    global beta eps ds
    
    % load in values
    X = body.X; %parametrization
    Xcm = centroid(body,iks); %center of mass
    Ns = length(X); %number of grid points
    kss = iks.*iks;

    % curvature
    kappa = body.kappa;
    kappa_m1 = body_m1.kappa;
    kappa_h = fft(kappa);    
    kappa_hm1 = fft(kappa_m1);

    L = body.L; %arclength    
    theta0 = body.theta0; %tangent angle at s = 0
    
    ohm = body.ohm; ohm_m1 = body_m1.ohm; %angular velocity
    fkappa_hm1 = f_m1.kappa; %rhs of kappa eqn
    fL_m1 = f_m1.L; %rhs of L eqn

    %% Update curvature and arclength
    
    % normal velocity
    L0 = body.L0;
    Vn = beta*(L0/L)^(1/4)*ones(Ns,1);
    Vn = Vn.*sigmoid(X(:,2),1,0);  
    
    % Compute tangential velocity
    Iarg = sum(kappa.*Vn)*ds;
    P_h = 1./kss;
    P_h(1,1) = 0;
    Vs = real(ifft(P_h.*iks.*fft(kappa.*Vn)));    
        
    Vne = [Vn(end);Vn;Vn(1)];
    Vn_ss = diff(Vne,2)/ds^2;

    kVs = kappa.*Vs;
    kVse = [kVs(end);kVs;kVs(1)];
    kVs_s = (kVse(3:end) - kVse(1:end-2))/(2*ds);
    
    % Update curvature
    fkappa_h = (1/L)*fft(Vn_ss + kVs_s);
    rhs = 4*kappa_h - kappa_hm1 + 2*dt*(2*fkappa_h - fkappa_hm1).*filter;
    kappa_h = rhs./(3 - 2*dt*kss*eps/L^2);
    kappa = real(ifft(kappa_h));
    
    % Update arclength
    fL = -Iarg;
    L = ab2(L,fL,fL_m1,dt);
    
    %% Update anchor and tangent angle
    theta0 = ab2(theta0,ohm,ohm_m1,dt);
    theta = real(ifft(kinv_h.*kappa_h)) + theta0 + 2*pi*s;
    
    % Get anchor coordinate and advance
    x0 = X(1,1);
    y0 = X(1,2);
    
    fx0 = Vn(1)*cos(theta0) + Vcm(1);
    fy0 = Vn(1)*sin(theta0) + Vcm(2);
    fx0_m1 = f_m1.x0; 
    fy0_m1 = f_m1.y0;
    
    x0_p1 = ab2(x0,fx0,fx0_m1,dt);
    y0_p1 = ab2(y0,fy0,fy0_m1,dt);
            
    % apply rotation to anchor point
    rx = x0_p1 - Xcm(1);
    ry = y0_p1 - Xcm(2);    
    
    x0_p1 = cos(dt*ohm)*rx - sin(dt*ohm)*ry + Xcm(1);
    y0_p1 = sin(dt*ohm)*rx + cos(dt*ohm)*ry + Xcm(2);
        
    % update body
    rhsx_h = fft(cos(theta));
    rhsy_h = fft(sin(theta));
    x = real(ifft(kinv_h.*rhsx_h));
    y = real(ifft(kinv_h.*rhsy_h));
    X = L*[x y];
    X = X - X(1,:) + [x0_p1,y0_p1];
    
    body_p1 = struct('X',X,'kappa',kappa,'L',L,'theta',theta,'theta0',theta0,'ohm',ohm,'L0',L0);
    f_m1 = struct('x0',fy0,'y0',fx0,'kappa',fkappa_h,'L',fL);
    
end
