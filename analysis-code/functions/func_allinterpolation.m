%--------------------------------------------------------------------------
% Conduct interpolation based on the chosen interpolation_method
% Functions are used from func_itp.m
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------

function [pt,res] = func_allinterpolation(intp,natpara,ex,ey,s,xx)
    % interpolation_method: Fourier, Linear, or Polynomial
    % intp: function register
    % numpt: number of interpolation point
    % ex: x value
    % ey: y value 
    % s: arclength in summation
    % xx: equidistributed arclength
    % variationset: range of order

    global interpolation_method

    error_var_set = [];
    res = cell(1,5);

    % Linear
    if interpolation_method(1) == 'L'
        pt = intp.interparc(natpara.numpt,ex,ey,'linear');
        pt(:,2) = -pt(:,2);

    % Fourier
    elseif interpolation_method(1) == 'F' 
        for i=1:length(natpara.vrr)
            Nk = natpara.vrr(i); % Fourier mode; 
            [~,~,~,err1,~,~] = intp.interpft(ex,s,Nk,natpara.numpt);
            [~,~,~,err2,~,~] = intp.interpft(ey,s,Nk,natpara.numpt);
            allerr = sqrt(err1.^2+err2.^2);
            error_var_set(end+1) = allerr;
        end
        
        minerror = min(error_var_set);
        % restrict the fourier node of being too large and prevent
        % of overfitting
        upper = minerror; 
        ktt = (error_var_set <= upper);
        for i=1:length(ktt)
            if(ktt(i) == 1)
                firstind = i; % optimal solution 
                break
            end
        end

        final_deg = natpara.vrr(firstind);
        [x_e,xs_e,xss_e,~,x_h,k1] = intp.interpft(ex,s,final_deg,natpara.numpt);
        [y_e,ys_e,yss_e,~,y_h,~] = intp.interpft(ey,s,final_deg,natpara.numpt);
        % interpolation point
        pt = [x_e,-y_e];
        fd = [xs_e,-ys_e]; % first derivative
        sd = [xss_e,yss_e]; % second derivative

        % calculate curve using second derivatives
        curv = sqrt(sd(:,1).^2 + sd(:,2).^2);

        % Fourier nodes
        k1 = k1';
        xss = real(ifft(-k1.^2.*x_h));
        yss = real(ifft(-k1.^2.*y_h));
        kappa = sqrt(xss.^2 + yss.^2);
        kappa_h = fft(kappa);

        res{1} = kappa; res{2} = kappa_h; res{3} = curv; res{4} = x_e;
        res{5} = y_e; res{6} = fd; res{7} = sd;

    % Polynomial
    elseif interpolation_method(1) == 'P' 
        % iterate through the whole set of degree
        for i=1:length(natpara.vrr)
            degree = natpara.vrr(i);
            [px,Sx] = polyfit(s,ex,degree);
            [py,Sy] = polyfit(s,ey,degree);
            % update the range of parameterization [0,finalL]
            [~,errx] = polyval(px,xx,Sx);
            [~,erry] = polyval(py,xx,Sy);
            errort = [errx;erry];
            
            sqerror = sqrt(errort(1,:).^2+errort(2,:).^2);
            totalerror = sqrt(sum(sqerror.^2));

            error_var_set(end+1) = totalerror;
        end

        minerror = min(error_var_set);
        upper = minerror;
        ktt = (error_var_set <= upper);
        for i=1:length(ktt)
            if(ktt(i) == 1)
                firstind = i; % optimal solution 
                break
            end
        end

        final_deg = natpara.vrr(firstind);

        [px,Sx] = polyfit(s,ex,final_deg);
        [py,Sy] = polyfit(s,ey,final_deg);
        % update the range of parameterization [0,finalL]
        [polyx,~] = polyval(px,xx,Sx);
        [polyy,~] = polyval(py,xx,Sy);
        pt = [polyx;polyy]'; % coordinates of interpolation

    end

end