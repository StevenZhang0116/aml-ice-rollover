%--------------------------------------------------------------------------
% Calculate harmonic moment of the boundary curve [outdated]
%
% Steven Zhang, Courant Institute
% Updated July 2022
%--------------------------------------------------------------------------

function [aa,mset] = func_harmonic(krange,cx,cy,x_e,y_e,xx,fd)
    M0 = 0;
    mset = [];
    aa = [];
    for kk=1:length(krange)
        kkk = krange(kk);
        z_h = (x_e-cx)+1i*(y_e+cy);
        conjz_h = conj(z_h);
        z_h_p = z_h.^(-kkk);
        ssdi = diff(xx);
        ds = ssdi(1);
        dz = fd(:,1).*ds+1i.*fd(:,2).*ds;  
        m_kkk = -1*z_h_p.*conjz_h.*dz/(2*1i*pi);
        m_kkk = sum(abs(m_kkk));
        aa(end+1) = m_kkk;
        if kkk == 0
            M0 = m_kkk;
        else
            norm_m = m_kkk/M0^(2/(2-kkk));
            mset(end+1) = norm_m;
        end
    end

end