%--------------------------------------------------------------------------
% Iteration to find the optimal alphashape given the detected points
%
% Crafted by Jayda Robison
% Adapted by Steven Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------

function [prevarea,prevan,edges,bestalpha_r,ashape] ...
    = ashape_iteration(f,bw,bw3,prevarea,alphanum,areaopen,timeInterval,adjw,wst)
    %% create Alpha Shape
    [ptheta, py] = find(bw3==1);
    origp = [py,-ptheta];
    shp = alphaShape(origp);
    
    % specify properties of shape
    thisa = criticalAlpha(shp,'all-points'); %% key part
    shp2 = alphaShape(origp,thisa,'HoleThreshold',500000000);
    shp2.Alpha = alphanum; 
    areanow = area(shp2);
    
    %% edge case consideration
    % first base case
    if f == timeInterval(1)
        prevarea = areanow;
        prevan = 0;
    else
        % adjust bwareaopen if area gets bigger
        prevao = areaopen;

        for i=1:100
            if areanow > 1.1*prevarea
                areaopen = areaopen + 1;
                bw3 = bwareaopen(bw, areaopen);
                % stick together
                for ii=1:size(adjw,1)
                    for j=1:size(adjw, 2)
                        bw3(wst+ii,j) = adjw(ii,j);
                    end
                end

                % create Alpha Shape
                [ptheta, py] = find(bw3==1);
                p = [py,-ptheta];
                shp = alphaShape(p);
        
                % specify properties of shape
                thisa = criticalAlpha(shp, 'all-points');
                shp2 = alphaShape(p,thisa,'HoleThreshold',500000000);
                shp2.Alpha = alphanum;
                areanow = area(shp2);
            end
        end
    
        areaopen = prevao;
    
        % adjust alpha number if area gets smaller
        prevan = alphanum;
        
        for i=1:100
            if areanow < 0.7*prevarea
                alphanum = alphanum + 5;
                bw3 = bwareaopen(bw, areaopen);
                % stick together
                for ii=1:size(adjw,1)
                    for j=1:size(adjw, 2)
                        bw3(wst+ii,j) = adjw(ii,j);
                    end
                end

                % create Alpha Shape
                [ptheta, py] = find(bw3==1);
                p = [py,-ptheta]; 
                shp = alphaShape(p);
        
                % specify properties of shape
                thisa = criticalAlpha(shp, 'all-points');
                shp2 = alphaShape(p,thisa,'HoleThreshold',500000000);
                shp2.Alpha = alphanum;
                areanow = area(shp2);
            end
        end
    
        alphanum = prevan;
        prevarea = areanow;
    end
    
    % find and plot edge of alpha shape
    [~,edges] = boundaryFacets(shp2);
    bestalpha_r = shp2.Alpha; % best alpha radius
    ashape = shp2; % store the best alpha shape
end


