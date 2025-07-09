function [mua_RMSE_O, musp_RMSE_O, mua0, musp0] = ...
    posOrbitRMSE(optChars, ...
    rsPen_ass, rd_ass, rOrb, optProp_act, dmua_act, omega, ...
    recMeth)
arguments
    optChars (:,1) cell;

    rsPen_ass (2,3) double;
    rd_ass (2,3) double;

    rOrb (1,1) double;
    
    optProp_act struct;
    % Contains: mua, musp, nin, nout
    
    dmua_act (1,1) double = 0.0001; %1/mm
    omega (1,1) double = 100e6*2*pi; %rad/s

    recMeth string = "iter"; %"iter" or "lin"
end
    
    rd_act=rd_ass;
    rsPen_act=rsPen_ass;
    
    [~, ~, ~, ~, optProp_rec]=DSpos2err( ...
        rsPen_act, rd_act, rsPen_ass, rd_ass, ...
        optProp_act, dmua_act, omega, recMeth);
    mua0=optProp_rec.mua;
    musp0=optProp_rec.musp;
    
    switch length(optChars)
        case 1
            theta=linspace(0, 2*pi, 1001);
            theta(end)=[];
            dx=rOrb*cos(theta);
            dy=rOrb*sin(theta);
            
            mua=NaN(length(theta), 1);
            musp=NaN(length(theta), 1);
            for i=1:length(theta)
                optChar=optChars{1};
                ind=str2double(optChar(2));
                switch optChar(1)
                    case 'S'
                        rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(i), dy(i), 0];
                    case 'D'
                        rd_act(ind, :)=rd_ass(ind, :)+[dx(i), dy(i), 0];
                    otherwise
                        error('optChar(1) should be ''S'' or ''D''')
                end
                
                [~, ~, ~, ~, optProp_rec]=DSpos2err( ...
                    rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                    optProp_act, dmua_act, omega, recMeth);
                mua(i)=optProp_rec.mua;
                musp(i)=optProp_rec.musp;
            end

        case 2
            theta=linspace(0, 2*pi, 101);
            theta(end)=[];
            dx=rOrb*cos(theta);
            dy=rOrb*sin(theta);
            
            mua=NaN(length(theta), length(theta));
            musp=NaN(length(theta), length(theta));
            for i=1:length(theta)
                optChar=optChars{1};
                ind=str2double(optChar(2));
                switch optChar(1)
                    case 'S'
                        rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(i), dy(i), 0];
                    case 'D'
                        rd_act(ind, :)=rd_ass(ind, :)+[dx(i), dy(i), 0];
                    otherwise
                        error('optChar(1) should be ''S'' or ''D''')
                end
                
                for j=1:length(theta)
                    optChar=optChars{2};
                    ind=str2double(optChar(2));
                    switch optChar(1)
                        case 'S'
                            rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(j), dy(j), 0];
                        case 'D'
                            rd_act(ind, :)=rd_ass(ind, :)+[dx(j), dy(j), 0];
                        otherwise
                            error('optChar(1) should be ''S'' or ''D''')
                    end
                
                    [~, ~, ~, ~, optProp_rec]=DSpos2err( ...
                        rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                        optProp_act, dmua_act, omega, recMeth);
                    mua(i, j)=optProp_rec.mua;
                    musp(i, j)=optProp_rec.musp;
                end
            end
            
        case 3
            theta=linspace(0, 2*pi, 21);
            theta(end)=[];
            dx=rOrb*cos(theta);
            dy=rOrb*sin(theta);
            
            mua=NaN(length(theta), length(theta), length(theta));
            musp=NaN(length(theta), length(theta), length(theta));
            for i=1:length(theta)
                optChar=optChars{1};
                ind=str2double(optChar(2));
                switch optChar(1)
                    case 'S'
                        rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(i), dy(i), 0];
                    case 'D'
                        rd_act(ind, :)=rd_ass(ind, :)+[dx(i), dy(i), 0];
                    otherwise
                        error('optChar(1) should be ''S'' or ''D''')
                end
                
                for j=1:length(theta)
                    optChar=optChars{2};
                    ind=str2double(optChar(2));
                    switch optChar(1)
                        case 'S'
                            rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(j), dy(j), 0];
                        case 'D'
                            rd_act(ind, :)=rd_ass(ind, :)+[dx(j), dy(j), 0];
                        otherwise
                            error('optChar(1) should be ''S'' or ''D''')
                    end

                    for k=1:length(theta)
                        optChar=optChars{3};
                        ind=str2double(optChar(2));
                        switch optChar(1)
                            case 'S'
                                rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(k), dy(k), 0];
                            case 'D'
                                rd_act(ind, :)=rd_ass(ind, :)+[dx(k), dy(k), 0];
                            otherwise
                                error('optChar(1) should be ''S'' or ''D''')
                        end
                    
                        [~, ~, ~, ~, optProp_rec]=DSpos2err( ...
                            rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                            optProp_act, dmua_act, omega, recMeth);
                        mua(i, j, k)=optProp_rec.mua;
                        musp(i, j, k)=optProp_rec.musp;
                    end
                end
            end

        case 4
            theta=linspace(0, 2*pi, 11);
            theta(end)=[];
            dx=rOrb*cos(theta);
            dy=rOrb*sin(theta);
            
            mua=NaN(length(theta), length(theta), length(theta), length(theta));
            musp=NaN(length(theta), length(theta), length(theta), length(theta));
            for i=1:length(theta)
                optChar=optChars{1};
                ind=str2double(optChar(2));
                switch optChar(1)
                    case 'S'
                        rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(i), dy(i), 0];
                    case 'D'
                        rd_act(ind, :)=rd_ass(ind, :)+[dx(i), dy(i), 0];
                    otherwise
                        error('optChar(1) should be ''S'' or ''D''')
                end
                
                for j=1:length(theta)
                    optChar=optChars{2};
                    ind=str2double(optChar(2));
                    switch optChar(1)
                        case 'S'
                            rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(j), dy(j), 0];
                        case 'D'
                            rd_act(ind, :)=rd_ass(ind, :)+[dx(j), dy(j), 0];
                        otherwise
                            error('optChar(1) should be ''S'' or ''D''')
                    end

                    for k=1:length(theta)
                        optChar=optChars{3};
                        ind=str2double(optChar(2));
                        switch optChar(1)
                            case 'S'
                                rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(k), dy(k), 0];
                            case 'D'
                                rd_act(ind, :)=rd_ass(ind, :)+[dx(k), dy(k), 0];
                            otherwise
                                error('optChar(1) should be ''S'' or ''D''')
                        end

                        for l=1:length(theta)
                            optChar=optChars{4};
                            ind=str2double(optChar(2));
                            switch optChar(1)
                                case 'S'
                                    rsPen_act(ind, :)=rsPen_ass(ind, :)+[dx(l), dy(l), 0];
                                case 'D'
                                    rd_act(ind, :)=rd_ass(ind, :)+[dx(l), dy(l), 0];
                                otherwise
                                    error('optChar(1) should be ''S'' or ''D''')
                            end
                        
                            [~, ~, ~, ~, optProp_rec]=DSpos2err( ...
                                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                                optProp_act, dmua_act, omega, recMeth);
                            mua(i, j, k, l)=optProp_rec.mua;
                            musp(i, j, k, l)=optProp_rec.musp;
                        end
                    end
                end
            end
            
        otherwise
            error('Must be at most four optodes');
    end
    
    mua_RMSE_O=sqrt( ...
        sum((mua(:)-mua0).^2) ...
        /numel(mua));
    musp_RMSE_O=sqrt( ...
        sum((musp(:)-musp0).^2) ...
        /numel(musp));

end