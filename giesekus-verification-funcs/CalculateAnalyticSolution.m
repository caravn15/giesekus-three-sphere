function uAn = CalculateAnalyticSolution(params,L,y,type)

pin   = params.pin;
pout  = params.pout;
alpha = params.alpha;
De    = params.De;
G     = (pin-pout)/L;

switch type
    case 'Newtonian'

        uAn = -G/2*y.*(1-y)+y;

    case 'Giesekus'
        if De==0
            uAn = -G/2*y.*(1-y)+y;
        else

            options = optimset('Display','off');

            t00v = linspace(-1,0,500);
            for ii=1:length(t00v)

                % initial condition
                t00 = t00v(ii);

                % bc function
                bcfunc = @(t) (2*alpha-1)/(2*alpha*De^2*G)*log((1-alpha*De^2*t^2)/...
                    (1-alpha*De^2*(t+G)^2))-((alpha-1)*(G+2*t))/((1-alpha*De^2*t^2)*...
                    (1-alpha*De^2*(t+G)^2))-1;

                % solve for t0
                t0 = fsolve(bcfunc,t00,options);

                % velocity function
                vel = @(y) -1/(2*alpha*G*De^2)*((2*(alpha-1))./...
                    (1-alpha*De^2*(t0+G*y).^2)+(2*alpha-1)*...
                    log(1-alpha*De^2*(t0+G*y).^2));

                % integration constant C
                C = -vel(0);

                % velocity profile
                ue = vel(y) + C;

                % checks
                if isreal(ue(end)) && abs(ue(end)-1) < 0.1
                    solinds(ii) = ii;
                else
                    solinds(ii) = 0;
                end

                ueEnd(ii) = ue(end);
            end

            isnz       = solinds == 0;
            evals      = abs(ueEnd(:)-1)+1000*isnz';
            [~,indmin] = min(evals);


            t00v = linspace(t00v(indmin)-0.1,t00v(indmin)+0.1,500);
            for ii=1:length(t00v)

                % initial condition
                t00 = t00v(ii);

                % bc function
                bcfunc = @(t) (2*alpha-1)/(2*alpha*De^2*G)*log((1-alpha*De^2*t^2)/...
                    (1-alpha*De^2*(t+G)^2))-((alpha-1)*(G+2*t))/((1-alpha*De^2*t^2)*...
                    (1-alpha*De^2*(t+G)^2))-1;

                % solve for t0
                t0 = fsolve(bcfunc,t00,options);

                % velocity function
                vel = @(y) -1/(2*alpha*G*De^2)*((2*(alpha-1))./...
                    (1-alpha*De^2*(t0+G*y).^2)+(2*alpha-1)*...
                    log(1-alpha*De^2*(t0+G*y).^2));

                % integration constant C
                C = -vel(0);

                % velocity profile
                ue = vel(y) + C;

                % checks
                if isreal(ue(end)) == 1 && abs(ue(end)-1) < 0.01
                    solinds(ii) = ii;
                else
                    solinds(ii) = 0;
                end

                ueEnd(ii) = ue(end);
            end

            isnz       = solinds == 0;
            evals      = abs(ueEnd(:)-1)+1000*isnz';
            [~,indmin] = min(evals);

            t00 = t00v(indmin);

            % bc function
            bcfunc = @(t) (2*alpha-1)/(2*alpha*De^2*G)*log((1-alpha*De^2*t^2)/...
                (1-alpha*De^2*(t+G)^2))-((alpha-1)*(G+2*t))/((1-alpha*De^2*t^2)*...
                (1-alpha*De^2*(t+G)^2))-1;

            % solve for t0
            t0 = fsolve(bcfunc,t00,options);

            % velocity function
            vel = @(y) -1/(2*alpha*G*De^2)*((2*(alpha-1))./...
                (1-alpha*De^2*(t0+G*y).^2)+(2*alpha-1)*...
                log(1-alpha*De^2*(t0+G*y).^2));

            % integration constant C
            C = -vel(0);

            % velocity profile
            uAn = vel(y) + C;
        end

end

end