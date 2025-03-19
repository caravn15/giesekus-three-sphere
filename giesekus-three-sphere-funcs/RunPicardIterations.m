function [U,P,T] = RunPicardIterations(model,U0,P0,T0,FixedTime,FixedIts)

%% setup

% unpack model structure
tol = model.num.tol;

% initilise
Uit(:,1) = U0;
Pit(:,1) = P0;
Tit(:,1) = T0;

%% iterate

slowsp = 0; FinalWi = model.Wi; itconverged = 0;


while true

    if slowsp ~= 0
        Wivec = linspace(0.25,FinalWi,slowsp);
        fprintf(['Trying again and stepping with ',num2str(slowsp),...
            ' steps...\n']);
    end

    r=2; solgrowing = 0;
    while true
        if slowsp ~= 0
            try
                model.Wi = Wivec(r);
            catch
                model.Wi = Wivec(end);
            end
        end

        % solve non-linear (but linearised) probem
        [Uit(:,r),Pit(:,r),Tit(:,r)] = SolveNonLinearProblem...
            (model,Tit(:,r-1),FixedTime,FixedIts);

        % calculate errors (this step)
        Uerr = norm(Uit(:,r)-Uit(:,r-1));
        Terr = norm(Tit(:,r)-Tit(:,r-1));

        if r>slowsp+10

            % calculate error at previous step
            errvp = norm(Uit(:,r-1)-Uit(:,r-2));
            errtp = norm(Tit(:,r-1)-Tit(:,r-2));

            if (Uerr > errvp || Terr > errtp)
                solgrowing = 1;
            end

        end

        % print errors
        fprintf(['Iteration ',num2str(r-1),' Uerr = ',num2str(Uerr), ...
            ', Terr = ',num2str(Terr),'\n']);

        if r>slowsp+1 && Uerr < tol && Terr < tol
            itconverged = 1;
            break;
        end

        if r > slowsp+100 || Terr > 1000 || solgrowing == 1
            slowsp = slowsp + 10;
            % error('Something is not working...')
            break;
        end

        if slowsp > 100
            error('Something is not working...')
        end
        r = r + 1;

    end

    if itconverged == 1
        break;
    end

end

%% sort outputs


% get final solution
U = Uit(:,end);
P = Pit(:,end);
T = Tit(:,end);

fprintf(['Picard iterations converged after ',...
    num2str(r-1),' iterations.\n']);

end