function  omg = AdjustIterations(Uit,Tit,omg,r)

while 1
    Uit(:,r) = omg*Uit(:,r) + (1-omg)*Uit(:,r-1);
    Tit(:,r) = omg*Tit(:,r) + (1-omg)*Tit(:,r-1);
    
    % calculate error at this step
    errvn = norm(Uit(:,r)-Uit(:,r-1));
    errtn = norm(Tit(:,r)-Tit(:,r-1));
    
    % calculate error at previous step
    errvp = norm(Uit(:,r-1)-Uit(:,r-2));
    errtp = norm(Tit(:,r-1)-Tit(:,r-2));
    
    if (errvn > errvp || errtn > errtp)
        if omg >= 0.1
            omg = omg-0.05; 
            fprintf(['Omega changed to ',num2str(omg),'\n']);
        else
            break;
        end
    else
        break;
    end
    
end

end