function epp = RunIterations(uit,tit,pit,r,epp)

while 1
    uit(:,r) = epp*uit(:,r)+(1-epp)*uit(:,r-1);
    pit(:,r) = epp*pit(:,r)+(1-epp)*pit(:,r-1);
    tit(:,r) = epp*tit(:,r)+(1-epp)*tit(:,r-1);
    
    % calculate error at this step
    errvn = norm(uit(:,r)-uit(:,r-1));
    errtn = norm(tit(:,r)-tit(:,r-1));
    
    % calculate error at previous step
    errvp = norm(uit(:,r-1)-uit(:,r-2));
    errtp = norm(tit(:,r-1)-tit(:,r-2));
   
    if (errvn > errvp || errtn > errtp)
        if round(epp,1) <= 0.1
            epp=epp/2;
            fprintf(['Epsilon changed to ',num2str(epp),'\n']);
            break;
        else
 epp = epp-0.1;
 
        fprintf(['Epsilon changed to ',num2str(epp),'\n']);
        end
    else
%         if (r==50 || r==100 || r==150 || r==200)
%            epp = 1;
%            fprintf(['Epsilon changed to ',num2str(epp),'\n']);
%         end
    break;
    end
    
    
    
end


end