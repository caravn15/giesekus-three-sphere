function PrintParameters(model)

fprintf('###############################################\n');
fprintf('#  RUNNING GIESEKUS THREE SPHERE SWIMMER SIM  #\n');
fprintf('-----------------------------------------------\n');
fprintf(['   PARAMETERS: alpha  = ',num2str(model.alpha),'\n']);
fprintf(['               Wi     = ',num2str(model.Wi),'\n']);
fprintf(['               kappa  = ',num2str(model.kappa),'\n']);
fprintf(['               Amp    = ',num2str(model.Amp),'\n']);
fprintf(['               chi    = ',num2str(model.chi),'\n']);
fprintf(['               tmin   = ',num2str(model.tmin),'\n']);
fprintf(['               tmax   = ',num2str(model.tmax),'\n']);
fprintf('-----------------------------------------------\n');

end