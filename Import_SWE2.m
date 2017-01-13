%% Generates all swe otpions

OPTIONS

for t = 2:9
run OPTIONS.m
options.DensitySWE  = t;
run MAIN

  for i = 1:3
    glacier = char(options.glacier(i)); 
    sweOPT(t).(glacier) = [SWE(i).swe, SWE(i).utm(:,1:2)];
  end
end

    clear i t glacier