function [msd,mse,emse,w_n] = fir_lms(u,d,yo,wo,M,N,mu,normalize)
    % Preallocate variables to improve performance
    w_n = zeros(N,M);
    mse = zeros(N,1);
    emse = zeros(N,1);
    msd = zeros(N,1);
    w = zeros(M,1);        
    u_reg = randn(1,M);

    % Assert dimensions to make sure arguments have the right ones.
    assertDimensions(w,wo);
    assertDimensions(u,zeros(N,1));
    assertDimensions(d,zeros(N,1));
    assertDimensions(yo,zeros(N,1));
    assertDimensions(wo,zeros(M,1));
    assertDimensions(M,1);
    assertDimensions(N,1);
    assert(islogical(normalize));
    
    % Define if it will run LMS or NMLS
    if normalize
       normalize_function =  @(u_reg) (1e-6 + norm(u_reg)^2);    
    else
       normalize_function =  @(u_reg) (1);    
    end
    
    % iterate over all samples
    for n=1:N
        % store current filter estimate
        w_n(n,:) = w';
        
        % Update  regressor
        u_reg = [u(n) u_reg(1:M-1)];
        
        % Current output        
        y = u_reg * w;    
        
        % Estimation Error
        e = d(n) - y;   
        
        % A priori error
        ea = yo(n) - y;             
                
        % Mean Square error
        mse(n) = e^2;             
        
        % Excess Mean Square error
        emse(n) = ea^2;        
        
        % Mean Square deviation
        msd(n) = norm(wo - w)^2;    
        
        % Update the coefficients
        w = w + mu * e * u_reg' / normalize_function(u_reg);    
        
    end
end
